//! The search driver: spectrum → candidate peptides → RawScore → SpecEValue.
//!
//! The generating function depends only on `(spectrum, precursor mass, isotope range, amino-acid
//! alphabet)` — never on a candidate peptide — so it is built **once per `(spectrum, charge)`** and
//! every candidate in the precursor window becomes a RawScore plus a tail lookup. That is the whole
//! reason a generating-function search is affordable, and it is why candidate generation (this
//! module) and the DP (`msgf-genfunc`) stay separate.
//!
//! ## Deliberate divergences from MS-GF+ (see `CLAUDE.md`, "Fidelity is the contract")
//!
//! - **Cleavage scoring is modelled for C-terminal enzymes only.** The de novo graph is built in
//!   the reverse (C-terminal) direction, which is where MS-GF+ puts the peptide-cleavage credit for
//!   trypsin-like enzymes. For an N-terminal enzyme (Lys-N, Asp-N) MS-GF+ builds the graph in the
//!   opposite direction; rather than silently applying the wrong credit we disable cleavage scoring
//!   for those enzymes and for unspecific ones. [`SearchEngine::warnings`] reports it.
//! - **`EValue = SpecEValue × database size`** (the number of candidates in the index), or the
//!   explicit [`SearchParams::db_size`] when given. MS-GF+ derives its own candidate-count estimate
//!   internally, so E-values are the same order of magnitude but not directly comparable.
//!   **Q-values are computed from SpecEValue**, so this scaling does not affect FDR at all.

use msgf_chem::Tolerance;
use msgf_db::enzyme::DigestParams;
use msgf_db::fasta::ProteinDb;
use msgf_genfunc::{score_distribution, AlphabetEntry, NullModel};
use msgf_io::Spectrum;
use msgf_scorer::{
    match_and_terminal_score, prepare_with_cache, Candidate as ScoredCandidate, Cleavage, Ms2,
    PreparedSpectrum, Residue, ScoreModel, RESIDUE_ORDER,
};
use rayon::prelude::*;
use std::collections::HashMap;

use crate::index::{Candidate, PeptideIndex};
use crate::mods::{ModPosition, ModSet};

/// Mass difference between the ¹³C and ¹²C isotopes — one isotope-error step on the precursor.
pub const ISOTOPE_STEP: f64 = 1.003_354_838;

/// Cleavage credit and penalty applied at an enzymatic / non-enzymatic terminus (the MS-GF+
/// convention for trypsin).
pub const CLEAVAGE_CREDIT: i32 = 2;
pub const CLEAVAGE_PENALTY: i32 = -11;

/// Tunables for one search.
#[derive(Debug, Clone)]
pub struct SearchParams {
    /// Precursor mass tolerance.
    pub precursor_tol: Tolerance,
    /// Isotope-error range `(lo, hi)`, like MS-GF+ `-ti`. `(0, 1)` allows the precursor to have
    /// been picked one ¹³C isotope high.
    pub isotope_errors: (i32, i32),
    /// How many matches to report per spectrum (MS-GF+ `-n`).
    pub num_matches: usize,
    /// Charges to try when the spectrum does not declare one.
    pub charge_range: (i32, i32),
    /// Override for the E-value database-size multiplier. `None` = the candidate-index size.
    pub db_size: Option<f64>,
}

impl Default for SearchParams {
    fn default() -> Self {
        Self {
            precursor_tol: Tolerance::ppm(10.0),
            isotope_errors: (0, 1),
            num_matches: 1,
            charge_range: (2, 3),
            db_size: None,
        }
    }
}

/// One peptide-spectrum match.
#[derive(Debug, Clone)]
pub struct Psm {
    pub spec_index: usize,
    pub scan: String,
    pub title: String,
    pub precursor_mz: f64,
    pub charge: i32,
    pub isotope_error: i32,
    /// `(observed − theoretical) / theoretical × 1e6`, after removing the isotope error.
    pub precursor_error_ppm: f64,
    /// Peptide with flanking context and inline mod deltas, e.g. `K.SAM+15.995PLER.A`.
    pub peptide: String,
    /// Mod-bearing sequence with flanks stripped and upper-cased — the peptide identity used for
    /// peptide-level FDR (`plans/PLAN2.md` §1.4).
    pub peptide_key: String,
    /// Every protein occurrence of this peptide, in database order.
    pub proteins: Vec<String>,
    /// `true` only when **every** occurrence is a decoy (`plans/PLAN2.md` §1.3).
    pub is_decoy: bool,
    pub raw_score: i32,
    pub denovo_score: i32,
    pub spec_evalue: f64,
    pub evalue: f64,
    /// PSM-level q-value, filled in by [`crate::assign_q_values`].
    pub q_value: f32,
    /// Peptide-level q-value, filled in by [`crate::assign_q_values`].
    pub pep_q_value: f32,
}

/// How the enzyme's specificity maps onto cleavage scoring.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CleavageMode {
    /// C-terminal enzyme (trypsin-like): full credit/penalty at both termini.
    CTerminal,
    /// N-terminal or unspecific enzyme: cleavage scoring off (see the module docs).
    Off,
}

/// Everything a search needs, assembled once and shared across spectra.
pub struct SearchEngine<'a> {
    model: &'a ScoreModel,
    db: &'a ProteinDb,
    index: &'a PeptideIndex,
    mods: &'a ModSet,
    params: SearchParams,
    /// The generating function's null model: the 20 residues (with fixed mods folded in) plus one
    /// entry per variable-mod variant, each weighted by its database frequency; the enzyme's
    /// terminus rule; the isotope-error range.
    null: NullModel,
    /// Residues the enzyme cleaves at, for terminal cleavage scoring.
    cleave_at: Vec<u8>,
    warnings: Vec<String>,
}

/// Per-thread reusable buffers for scoring one spectrum's candidates.
#[derive(Default)]
pub struct SearchScratch {
    peaks: Vec<(f64, f64)>,
    buf: ScoreBuffers,
}

impl<'a> SearchEngine<'a> {
    /// Assemble an engine. The amino-acid background frequencies are taken from `db`, which is what
    /// makes the SpecEValue reflect the composition of the database actually being searched.
    pub fn new(
        model: &'a ScoreModel,
        db: &'a ProteinDb,
        index: &'a PeptideIndex,
        mods: &'a ModSet,
        digest_params: &DigestParams,
        params: SearchParams,
    ) -> SearchEngine<'a> {
        let mut warnings = Vec::new();
        let probs: HashMap<u8, f64> = db.aa_probabilities().into_iter().collect();
        let prob_of = |r: u8| probs.get(&r).copied().unwrap_or(0.0);

        // Base alphabet: each standard residue at its fixed-modified mass.
        let mut alphabet: Vec<AlphabetEntry> = RESIDUE_ORDER
            .iter()
            .map(|&residue| AlphabetEntry {
                letter: residue,
                mass: msgf_chem::residue_mass(residue).expect("standard residue")
                    + mods.fixed_residue_delta(residue),
                prob: prob_of(residue),
            })
            .collect();

        // One extra edge per variable-mod variant, at the same background frequency as the
        // unmodified residue — the convention `msgf-cli`'s `--ox-m` was validated with.
        let base_residues: Vec<u8> = alphabet.iter().map(|a| a.letter).collect();
        for (_, spec) in mods.variable() {
            if spec.position != ModPosition::Any {
                warnings.push(format!(
                    "variable mod `{}` is position-restricted ({:?}); it is searched but not represented \
                     in the de novo graph alphabet, so its SpecEValue is slightly conservative",
                    spec.name, spec.position
                ));
                continue;
            }
            let targets: &[u8] = if spec.residues.is_empty() {
                &base_residues
            } else {
                &spec.residues
            };
            for &r in targets {
                let Some(base) = msgf_chem::residue_mass(r) else {
                    continue;
                };
                let m = base + mods.fixed_residue_delta(r) + spec.mass;
                if m <= 0.0 {
                    continue;
                }
                alphabet.push(AlphabetEntry {
                    letter: r,
                    mass: m,
                    prob: prob_of(r),
                });
            }
        }

        let enzyme = &digest_params.enzyme;
        let cleavage_mode = if enzyme.is_unspecific() {
            CleavageMode::Off
        } else if enzyme.c_term {
            CleavageMode::CTerminal
        } else {
            warnings.push(format!(
                "enzyme `{}` cleaves N-terminal to its residues; the de novo graph is built in the \
                 C-terminal direction, so cleavage credit/penalty is disabled for this search \
                 (scores and SpecEValue are therefore not comparable to MS-GF+ for this enzyme)",
                enzyme.name
            ));
            CleavageMode::Off
        };
        let cleave_at = enzyme.cleave_at.clone();
        // The terminus rule of both the RawScore and the null distribution. The neighbouring
        // residue's cleavage is not known a priori, so the null weights it by the summed database
        // frequency of the cleavage residues.
        let cleavage = Cleavage {
            sites: cleave_at.clone(),
            credit: CLEAVAGE_CREDIT,
            penalty: CLEAVAGE_PENALTY,
            enabled: cleavage_mode == CleavageMode::CTerminal,
        };
        let null = NullModel {
            alphabet,
            cleavage,
            isotope: params.isotope_errors,
        };

        SearchEngine {
            model,
            db,
            index,
            mods,
            params,
            null,
            cleave_at,
            warnings,
        }
    }

    /// Non-fatal configuration notes (unsupported combinations that were degraded, not rejected).
    pub fn warnings(&self) -> &[String] {
        &self.warnings
    }

    /// The multiplier turning a SpecEValue into an E-value.
    fn db_size(&self) -> f64 {
        self.params.db_size.unwrap_or(self.index.len() as f64)
    }

    /// Search every spectrum, in parallel. Results are ordered by `(spec_index, spec_evalue)`.
    /// Q-values are **not** filled in — FDR is global, so run [`crate::assign_q_values`] on the
    /// whole result set afterwards.
    pub fn run(&self, spectra: &[Spectrum]) -> Vec<Psm> {
        let mut out: Vec<Psm> = spectra
            .par_iter()
            .enumerate()
            .map_init(SearchScratch::default, |scratch, (i, spec)| {
                self.search_spectrum(scratch, i, spec)
            })
            .flatten()
            .collect();
        out.sort_by(|a, b| {
            a.spec_index.cmp(&b.spec_index).then(
                a.spec_evalue
                    .partial_cmp(&b.spec_evalue)
                    .expect("finite e-values"),
            )
        });
        out
    }

    /// Search one spectrum, trying every plausible charge and keeping the best matches overall.
    pub fn search_spectrum(
        &self,
        scratch: &mut SearchScratch,
        spec_index: usize,
        spec: &Spectrum,
    ) -> Vec<Psm> {
        let Some(mz) = spec.precursor_mz else {
            return Vec::new();
        };
        let charges: Vec<i32> = match spec.charge {
            Some(z) if z > 0 => vec![z],
            _ => (self.params.charge_range.0..=self.params.charge_range.1).collect(),
        };
        let mut best: Vec<Psm> = Vec::new();
        for z in charges {
            best.extend(self.search_at_charge(scratch, spec_index, spec, mz, z));
        }
        best.sort_by(|a, b| {
            a.spec_evalue
                .partial_cmp(&b.spec_evalue)
                .expect("finite e-values")
                .then(b.raw_score.cmp(&a.raw_score))
        });
        best.truncate(self.params.num_matches);
        best
    }

    fn search_at_charge(
        &self,
        scratch: &mut SearchScratch,
        spec_index: usize,
        spec: &Spectrum,
        mz: f64,
        charge: i32,
    ) -> Vec<Psm> {
        // Neutral precursor mass (= candidate peptide mass, water included).
        let parent_mass = mz as f32 * charge as f32 - charge as f32 * msgf_scorer::PROTON as f32;
        let (ti_lo, ti_hi) = self.params.isotope_errors;

        // --- the per-spectrum half: prepare the peaks (None = implausible precursor mass) ---
        scratch.peaks.clear();
        scratch
            .peaks
            .extend(spec.peaks.iter().map(|p| (p.mz, p.intensity)));
        let ms2 = Ms2 {
            peaks: &scratch.peaks,
            precursor_mz: mz,
            charge,
        };
        let Some(prep) = prepare_with_cache(self.model, &ms2, (-ti_lo).max(2)) else {
            return Vec::new();
        };

        // --- the per-candidate half: every peptide in the precursor window gets a RawScore ---
        // Identical peptides occurring in several proteins score identically, so they are grouped
        // into one match carrying every protein occurrence — that is what decides decoy status
        // (`plans/PLAN2.md` §1.3) and stops a repeated peptide consuming the whole top-N list.
        //
        // This runs *before* the generating function, which needs nothing from it but gains a great
        // deal: the RawScore of the worst PSM we will report is the tail threshold, and the DP can
        // then skip every score cell that provably cannot reach it.
        let mut grouped: HashMap<String, Hit> = HashMap::new();
        let buf = &mut scratch.buf;
        for k in ti_lo..=ti_hi {
            let target = parent_mass as f64 - k as f64 * ISOTOPE_STEP;
            let win = self.params.precursor_tol.window_da(target);
            for cand in self.index.window(target - win, target + win) {
                let key = self.peptide_string(cand, false);
                match grouped.get_mut(&key) {
                    Some(hit) => hit.proteins.push(cand.protein),
                    None => {
                        let raw_score = self.raw_score(&prep, cand, buf);
                        grouped.insert(
                            key,
                            Hit {
                                raw_score,
                                isotope_error: k,
                                candidate: *cand,
                                proteins: vec![cand.protein],
                            },
                        );
                    }
                }
            }
        }
        if grouped.is_empty() {
            return Vec::new(); // no candidates — the whole generating function is skipped
        }
        let mut hits: Vec<(String, Hit)> = grouped.into_iter().collect();
        // Highest RawScore first; the SpecEValue tail is monotone in RawScore, so for a single
        // spectrum this is also best-SpecEValue order. The peptide key breaks ties deterministically.
        hits.sort_by(|a, b| b.1.raw_score.cmp(&a.1.raw_score).then(a.0.cmp(&b.0)));
        hits.truncate(self.params.num_matches);
        // Every reported PSM is looked up at or above this score, so nothing below it is needed.
        let threshold = hits.iter().map(|h| h.1.raw_score).min().unwrap_or(i32::MIN);

        // --- the generating function: built once over the isotope-error sinks, tail-pruned to the
        // reported PSMs' scores. An isotope error of +k means the measured precursor is ~k Da high,
        // so the true peptide mass is k nominal bins lower.
        let Some(tail) = score_distribution(&prep, &self.null, Some(threshold)) else {
            return Vec::new();
        };
        let denovo = tail.best_possible();
        let _ = ti_hi;

        let db_size = self.db_size();
        hits.into_iter()
            .map(|(peptide_key, hit)| {
                // An exact-zero SpecEValue (RawScore above the support) is carried as -0.0, the
                // established output of this tool (`-0.000000e0`).
                let spec_evalue = match tail.tail_mass(hit.raw_score) {
                    v if v == 0.0 => -0.0,
                    v => v,
                };
                let cand = &hit.candidate;
                let observed = parent_mass as f64 - hit.isotope_error as f64 * ISOTOPE_STEP;
                let proteins: Vec<String> = hit
                    .proteins
                    .iter()
                    .map(|&p| self.db.proteins[p as usize].name.clone())
                    .collect();
                let is_decoy = hit
                    .proteins
                    .iter()
                    .all(|&p| self.db.proteins[p as usize].is_decoy);
                Psm {
                    spec_index,
                    scan: spec
                        .scan
                        .clone()
                        .unwrap_or_else(|| (spec_index + 1).to_string()),
                    title: spec.title.clone().unwrap_or_default(),
                    precursor_mz: mz,
                    charge,
                    isotope_error: hit.isotope_error,
                    precursor_error_ppm: (observed - cand.mass) / cand.mass * 1e6,
                    peptide: self.peptide_string(cand, true),
                    peptide_key: peptide_key.to_ascii_uppercase(),
                    proteins,
                    is_decoy,
                    raw_score: hit.raw_score,
                    denovo_score: denovo,
                    spec_evalue,
                    evalue: spec_evalue * db_size,
                    q_value: f32::NAN,
                    pep_q_value: f32::NAN,
                }
            })
            .collect()
    }

    /// RawScore for one candidate: the node + edge match score plus the terminal cleavage
    /// credit/penalty at both termini.
    fn raw_score(&self, prep: &PreparedSpectrum, cand: &Candidate, buf: &mut ScoreBuffers) -> i32 {
        buf.fill(cand, self.db, self.mods);
        let scored = ScoredCandidate {
            residues: &buf.residues,
            n_term_credit: self.n_term_credit(cand),
        };
        match_and_terminal_score(prep, &scored, &self.null.cleavage)
            .expect("candidates hold standard residues")
    }

    /// Whether the candidate's N-terminus counts as enzymatic, from its real protein context.
    ///
    /// Credit when the preceding residue is a cleavage residue, when the peptide starts the
    /// protein (there is no preceding residue to fail the test), or when it starts just after an
    /// excised initiator methionine. The C-terminal term (credit only when the last residue is a
    /// cleavage residue; ending the protein does **not** substitute) is applied by the scorer.
    /// Irrelevant when cleavage scoring is off.
    fn n_term_credit(&self, cand: &Candidate) -> bool {
        let protein = &self.db.proteins[cand.protein as usize];
        let start = cand.start as usize;
        let at_prot_n = start == protein.start;
        let after_initiator_met = start == protein.start + 1 && self.db.seq[protein.start] == b'M';
        at_prot_n || after_initiator_met || self.cleave_at.contains(&self.db.seq[start - 1])
    }

    /// Format a candidate as a peptide string: `K.SAM+15.995PLER.A` (with flanking protein context)
    /// or `SAM+15.995PLER` (without). Mod deltas use three decimals, matching MS-GF+'s TSV.
    pub fn peptide_string(&self, cand: &Candidate, with_context: bool) -> String {
        let protein = &self.db.proteins[cand.protein as usize];
        let (start, len) = (cand.start as usize, cand.len as usize);
        let seq = &self.db.seq[start..start + len];
        let at_prot_n = start == protein.start;
        let at_prot_c = start + len == protein.start + protein.len;

        let mut s = String::with_capacity(len + 16);
        if with_context {
            s.push(if at_prot_n {
                '-'
            } else {
                self.db.seq[start - 1] as char
            });
            s.push('.');
        }
        for (i, &r) in seq.iter().enumerate() {
            s.push(r as char);
            let delta = self.mods.fixed_delta(r, i, len, at_prot_n, at_prot_c)
                + cand.placement.delta_at(i, self.mods);
            if delta != 0.0 {
                s.push_str(&format!("{delta:+.3}"));
            }
        }
        if with_context {
            s.push('.');
            s.push(if at_prot_c {
                '-'
            } else {
                self.db.seq[start + len] as char
            });
        }
        s
    }
}

/// One grouped match while a spectrum is being searched.
struct Hit {
    raw_score: i32,
    isotope_error: i32,
    candidate: Candidate,
    proteins: Vec<u32>,
}

/// Reusable per-candidate buffer, so scoring millions of candidates does no repeated allocation.
#[derive(Default)]
struct ScoreBuffers {
    residues: Vec<Residue>,
}

impl ScoreBuffers {
    /// Materialise the candidate's residues with every fixed and variable modification delta.
    fn fill(&mut self, cand: &Candidate, db: &ProteinDb, mods: &ModSet) {
        let protein = &db.proteins[cand.protein as usize];
        let (start, len) = (cand.start as usize, cand.len as usize);
        let seq = &db.seq[start..start + len];
        let at_prot_n = start == protein.start;
        let at_prot_c = start + len == protein.start + protein.len;

        self.residues.clear();
        for (i, &r) in seq.iter().enumerate() {
            let delta = mods.fixed_delta(r, i, len, at_prot_n, at_prot_c)
                + cand.placement.delta_at(i, mods);
            self.residues.push(Residue { letter: r, delta });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mods::{ModPlacement, ModSet, ModSpec, PlacedMod};
    use msgf_db::fasta::Protein;

    #[test]
    fn isotope_step_is_the_c13_gap() {
        assert!((ISOTOPE_STEP - 1.00335).abs() < 1e-4);
    }

    #[test]
    fn score_buffers_match_the_string_parser() {
        let db = ProteinDb {
            seq: b"SAMPLERK".to_vec(),
            proteins: vec![Protein {
                name: "P".into(),
                desc: String::new(),
                start: 0,
                len: 8,
                is_decoy: false,
            }],
        };
        let mods = ModSet {
            mods: vec![ModSpec::parse("O1,M,opt,any,Oxidation").unwrap()],
            max_var_mods: 1,
        };
        let placement = ModPlacement {
            n: 1,
            slots: {
                let mut s = [PlacedMod::default(); crate::mods::MAX_PLACED_MODS];
                s[0] = PlacedMod { pos: 2, mod_idx: 0 };
                s
            },
        };
        let cand = Candidate {
            mass: 0.0,
            protein: 0,
            start: 0,
            len: 8,
            n_termini: 2,
            placement,
        };
        let mut buf = ScoreBuffers::default();
        buf.fill(&cand, &db, &mods);

        // The same peptide via the string parser.
        let delta = mods.mods[0].mass;
        let pep = format!("SAM{delta:+}PLERK");
        let residues = msgf_chem::peptide::parse(&pep).unwrap();
        let via_parser: Vec<Residue> = residues
            .iter()
            .map(|r| Residue {
                letter: r.aa,
                delta: r.mod_delta,
            })
            .collect();
        assert_eq!(buf.residues, via_parser);
        assert_eq!(
            msgf_scorer::cumulative_nominal(&buf.residues).unwrap(),
            msgf_chem::peptide::nominal_prefix_masses(&residues)
        );
        assert_eq!(buf.residues.len(), 8);
    }
}
