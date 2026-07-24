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

use msgf_chem::peptide::Residue;
use msgf_chem::{mass, scaling, Tolerance};
use msgf_db::enzyme::DigestParams;
use msgf_db::fasta::ProteinDb;
use msgf_genfunc::graph::{build_reverse_graph, Aa, PeptideCleavage};
use msgf_genfunc::{compute_into, merge_group, Cleavage, DpScratch, GenFunc};
use msgf_io::Spectrum;
use msgf_scorer::preprocess::preprocess;
use msgf_scorer::scored_spectrum::ScoredSpectrum;
use msgf_scorer::ScoringModel;
use rayon::prelude::*;
use std::collections::HashMap;

use crate::index::{Candidate, PeptideIndex};
use crate::mods::{ModPosition, ModSet};

/// Mass difference between the ¹³C and ¹²C isotopes — one isotope-error step on the precursor.
pub const ISOTOPE_STEP: f64 = 1.003_354_838;

/// Cleavage credit and penalty applied at an enzymatic / non-enzymatic terminus. These are the
/// values validated bit-exact against MS-GF+ for trypsin (see `msgf-cli`'s golden rescore test).
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
    /// peptide-level FDR (`PLAN2.md` §1.4).
    pub peptide_key: String,
    /// Every protein occurrence of this peptide, in database order.
    pub proteins: Vec<String>,
    /// `true` only when **every** occurrence is a decoy (`PLAN2.md` §1.3).
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
    model: &'a ScoringModel,
    db: &'a ProteinDb,
    index: &'a PeptideIndex,
    mods: &'a ModSet,
    params: SearchParams,
    /// Graph alphabet: the 20 residues (with fixed mods folded in) plus one entry per variable-mod
    /// variant, each weighted by its database frequency.
    alphabet: Vec<Aa>,
    /// Residues the enzyme cleaves at, for terminal cleavage scoring.
    cleave_at: Vec<u8>,
    cleavage_mode: CleavageMode,
    /// Summed database frequency of the cleavage residues — `probCleavageSites` in MS-GF+.
    prob_cleavage_sites: f64,
    warnings: Vec<String>,
}

impl<'a> SearchEngine<'a> {
    /// Assemble an engine. The amino-acid background frequencies are taken from `db`, which is what
    /// makes the SpecEValue reflect the composition of the database actually being searched.
    pub fn new(
        model: &'a ScoringModel,
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
        let mut alphabet: Vec<Aa> = msgf_genfunc::graph::standard_aa_nominal()
            .into_iter()
            .map(|(residue, _)| {
                let m = msgf_chem::residue_mass(residue).expect("standard residue")
                    + mods.fixed_residue_delta(residue);
                Aa {
                    residue,
                    nominal: scaling::nominal_bin(m as f32),
                    accurate_mass: m as f32,
                    prob: prob_of(residue),
                }
            })
            .collect();

        // One extra edge per variable-mod variant, at the same background frequency as the
        // unmodified residue — the convention `msgf-cli`'s `--ox-m` was validated with.
        let base_residues: Vec<u8> = alphabet.iter().map(|a| a.residue).collect();
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
                alphabet.push(Aa {
                    residue: r,
                    nominal: scaling::nominal_bin(m as f32),
                    accurate_mass: m as f32,
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
        let prob_cleavage_sites = cleave_at.iter().map(|&r| prob_of(r)).sum();

        SearchEngine {
            model,
            db,
            index,
            mods,
            params,
            alphabet,
            cleave_at,
            cleavage_mode,
            prob_cleavage_sites,
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
            .map_init(DpScratch::default, |scratch, (i, spec)| {
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
        scratch: &mut DpScratch,
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
        scratch: &mut DpScratch,
        spec_index: usize,
        spec: &Spectrum,
        mz: f64,
        charge: i32,
    ) -> Vec<Psm> {
        // Neutral precursor mass (= candidate peptide mass, water included).
        let parent_mass = mz as f32 * charge as f32 - charge as f32 * mass::PROTON as f32;
        let pep_nominal = scaling::nominal_bin(parent_mass - mass::WATER as f32);
        if !(50..=10_000).contains(&pep_nominal) {
            return Vec::new();
        }
        let (ti_lo, ti_hi) = self.params.isotope_errors;

        // --- the per-spectrum half: preprocess, score, and build the generating function once ---
        let peaks: Vec<(f32, f32)> = spec
            .peaks
            .iter()
            .map(|p| (p.mz as f32, p.intensity as f32))
            .collect();
        let ranked = preprocess(self.model, charge, parent_mass, &peaks);
        let scored = ScoredSpectrum::from_ranked_peaks(self.model, charge, parent_mass, ranked);

        // An isotope error of +k means the measured precursor is ~k Da high, so the true peptide
        // mass is k nominal bins lower.
        let sinks: Vec<i32> = (pep_nominal - ti_hi..=pep_nominal - ti_lo)
            .filter(|&p| p > 0)
            .collect();
        let Some(&max_p) = sinks.iter().max() else {
            return Vec::new();
        };
        let tables = scored.tables(max_p);
        let peptide_cleavage = match self.cleavage_mode {
            CleavageMode::CTerminal => PeptideCleavage {
                cleave_at: &self.cleave_at,
                credit: CLEAVAGE_CREDIT,
                penalty: CLEAVAGE_PENALTY,
            },
            CleavageMode::Off => PeptideCleavage::NONE,
        };
        let (mut graph, _) = build_reverse_graph(
            &scored,
            &tables,
            max_p,
            &[max_p],
            &self.alphabet,
            peptide_cleavage,
        );
        // The *neighbouring* residue's cleavage is probabilistic (we do not know it a priori), so
        // it weights the final distribution rather than an edge.
        let cleavage = match self.cleavage_mode {
            CleavageMode::CTerminal => Some(Cleavage {
                credit: CLEAVAGE_CREDIT,
                penalty: CLEAVAGE_PENALTY,
                prob_cleavage_sites: self.prob_cleavage_sites,
            }),
            CleavageMode::Off => None,
        };
        let mut gfs: Vec<GenFunc> = Vec::with_capacity(sinks.len());
        for &p in &sinks {
            graph.recompute_node_scores(&tables, p, &[p]);
            if let Some(gf) = compute_into(scratch, &graph, &[p as usize], cleavage) {
                gfs.push(gf);
            }
        }
        let Some(gf) = merge_group(&gfs) else {
            return Vec::new();
        };
        let denovo = gf.max_score();

        // --- the per-candidate half: every peptide in the precursor window is a tail lookup ---
        // Identical peptides occurring in several proteins score identically, so they are grouped
        // into one match carrying every protein occurrence — that is what decides decoy status
        // (`PLAN2.md` §1.3) and stops a repeated peptide consuming the whole top-N list.
        let mut grouped: HashMap<String, Hit> = HashMap::new();
        let mut buf = ScoreBuffers::default();
        for k in ti_lo..=ti_hi {
            let target = parent_mass as f64 - k as f64 * ISOTOPE_STEP;
            let win = self.params.precursor_tol.window_da(target);
            for cand in self.index.window(target - win, target + win) {
                let key = self.peptide_string(cand, false);
                match grouped.get_mut(&key) {
                    Some(hit) => hit.proteins.push(cand.protein),
                    None => {
                        let raw_score = self.raw_score(&scored, cand, &mut buf);
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
            return Vec::new();
        }
        let mut hits: Vec<(String, Hit)> = grouped.into_iter().collect();
        // Highest RawScore first; the SpecEValue tail is monotone in RawScore, so for a single
        // spectrum this is also best-SpecEValue order. The peptide key breaks ties deterministically.
        hits.sort_by(|a, b| b.1.raw_score.cmp(&a.1.raw_score).then(a.0.cmp(&b.0)));
        hits.truncate(self.params.num_matches);

        let db_size = self.db_size();
        hits.into_iter()
            .map(|(peptide_key, hit)| {
                let spec_evalue = gf.spectral_probability(hit.raw_score);
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

    /// MS-GF+ RawScore for one candidate: the node+edge match score (`DBScanScorer.getScore`) plus
    /// the terminal cleavage credit/penalty `DBScanner` adds on top.
    fn raw_score(&self, scored: &ScoredSpectrum, cand: &Candidate, buf: &mut ScoreBuffers) -> i32 {
        buf.fill(cand, self.db, self.mods);
        scored.raw_score(&buf.nominal, &buf.accurate, buf.n_mods) + self.terminal_cleavage(cand)
    }

    /// Cleavage contribution at the two termini, using the candidate's real protein context.
    ///
    /// - **C-terminal (peptide) cleavage:** credit only when the last residue is an enzyme cleavage
    ///   residue. Ending the protein does **not** substitute — verified against MS-GF+ on the F13
    ///   corpus, where protein-C-terminal peptides ending in a non-K/R residue (`R.PILVPL.-`,
    ///   `R.GCAFTM+15.995.-`) take the penalty.
    /// - **N-terminal (neighbouring) cleavage:** credit when the preceding residue is a cleavage
    ///   residue, when the peptide starts the protein (there is no preceding residue to fail the
    ///   test), or when it starts just after an excised initiator methionine — all verified
    ///   against F13.
    fn terminal_cleavage(&self, cand: &Candidate) -> i32 {
        if self.cleavage_mode == CleavageMode::Off {
            return 0;
        }
        let protein = &self.db.proteins[cand.protein as usize];
        let (start, len) = (cand.start as usize, cand.len as usize);
        let at_prot_n = start == protein.start;
        let after_initiator_met = start == protein.start + 1 && self.db.seq[protein.start] == b'M';
        let credit_if = |b: bool| if b { CLEAVAGE_CREDIT } else { CLEAVAGE_PENALTY };

        let n_term = credit_if(
            at_prot_n || after_initiator_met || self.cleave_at.contains(&self.db.seq[start - 1]),
        );
        let c_term = credit_if(self.cleave_at.contains(&self.db.seq[start + len - 1]));
        n_term + c_term
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

/// Reusable per-candidate buffers, so scoring millions of candidates does no repeated allocation.
#[derive(Default)]
struct ScoreBuffers {
    residues: Vec<Residue>,
    nominal: Vec<i32>,
    accurate: Vec<f64>,
    n_mods: i32,
}

impl ScoreBuffers {
    /// Materialise the candidate's residues (with all mod deltas) and its cumulative nominal and
    /// accurate prefix masses. The arithmetic mirrors `msgf_chem::peptide::{nominal_prefix_masses,
    /// accurate_prefix_masses}` operation-for-operation, so scores stay bit-identical to the
    /// allocating path those functions provide.
    fn fill(&mut self, cand: &Candidate, db: &ProteinDb, mods: &ModSet) {
        let protein = &db.proteins[cand.protein as usize];
        let (start, len) = (cand.start as usize, cand.len as usize);
        let seq = &db.seq[start..start + len];
        let at_prot_n = start == protein.start;
        let at_prot_c = start + len == protein.start + protein.len;

        self.residues.clear();
        self.nominal.clear();
        self.accurate.clear();
        self.n_mods = 0;

        let (mut cum_nominal, mut cum_accurate) = (0i32, 0.0f64);
        for (i, &r) in seq.iter().enumerate() {
            let delta = mods.fixed_delta(r, i, len, at_prot_n, at_prot_c)
                + cand.placement.delta_at(i, mods);
            if delta != 0.0 {
                self.n_mods += 1;
            }
            let m = msgf_chem::residue_mass(r).expect("candidates hold standard residues") + delta;
            cum_nominal += scaling::nominal_bin(m as f32);
            cum_accurate += m;
            self.residues.push(Residue {
                aa: r,
                mod_delta: delta,
            });
            self.nominal.push(cum_nominal);
            self.accurate.push(cum_accurate);
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
    fn score_buffers_match_the_allocating_helpers() {
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

        // The same peptide via the string parser + the allocating prefix-mass helpers.
        let delta = mods.mods[0].mass;
        let pep = format!("SAM{delta:+}PLERK");
        let residues = msgf_chem::peptide::parse(&pep).unwrap();
        assert_eq!(
            buf.nominal,
            msgf_chem::peptide::nominal_prefix_masses(&residues)
        );
        assert_eq!(
            buf.accurate,
            msgf_chem::peptide::accurate_prefix_masses(&residues)
        );
        assert_eq!(buf.n_mods, msgf_chem::peptide::num_mods(&residues) as i32);
        assert_eq!(buf.residues.len(), 8);
    }
}
