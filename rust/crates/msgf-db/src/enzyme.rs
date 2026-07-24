//! Enzymes and in-silico digestion.
//!
//! An enzyme is `(residues it cleaves at, which side of them it cuts)`. Digestion walks a protein,
//! marks the cleavage sites, and emits the peptides spanning consecutive sites (plus up to
//! `max_missed_cleavages` interior ones). Semi- and non-enzymatic searches are expressed through
//! [`DigestParams::min_termini`] (MS-GF+'s "number of tolerable termini", `-ntt`).
//!
//! Enzyme definitions follow the documented `enzymes.txt` interchange format
//! (`ShortName,CleaveAt,Terminus,Description`); the numbering of the built-ins matches MS-GF+'s
//! `-e` argument so a command line transfers over unchanged.

use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// A cleavage specificity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Enzyme {
    pub name: String,
    /// Residues the enzyme cleaves at. **Empty means unspecific** (cleaves anywhere).
    pub cleave_at: Vec<u8>,
    /// `true` = cuts C-terminal to `cleave_at` (trypsin), `false` = N-terminal (Lys-N, Asp-N).
    pub c_term: bool,
    pub description: String,
}

impl Enzyme {
    /// The built-in enzymes, indexed by MS-GF+'s `-e` number (`0` = unspecific cleavage).
    pub fn builtin(id: u32) -> Option<Enzyme> {
        let (name, cleave_at, c_term, description) = match id {
            0 => ("UnspecificCleavage", "", true, "unspecific cleavage"),
            1 => ("Tryp", "KR", true, "Trypsin"),
            2 => ("Chymotrypsin", "FYWL", true, "Chymotrypsin"),
            3 => ("LysC", "K", true, "Lys-C"),
            4 => ("LysN", "K", false, "Lys-N"),
            5 => ("GluC", "E", true, "glutamyl endopeptidase"),
            6 => ("ArgC", "R", true, "Arg-C"),
            7 => ("AspN", "D", false, "Asp-N"),
            8 => ("aLP", "", true, "alphaLP"),
            9 => ("NoCleavage", "", true, "no cleavage (endogenous peptides)"),
            10 => ("TrypPlusC", "KRC", true, "Cleave after K, R, or C"),
            _ => return None,
        };
        Some(Enzyme {
            name: name.to_string(),
            cleave_at: cleave_at.as_bytes().to_vec(),
            c_term,
            description: description.to_string(),
        })
    }

    /// Resolve an enzyme from a command-line spec: either a built-in number (`1`), a built-in short
    /// name (`Tryp`, case-insensitive), or a full inline definition
    /// (`ShortName,CleaveAt,Terminus,Description`).
    pub fn parse(spec: &str) -> Result<Enzyme, String> {
        let spec = spec.trim();
        if spec.contains(',') {
            return Self::parse_definition(spec);
        }
        if let Ok(n) = spec.parse::<u32>() {
            return Self::builtin(n)
                .ok_or_else(|| format!("unknown enzyme number `{n}` (valid: 0-10)"));
        }
        (0..=10)
            .filter_map(Self::builtin)
            .find(|e| e.name.eq_ignore_ascii_case(spec))
            .ok_or_else(|| format!("unknown enzyme `{spec}` (use a number 0-10, a short name, or a full definition)"))
    }

    /// Parse one `ShortName,CleaveAt,Terminus,Description` line. `CleaveAt` of `null` (or empty)
    /// means unspecific; `Terminus` is `C` or `N`.
    pub fn parse_definition(line: &str) -> Result<Enzyme, String> {
        let f: Vec<&str> = line.splitn(4, ',').map(str::trim).collect();
        if f.len() < 3 {
            return Err(format!(
                "enzyme `{line}`: expected ShortName,CleaveAt,Terminus[,Description]"
            ));
        }
        let cleave_at = if f[1].eq_ignore_ascii_case("null") || f[1].is_empty() {
            Vec::new()
        } else {
            if !f[1].bytes().all(|b| b.is_ascii_uppercase()) {
                return Err(format!(
                    "enzyme `{}`: CleaveAt must be upper-case residues",
                    f[0]
                ));
            }
            f[1].as_bytes().to_vec()
        };
        let c_term = match f[2].to_ascii_uppercase().as_str() {
            "C" => true,
            "N" => false,
            other => {
                return Err(format!(
                    "enzyme `{}`: Terminus must be C or N, got `{other}`",
                    f[0]
                ))
            }
        };
        Ok(Enzyme {
            name: f[0].to_string(),
            cleave_at,
            c_term,
            description: f.get(3).unwrap_or(&"").to_string(),
        })
    }

    /// Read an `enzymes.txt`-format file (`#` comments, one definition per line).
    pub fn read_file(path: &Path) -> io::Result<Vec<Enzyme>> {
        let mut out = Vec::new();
        for line in BufReader::new(File::open(path)?).lines() {
            let line = line?;
            let body = line.split('#').next().unwrap_or("").trim();
            if body.is_empty() {
                continue;
            }
            match Enzyme::parse_definition(body) {
                Ok(e) => out.push(e),
                Err(e) => return Err(io::Error::new(io::ErrorKind::InvalidData, e)),
            }
        }
        Ok(out)
    }

    /// Whether this enzyme has no specificity (cleaves anywhere) — `aLP`, `NoCleavage`, and
    /// `UnspecificCleavage`. Such a search has no enzymatic termini, so cleavage scoring is off.
    #[inline]
    pub fn is_unspecific(&self) -> bool {
        self.cleave_at.is_empty()
    }

    #[inline]
    fn cleaves(&self, residue: u8) -> bool {
        self.cleave_at.contains(&residue)
    }

    /// Whether the bond **between** `seq[i - 1]` and `seq[i]` is cut by this enzyme. `i` ranges over
    /// `1..seq.len()`; the protein termini are handled separately by [`digest`].
    #[inline]
    fn cuts_before(&self, seq: &[u8], i: usize) -> bool {
        if self.is_unspecific() {
            return true;
        }
        if self.c_term {
            self.cleaves(seq[i - 1])
        } else {
            self.cleaves(seq[i])
        }
    }
}

/// "No limit on missed cleavages" — MS-GF+'s default. Peptide length still bounds the search.
pub const UNLIMITED_MISSED_CLEAVAGES: usize = usize::MAX;

/// Digestion settings.
#[derive(Debug, Clone)]
pub struct DigestParams {
    pub enzyme: Enzyme,
    /// Maximum internal cleavage sites a peptide may span. **MS-GF+ does not limit this by
    /// default** ([`UNLIMITED_MISSED_CLEAVAGES`]) — its `-maxMissedCleavages` defaults to "no
    /// limit", which is why its own results are full of long K/R-rich peptides. Reproducing a
    /// MS-GF+ run therefore needs the unlimited setting; a small number (2) is the conventional
    /// choice for a normal search and keeps the index far smaller.
    pub max_missed_cleavages: usize,
    pub min_len: usize,
    pub max_len: usize,
    /// Minimum number of enzymatic termini a peptide must have: `2` = fully enzymatic (MS-GF+
    /// default), `1` = semi-enzymatic, `0` = non-enzymatic. Protein termini always count as
    /// enzymatic, as in MS-GF+.
    pub min_termini: u8,
    /// Allow a peptide to start at residue 2 of a protein whose first residue is `M`, treating that
    /// as an enzymatic N-terminus. This is initiator-methionine excision, which MS-GF+ models: its
    /// F13 results contain peptides flanked `M.` that are otherwise unreachable.
    pub cleave_initiator_met: bool,
}

impl Default for DigestParams {
    fn default() -> Self {
        Self {
            enzyme: Enzyme::builtin(1).expect("trypsin is built in"),
            max_missed_cleavages: UNLIMITED_MISSED_CLEAVAGES,
            min_len: 6,
            max_len: 40,
            min_termini: 2,
            cleave_initiator_met: true,
        }
    }
}

/// One digested peptide, addressed as an offset+length within its protein's sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DigestedPeptide {
    /// Offset of the peptide's first residue within the protein sequence.
    pub start: usize,
    pub len: usize,
    /// How many enzymatic termini this peptide has (0, 1, or 2).
    pub n_termini: u8,
}

/// Digest one protein sequence, calling `emit` for every peptide that satisfies `params`.
///
/// Peptides containing a non-standard residue are skipped — their mass is undefined, so they can
/// neither be matched to a precursor nor scored.
pub fn digest(seq: &[u8], params: &DigestParams, mut emit: impl FnMut(DigestedPeptide)) {
    let n = seq.len();
    if n == 0 {
        return;
    }
    // Cut points, as offsets where a peptide may start: 0, every enzymatic bond, and n (the end).
    let mut cuts: Vec<usize> = Vec::with_capacity(n / 8 + 2);
    cuts.push(0);
    // Initiator-methionine excision: a protein starting with M may also be cleaved after it.
    if params.cleave_initiator_met && seq[0] == b'M' && n > 1 && !params.enzyme.cuts_before(seq, 1)
    {
        cuts.push(1);
    }
    for i in 1..n {
        if params.enzyme.cuts_before(seq, i) {
            cuts.push(i);
        }
    }
    cuts.sort_unstable();
    cuts.dedup();
    cuts.push(n);

    let standard = |s: &[u8]| s.iter().all(|&b| msgf_chem::residue_mass(b).is_some());

    if params.enzyme.is_unspecific() {
        // No specificity: every bond is a cut, so "missed cleavages" and the termini count carry no
        // information. Enumerate every substring in the length range, as MS-GF+ `-e 0` does.
        for start in 0..n {
            for end in (start + params.min_len.max(1))..=n.min(start + params.max_len) {
                if standard(&seq[start..end]) {
                    emit(DigestedPeptide {
                        start,
                        len: end - start,
                        n_termini: 2,
                    });
                }
            }
        }
        return;
    }

    if params.min_termini >= 2 {
        // Fully enzymatic: both ends are cut points, with at most `max_missed_cleavages` interior ones.
        for (ci, &start) in cuts.iter().enumerate().take(cuts.len().saturating_sub(1)) {
            // `max_len` bounds the walk, so an unlimited missed-cleavage setting is still finite.
            for &end in cuts
                .iter()
                .skip(ci + 1)
                .take(params.max_missed_cleavages.saturating_add(1))
            {
                let len = end - start;
                if len < params.min_len {
                    continue;
                }
                if len > params.max_len {
                    break; // `cuts` is ascending, so every later end is longer too
                }
                if standard(&seq[start..end]) {
                    emit(DigestedPeptide {
                        start,
                        len,
                        n_termini: 2,
                    });
                }
            }
        }
        return;
    }

    // Semi- (min_termini == 1) or non-enzymatic (0): every substring in the length range whose
    // enzymatic-terminus count reaches the threshold. `is_cut[i]` = "a peptide may start at i".
    let mut is_cut = vec![false; n + 1];
    for &c in &cuts {
        is_cut[c] = true;
    }
    for start in 0..n {
        // Bound the walk by missed cleavages so a semi-tryptic search stays proportional to the
        // enzymatic one rather than quadratic in protein length.
        let mut missed = 0usize;
        for end in (start + 1)..=n.min(start + params.max_len) {
            let len = end - start;
            if end < n && is_cut[end] {
                missed += 1;
                if missed > params.max_missed_cleavages.saturating_add(1) {
                    break;
                }
            }
            if len < params.min_len {
                continue;
            }
            let n_termini = u8::from(is_cut[start]) + u8::from(is_cut[end]);
            if n_termini < params.min_termini {
                continue;
            }
            if standard(&seq[start..end]) {
                emit(DigestedPeptide {
                    start,
                    len,
                    n_termini,
                });
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn peptides(seq: &[u8], p: &DigestParams) -> Vec<String> {
        let mut out = Vec::new();
        digest(seq, p, |d| {
            out.push(String::from_utf8(seq[d.start..d.start + d.len].to_vec()).unwrap())
        });
        out
    }

    #[test]
    fn tryptic_cuts_after_k_and_r() {
        let p = DigestParams {
            max_missed_cleavages: 0,
            min_len: 1,
            max_len: 100,
            ..Default::default()
        };
        assert_eq!(
            peptides(b"PEPTIDEKSAMPLERTAIL", &p),
            ["PEPTIDEK", "SAMPLER", "TAIL"]
        );
    }

    #[test]
    fn tryptic_cuts_before_proline_too() {
        // MS-GF+ does not apply the "no cleavage before proline" rule.
        let p = DigestParams {
            max_missed_cleavages: 0,
            min_len: 1,
            max_len: 100,
            ..Default::default()
        };
        assert_eq!(peptides(b"AAKPEEK", &p), ["AAK", "PEEK"]);
    }

    #[test]
    fn missed_cleavages_extend_peptides() {
        let p = DigestParams {
            max_missed_cleavages: 1,
            min_len: 1,
            max_len: 100,
            ..Default::default()
        };
        assert_eq!(
            peptides(b"AAKGGKCCK", &p),
            ["AAK", "AAKGGK", "GGK", "GGKCCK", "CCK"]
        );
    }

    #[test]
    fn n_terminal_enzyme_cuts_before_residue() {
        let p = DigestParams {
            enzyme: Enzyme::builtin(7).unwrap(), // AspN, cleaves N-terminal to D
            max_missed_cleavages: 0,
            min_len: 1,
            max_len: 100,
            ..Default::default()
        };
        assert_eq!(peptides(b"AADEEDFF", &p), ["AA", "DEE", "DFF"]);
    }

    #[test]
    fn nonstandard_residues_are_skipped() {
        let p = DigestParams {
            max_missed_cleavages: 0,
            min_len: 1,
            max_len: 100,
            ..Default::default()
        };
        // BXZ are not scorable residues, so AAXK never appears.
        assert_eq!(peptides(b"AAXKGGK", &p), ["GGK"]);
    }

    #[test]
    fn length_bounds_are_respected() {
        let p = DigestParams {
            max_missed_cleavages: 2,
            min_len: 6,
            max_len: 8,
            ..Default::default()
        };
        // AAK(3) too short; AAKGGK(6) ok; AAKGGKCCK(9) too long.
        assert_eq!(peptides(b"AAKGGKCCK", &p), ["AAKGGK", "GGKCCK"]);
    }

    #[test]
    fn semi_tryptic_keeps_one_enzymatic_terminus() {
        let p = DigestParams {
            max_missed_cleavages: 0,
            min_len: 3,
            max_len: 5,
            min_termini: 1,
            ..Default::default()
        };
        // AAGGKTT: cut points are the protein start (0), after K (5), and the protein end (7).
        let got = peptides(b"AAGGKTT", &p);
        assert!(got.contains(&"AAGGK".to_string())); // fully enzymatic (both termini)
        assert!(got.contains(&"AAGG".to_string())); // enzymatic N-term only (protein start)
        assert!(got.contains(&"AGGK".to_string())); // enzymatic C-term only (after K)
        assert!(!got.contains(&"AGG".to_string())); // neither terminus enzymatic -> excluded
                                                    // With min_termini = 2 only the fully enzymatic peptide survives.
        let strict = DigestParams {
            min_termini: 2,
            ..p
        };
        assert_eq!(peptides(b"AAGGKTT", &strict), ["AAGGK"]);
    }

    #[test]
    fn unlimited_missed_cleavages_is_the_default_and_is_bounded_by_length() {
        let p = DigestParams {
            min_len: 1,
            max_len: 9,
            ..Default::default()
        };
        assert_eq!(p.max_missed_cleavages, UNLIMITED_MISSED_CLEAVAGES);
        // Every tryptic sub-run up to 9 residues, however many K's it spans.
        let got = peptides(b"AAKGGKCCK", &p);
        assert!(got.contains(&"AAKGGKCCK".to_string()), "{got:?}");
        assert!(got.contains(&"AAKGGK".to_string()));
        assert!(got.contains(&"AAK".to_string()));
        // The length bound still applies.
        let short = DigestParams { max_len: 6, ..p };
        assert!(!peptides(b"AAKGGKCCK", &short).contains(&"AAKGGKCCK".to_string()));
    }

    #[test]
    fn initiator_methionine_is_excised() {
        let p = DigestParams {
            min_len: 3,
            max_len: 20,
            ..Default::default()
        };
        // MAAGGK: the initiator M may be removed, so AAGGK is a valid fully-enzymatic peptide.
        let got = peptides(b"MAAGGKTTR", &p);
        assert!(got.contains(&"MAAGGK".to_string()), "{got:?}");
        assert!(got.contains(&"AAGGK".to_string()), "{got:?}");
        // Switched off, only the untruncated form survives.
        let off = DigestParams {
            cleave_initiator_met: false,
            ..p
        };
        let got = peptides(b"MAAGGKTTR", &off);
        assert!(got.contains(&"MAAGGK".to_string()));
        assert!(!got.contains(&"AAGGK".to_string()));
    }

    #[test]
    fn initiator_excision_does_not_duplicate_a_real_cut() {
        // A protein starting `MK...` already has a tryptic cut after position 1; adding the
        // initiator cut must not produce duplicate cut points (and so duplicate peptides).
        let p = DigestParams {
            min_len: 1,
            max_len: 20,
            max_missed_cleavages: 0,
            ..Default::default()
        };
        let got = peptides(b"MKAAGGR", &p);
        assert_eq!(
            got.len(),
            got.iter().collect::<std::collections::HashSet<_>>().len(),
            "{got:?}"
        );
    }

    #[test]
    fn unspecific_enzyme_cuts_everywhere() {
        let p = DigestParams {
            enzyme: Enzyme::builtin(0).unwrap(),
            max_missed_cleavages: 0,
            min_len: 3,
            max_len: 3,
            ..Default::default()
        };
        // Every 3-mer of a 5-residue protein: AGC, GCD, CDE.
        assert_eq!(peptides(b"AGCDE", &p), ["AGC", "GCD", "CDE"]);
    }

    #[test]
    fn enzyme_specs_parse() {
        assert_eq!(Enzyme::parse("1").unwrap().name, "Tryp");
        assert_eq!(Enzyme::parse("tryp").unwrap().cleave_at, b"KR");
        assert!(Enzyme::parse("9").unwrap().is_unspecific());
        let e = Enzyme::parse("CNBr,M,C,CNBr cleavage").unwrap();
        assert_eq!(e.cleave_at, b"M");
        assert!(e.c_term);
        assert!(Enzyme::parse("Nope").is_err());
    }
}
