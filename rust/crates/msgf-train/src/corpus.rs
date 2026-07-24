//! Training-corpus ingest — annotated MGF (peptide-labelled spectra).
//!
//! The corpus format is an MGF whose every spectrum carries a peptide annotation in `SEQ=`
//! (MassIVE-KB peptide-library MGFs are exactly this), with inline modification deltas:
//!
//! ```text
//! BEGIN IONS
//! PEPMASS=581.3062501893
//! CHARGE=3
//! SEQ=+42.011AAAADSFSGGPAGVRLPR
//! 101.0710678100586  15761.0380859375
//! ...
//! END IONS
//! ```
//!
//! A leading `+d`/`-d` (before the first residue) is an N-terminal modification; it is folded into
//! the first residue's delta, which is equivalent for every prefix/suffix mass we count.

use msgf_chem::mass;
use msgf_chem::peptide::Residue;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// One training example: a spectrum with a confident peptide identification.
#[derive(Debug, Clone)]
pub struct TrainingPsm {
    pub charge: i32,
    /// De-charged neutral precursor mass `M` (the scorer's `parent_mass`).
    pub parent_mass: f32,
    pub residues: Vec<Residue>,
    /// Raw `(m/z, intensity)` peaks in file order (m/z ascending).
    pub peaks: Vec<(f32, f32)>,
}

/// Acceptance rules applied while reading a corpus. Every rejection is counted so a training run
/// can report exactly what it kept.
#[derive(Debug, Clone)]
pub struct CorpusFilter {
    pub charge_min: i32,
    pub charge_max: i32,
    pub len_min: usize,
    pub len_max: usize,
    pub min_peaks: usize,
    /// Require a tryptic C-terminus (K/R) — set from the enzyme identity being trained.
    pub require_tryptic_cterm: bool,
    /// Max |computed peptide mass − precursor mass| (Da) before the annotation is distrusted.
    pub precursor_mass_tol: f64,
}

impl Default for CorpusFilter {
    fn default() -> Self {
        Self {
            charge_min: 2,
            charge_max: 8,
            len_min: 6,
            len_max: 50,
            min_peaks: 20,
            require_tryptic_cterm: true,
            precursor_mass_tol: 0.1,
        }
    }
}

/// Per-reason rejection tally (reported by the CLI so corpus quality is visible).
#[derive(Debug, Default, Clone)]
pub struct CorpusStats {
    pub read: usize,
    pub kept: usize,
    pub no_annotation: usize,
    pub bad_residue: usize,
    pub charge: usize,
    pub length: usize,
    pub peaks: usize,
    pub cterm: usize,
    pub mass_mismatch: usize,
}

impl CorpusStats {
    fn merge(&mut self, o: &CorpusStats) {
        self.read += o.read;
        self.kept += o.kept;
        self.no_annotation += o.no_annotation;
        self.bad_residue += o.bad_residue;
        self.charge += o.charge;
        self.length += o.length;
        self.peaks += o.peaks;
        self.cterm += o.cterm;
        self.mass_mismatch += o.mass_mismatch;
    }
}

/// Parse a `SEQ=` peptide with inline `+d`/`-d` modification deltas. A delta before the first
/// residue (N-terminal mod) is folded into the first residue. `None` on a non-standard residue.
pub fn parse_seq(seq: &str) -> Option<Vec<Residue>> {
    let b = seq.as_bytes();
    let mut out: Vec<Residue> = Vec::with_capacity(b.len());
    let mut nterm = 0.0f64;
    let mut i = 0usize;
    while i < b.len() {
        let c = b[i];
        if c == b'+' || c == b'-' {
            let start = i;
            i += 1;
            while i < b.len() && (b[i].is_ascii_digit() || b[i] == b'.') {
                i += 1;
            }
            let d: f64 = seq[start..i].parse().ok()?;
            match out.last_mut() {
                Some(r) => r.mod_delta += d,
                None => nterm += d,
            }
        } else if c.is_ascii_alphabetic() {
            msgf_chem::residue_mass(c)?; // reject non-standard residues (B, X, Z, U, …)
            out.push(Residue {
                aa: c,
                mod_delta: 0.0,
            });
            i += 1;
        } else if c.is_ascii_whitespace() {
            i += 1;
        } else {
            return None;
        }
    }
    if out.is_empty() {
        return None;
    }
    out[0].mod_delta += nterm;
    Some(out)
}

/// Read one annotated MGF, applying `filter`. Spectra without `SEQ=` are skipped.
pub fn read_annotated_mgf<P: AsRef<Path>>(
    path: P,
    filter: &CorpusFilter,
    out: &mut Vec<TrainingPsm>,
    stats: &mut CorpusStats,
) -> io::Result<()> {
    let f = File::open(path)?;
    let mut rd = BufReader::with_capacity(1 << 20, f);
    let mut line = String::new();

    let mut seq: Option<String> = None;
    let mut charge: i32 = 0;
    let mut pepmass: f64 = 0.0;
    let mut peaks: Vec<(f32, f32)> = Vec::new();
    let mut local = CorpusStats::default();

    loop {
        line.clear();
        if rd.read_line(&mut line)? == 0 {
            break;
        }
        let l = line.trim_end_matches(['\r', '\n']);
        if l.is_empty() {
            continue;
        }
        if l == "BEGIN IONS" {
            seq = None;
            charge = 0;
            pepmass = 0.0;
            peaks.clear();
        } else if l == "END IONS" {
            local.read += 1;
            finish(
                seq.take(),
                charge,
                pepmass,
                &mut peaks,
                filter,
                out,
                &mut local,
            );
        } else if let Some(v) = l.strip_prefix("SEQ=") {
            seq = Some(v.trim().to_string());
        } else if let Some(v) = l.strip_prefix("CHARGE=") {
            let v = v.trim().trim_end_matches(['+', '-']);
            charge = v.parse().unwrap_or(0);
        } else if let Some(v) = l.strip_prefix("PEPMASS=") {
            pepmass = v
                .split_whitespace()
                .next()
                .and_then(|t| t.parse().ok())
                .unwrap_or(0.0);
        } else if l.as_bytes()[0].is_ascii_digit() {
            let mut it = l.split_whitespace();
            if let (Some(a), Some(b)) = (it.next(), it.next()) {
                if let (Ok(m), Ok(i)) = (a.parse::<f32>(), b.parse::<f32>()) {
                    peaks.push((m, i));
                }
            }
        }
    }
    stats.merge(&local);
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn finish(
    seq: Option<String>,
    charge: i32,
    pepmass: f64,
    peaks: &mut Vec<(f32, f32)>,
    filter: &CorpusFilter,
    out: &mut Vec<TrainingPsm>,
    st: &mut CorpusStats,
) {
    let Some(seq) = seq else {
        st.no_annotation += 1;
        return;
    };
    if charge < filter.charge_min || charge > filter.charge_max {
        st.charge += 1;
        return;
    }
    if peaks.len() < filter.min_peaks {
        st.peaks += 1;
        return;
    }
    let Some(residues) = parse_seq(&seq) else {
        st.bad_residue += 1;
        return;
    };
    if residues.len() < filter.len_min || residues.len() > filter.len_max {
        st.length += 1;
        return;
    }
    if filter.require_tryptic_cterm {
        let c = residues.last().unwrap().aa.to_ascii_uppercase();
        if c != b'K' && c != b'R' {
            st.cterm += 1;
            return;
        }
    }
    // Peptide neutral mass from the annotation vs the de-charged precursor.
    let pep_mass: f64 = residues
        .iter()
        .map(|r| msgf_chem::residue_mass(r.aa).unwrap_or(0.0) + r.mod_delta)
        .sum::<f64>()
        + mass::WATER;
    let obs_mass = pepmass * charge as f64 - charge as f64 * mass::PROTON;
    if (pep_mass - obs_mass).abs() > filter.precursor_mass_tol {
        st.mass_mismatch += 1;
        return;
    }

    let mut pk = std::mem::take(peaks);
    pk.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap()
            .then(a.1.partial_cmp(&b.1).unwrap())
    });
    out.push(TrainingPsm {
        charge,
        // Use the *annotated* peptide mass: it is the exact theoretical parent mass the scorer
        // would see for this identification, free of precursor-picking error.
        parent_mass: pep_mass as f32,
        residues,
        peaks: pk,
    });
    st.kept += 1;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_nterm_and_residue_mods() {
        let r = parse_seq("+42.011AAM+15.995K").unwrap();
        assert_eq!(r.len(), 4);
        assert!((r[0].mod_delta - 42.011).abs() < 1e-9); // N-term folded into residue 1
        assert!((r[2].mod_delta - 15.995).abs() < 1e-9);
        assert!(parse_seq("PEPXIDEK").is_none()); // X is not a standard residue
    }
}
