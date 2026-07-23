//! Peptide parsing and prefix-mass arrays for scoring.
//!
//! Handles MS-GF+'s peptide string form: optional `X.`/`.X` enzymatic context and inline
//! `+delta`/`-delta` modification masses (e.g. `R.RTLMARPM+15.995IKEAR.M`). Produces the
//! cumulative nominal and accurate prefix-mass arrays `FastScorer`/`DBScanScorer` consume.

use crate::{residue_mass, scaling};

/// One residue with an optional modification delta (Da).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Residue {
    pub aa: u8,
    pub mod_delta: f64,
}

/// Strip an optional single-residue `X.` prefix and `.X` suffix (the `-` terminus is allowed).
pub fn strip_context(pep: &str) -> &str {
    let b = pep.as_bytes();
    if b.len() >= 4 && b[1] == b'.' && b[b.len() - 2] == b'.' {
        &pep[2..pep.len() - 2]
    } else {
        pep
    }
}

/// Parse a peptide into residues (+mod deltas). `None` on a non-standard residue or bad mod.
pub fn parse(pep: &str) -> Option<Vec<Residue>> {
    let core = strip_context(pep);
    let b = core.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < b.len() {
        let c = b[i];
        if !c.is_ascii_alphabetic() {
            return None;
        }
        residue_mass(c)?; // reject non-standard residues
        i += 1;
        let mut delta = 0.0f64;
        while i < b.len() && (b[i] == b'+' || b[i] == b'-') {
            let start = i;
            i += 1;
            while i < b.len() && (b[i].is_ascii_digit() || b[i] == b'.') {
                i += 1;
            }
            delta += core[start..i].parse::<f64>().ok()?;
        }
        out.push(Residue {
            aa: c,
            mod_delta: delta,
        });
    }
    (!out.is_empty()).then_some(out)
}

/// Number of modified residues (inline mods), i.e. MS-GF+'s `numMods`.
pub fn num_mods(residues: &[Residue]) -> usize {
    residues.iter().filter(|r| r.mod_delta != 0.0).count()
}

/// Cumulative **nominal** prefix masses — per-residue nominal masses summed, last element = full
/// peptide nominal mass. Matches `AminoAcid.getNominalMass` accumulation in MS-GF+.
pub fn nominal_prefix_masses(residues: &[Residue]) -> Vec<i32> {
    let mut out = Vec::with_capacity(residues.len());
    let mut cum = 0i32;
    for r in residues {
        let m = residue_mass(r.aa).expect("standard residue") as f64 + r.mod_delta;
        cum += scaling::nominal_bin(m as f32);
        out.push(cum);
    }
    out
}

/// Cumulative **accurate** (monoisotopic) prefix residue masses, last = full peptide residue mass.
pub fn accurate_prefix_masses(residues: &[Residue]) -> Vec<f64> {
    let mut out = Vec::with_capacity(residues.len());
    let mut cum = 0.0f64;
    for r in residues {
        cum += residue_mass(r.aa).expect("standard residue") as f64 + r.mod_delta;
        out.push(cum);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strips_context() {
        assert_eq!(strip_context("K.RSRRRRKR.A"), "RSRRRRKR");
        assert_eq!(strip_context("-.MPKRK.S"), "MPKRK");
        assert_eq!(strip_context("RSRRRRKR"), "RSRRRRKR");
    }

    #[test]
    fn parses_mods() {
        let r = parse("R.RTLMARPM+15.995IKEAR.M").unwrap();
        assert_eq!(r.len(), 13); // RTLMARPMIKEAR
        assert_eq!(num_mods(&r), 1);
        let m = r.iter().find(|x| x.mod_delta != 0.0).unwrap();
        assert_eq!(m.aa, b'M');
        assert!((m.mod_delta - 15.995).abs() < 1e-9);
        assert!(parse("K.PEPXIDE.R").is_none()); // X not standard
    }

    #[test]
    fn nominal_prefix_is_cumulative() {
        let r = parse("PEPTIDE").unwrap();
        let n = nominal_prefix_masses(&r);
        assert_eq!(n.len(), 7);
        // strictly increasing, last = sum of all 7 residue nominal masses
        assert!(n.windows(2).all(|w| w[1] > w[0]));
        let expect: i32 = r
            .iter()
            .map(|x| scaling::nominal_bin(residue_mass(x.aa).unwrap() as f32))
            .sum();
        assert_eq!(*n.last().unwrap(), expect);
    }
}
