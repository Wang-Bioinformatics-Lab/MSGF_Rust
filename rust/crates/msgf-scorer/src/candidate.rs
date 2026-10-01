//! RawScore of one candidate peptide: vertex scores + edge terms + enzymatic-terminus terms.

use super::cleavage::Cleavage;
use super::spectrum::PreparedSpectrum;
use super::{nominal, residue_mass};

/// One residue of a candidate: its (unmodified) letter and total modification delta.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Residue {
    pub letter: u8,
    /// Fixed + variable modification mass (Da), 0 if unmodified.
    pub delta: f64,
}

/// A candidate peptide with its N-terminal cleavage context decided by the caller.
#[derive(Clone, Copy, Debug)]
pub struct Candidate<'a> {
    pub residues: &'a [Residue],
    /// True when the N-terminus is enzymatic (preceding residue is a cleavage site, protein
    /// start, or the caller's own rule such as initiator-Met excision).
    pub n_term_credit: bool,
}

/// RawScore = Match + T_N + T_C. `None` if the candidate is empty or has a non-standard letter.
pub fn match_and_terminal_score(
    prep: &PreparedSpectrum,
    cand: &Candidate,
    cleavage: &Cleavage,
) -> Option<i32> {
    let (m, t) = match_and_terminal_parts(prep, cand, cleavage)?;
    Some(m + t)
}

/// `(Match, T_N + T_C)` separately (for per-stage checks).
pub fn match_and_terminal_parts(
    prep: &PreparedSpectrum,
    cand: &Candidate,
    cleavage: &Cleavage,
) -> Option<(i32, i32)> {
    let res = cand.residues;
    if res.is_empty() {
        return None;
    }
    let l = res.len();
    // Real masses, nominal prefix sums N_j and binary64 cumulative real sums A_j.
    let mut cum_nom = Vec::with_capacity(l + 1);
    let mut cum_real = Vec::with_capacity(l + 1);
    cum_nom.push(0i32);
    cum_real.push(0f64);
    for r in res {
        let m = residue_mass(r.letter)? + r.delta;
        cum_nom.push(cum_nom.last().unwrap() + nominal(m));
        cum_real.push(cum_real.last().unwrap() + m);
    }
    let p = cum_nom[l];

    let mut score = 0i32;
    if l > 1 {
        for j in 1..l {
            score += prep.vertex(p, cum_nom[j]);
        }
        let suffix_anchor = !prep.anchor.map_or(false, |a| a.prefix);
        for j in 1..l {
            score += if suffix_anchor {
                let ma = (cum_real[j + 1] - cum_real[j]) as f32;
                prep.edge(p - cum_nom[j], p - cum_nom[j + 1], ma)
            } else {
                let ma = (cum_real[j] - cum_real[j - 1]) as f32;
                prep.edge(cum_nom[j], cum_nom[j - 1], ma)
            };
        }
    }

    let term = if cleavage.enabled {
        let c = if cleavage.sites.contains(&res[l - 1].letter) {
            cleavage.credit
        } else {
            cleavage.penalty
        };
        let n = if cand.n_term_credit {
            cleavage.credit
        } else {
            cleavage.penalty
        };
        c + n
    } else {
        0
    };
    Some((score, term))
}

/// Cumulative nominal masses N_1..N_L (per-stage checks).
pub fn cumulative_nominal(res: &[Residue]) -> Option<Vec<i32>> {
    let mut acc = 0;
    res.iter()
        .map(|r| {
            acc += nominal(residue_mass(r.letter)? + r.delta);
            Some(acc)
        })
        .collect()
}
