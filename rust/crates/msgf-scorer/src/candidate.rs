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
    // Prefix sums live on the stack for ordinary lengths; longer candidates use the heap.
    const STACK: usize = 128;
    if l < STACK {
        let mut cum_nom = [0i32; STACK];
        let mut cum_real = [0f64; STACK];
        parts_with(
            prep,
            cand,
            cleavage,
            &mut cum_nom[..=l],
            &mut cum_real[..=l],
        )
    } else {
        let mut cum_nom = vec![0i32; l + 1];
        let mut cum_real = vec![0f64; l + 1];
        parts_with(prep, cand, cleavage, &mut cum_nom, &mut cum_real)
    }
}

/// [`match_and_terminal_parts`] with caller-provided prefix-sum buffers of length `L + 1`.
#[inline(always)]
fn parts_with(
    prep: &PreparedSpectrum,
    cand: &Candidate,
    cleavage: &Cleavage,
    cum_nom: &mut [i32],
    cum_real: &mut [f64],
) -> Option<(i32, i32)> {
    let res = cand.residues;
    let l = res.len();
    // Real masses, nominal prefix sums N_j and binary64 cumulative real sums A_j.
    cum_nom[0] = 0;
    cum_real[0] = 0.0;
    for (j, r) in res.iter().enumerate() {
        let m = residue_mass(r.letter)? + r.delta;
        cum_nom[j + 1] = cum_nom[j] + nominal(m);
        cum_real[j + 1] = cum_real[j] + m;
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
