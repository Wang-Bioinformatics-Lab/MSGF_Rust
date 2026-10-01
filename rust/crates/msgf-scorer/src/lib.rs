//! # msgf-scorer — the `.param` scoring model and the RawScore of a peptide-spectrum match
//!
//! Two halves:
//!
//! - **The model file.** [`ScoringModel`] is the `.param` format as data: [`read_param`] decodes
//!   it, [`write_param`] encodes it (the trainer, `msgf-train`, emits models through it). The
//!   format is documented in `docs/param-format.md`. [`bundled`] embeds the MIT/CC0 model this
//!   project trained itself.
//! - **Scoring.** [`ScoreModel`] is a model decoded for scoring (derived log-likelihood tables).
//!   [`prepare`] turns one MS/MS spectrum at one charge into a [`PreparedSpectrum`] (peak order,
//!   precursor suppression, ranks, isotope-cluster reduction, per-segment partitions, node and
//!   edge tables), and [`match_and_terminal_score`] gives a candidate's **RawScore** against it.
//!   The generating function that turns a RawScore into a SpecEValue lives in `msgf-genfunc`.
//!
//! ## Provenance
//!
//! The scoring half (`model.rs`, `spectrum.rs`, `candidate.rs`, and the constants and rounding
//! conventions in this file) is a **clean-room** implementation: written from a functional
//! specification (`docs/cleanroom/SPEC.md`) distilled from the published MS-GF papers and
//! black-box behaviour, by an implementer who never saw MS-GF+ source or this repository's prior
//! code. It was written in DIA_Proteomics_Rust (`rust/src/dda/specprob/`, commit `b4485bb`, merged
//! as `b930875`; MIT OR Apache-2.0, same author) and brought here unchanged apart from module
//! layout and the [`preprocess`] entry point for the trainer; the node tables, rounding helpers
//! and RawScore buffers were later made faster without changing a bit of output (2026-10-01).
//! See `docs/cleanroom/PROVENANCE.md`.
//! The `.param` reader in `param.rs` was written from `docs/param-format.md`.
//!
//! ```no_run
//! use msgf_scorer::{bundled, match_and_terminal_score, prepare, Candidate, Cleavage, Ms2, Residue};
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let model = bundled::score_model()?;
//! let peaks: Vec<(f64, f64)> = vec![(175.119, 1000.0), (276.155, 800.0)];
//! let prep = prepare(&model, &Ms2 { peaks: &peaks, precursor_mz: 400.2, charge: 2 })
//!     .expect("plausible precursor");
//! let residues: Vec<Residue> =
//!     b"SAMPLER".iter().map(|&l| Residue { letter: l, delta: 0.0 }).collect();
//! let cand = Candidate { residues: &residues, n_term_credit: true };
//! let raw = match_and_terminal_score(&prep, &cand, &Cleavage::trypsin()).expect("standard residues");
//! # let _ = raw; Ok(()) }
//! ```
//!
//! The SpecEValue of `raw` comes from `msgf_genfunc::score_distribution(&prep, &null, Some(raw))`.

pub mod bundled;
mod candidate;
mod cleavage;
mod model;
mod param;
mod spectrum;
mod write;

pub use candidate::{
    cumulative_nominal, match_and_terminal_parts, match_and_terminal_score, Candidate, Residue,
};
pub use cleavage::Cleavage;
pub use model::{IonType, ModelError, Partition as ScorePartition, ScoreModel};
pub use param::{
    read_param, read_param_file, ErrorDist, FragOff, ParamError, Partition, PrecursorOff, RankDist,
    ScoringModel, TERMINATOR,
};
pub use spectrum::{
    peak_by_mass, prepare, prepare_with_cache, preprocess, Ms2, PreparedSpectrum, PreprocessParams,
    RankedPeak,
};
pub use write::{write_param, write_param_file};

// ---- physical constants (Da) -------------------------------------------------------------------

const M_H: f64 = 1.0078250319;
const M_C: f64 = 12.0;
const M_N: f64 = 14.0030740052;
const M_O: f64 = 15.9949146221;
const M_S: f64 = 31.9720707300;
/// Water, as the elemental sum.
pub const WATER: f64 = 2.0 * M_H + M_O;
/// Proton mass used to strip the precursor charge.
pub const PROTON: f64 = 1.0072764669;
/// Charge-carrier mass used by precursor suppression and isotope reduction (used as f32).
const CARRIER: f64 = 1.00727649;
/// First and second isotope spacings of the isotope-cluster reduction (used as f32).
const ISO_STEP1: f64 = 13.00335483 - 12.0;
const ISO_STEP2: f64 = 14.003241 - 13.00335483;
/// Real mass to nominal-grid scaler.
pub const NOMINAL_SCALE: f32 = 0.999497;

/// Elemental formulas `[C, H, N, O, S]` of the 20 residues (free amino acid minus water).
const FORMULAS: [(u8, [u32; 5]); 20] = [
    (b'G', [2, 3, 1, 1, 0]),
    (b'A', [3, 5, 1, 1, 0]),
    (b'S', [3, 5, 1, 2, 0]),
    (b'P', [5, 7, 1, 1, 0]),
    (b'V', [5, 9, 1, 1, 0]),
    (b'T', [4, 7, 1, 2, 0]),
    (b'C', [3, 5, 1, 1, 1]),
    (b'L', [6, 11, 1, 1, 0]),
    (b'I', [6, 11, 1, 1, 0]),
    (b'N', [4, 6, 2, 2, 0]),
    (b'D', [4, 5, 1, 3, 0]),
    (b'Q', [5, 8, 2, 2, 0]),
    (b'K', [6, 12, 2, 1, 0]),
    (b'E', [5, 7, 1, 3, 0]),
    (b'M', [5, 9, 1, 1, 1]),
    (b'H', [6, 7, 3, 1, 0]),
    (b'F', [9, 9, 1, 1, 0]),
    (b'R', [6, 12, 4, 1, 0]),
    (b'Y', [9, 9, 1, 2, 0]),
    (b'W', [11, 10, 2, 1, 0]),
];

/// The 20 residue letters in the conventional alphabet order used by the null model.
pub const RESIDUE_ORDER: &[u8; 20] = b"GASPVTCLINDQKEMHFRYW";

/// Monoisotopic residue mass (binary64 elemental sum), or `None` for a non-standard letter.
#[inline]
pub fn residue_mass(letter: u8) -> Option<f64> {
    static TABLE: std::sync::OnceLock<[Option<f64>; 256]> = std::sync::OnceLock::new();
    TABLE.get_or_init(|| {
        let mut t = [None; 256];
        for l in 0..=255u8 {
            t[l as usize] = residue_mass_from_formula(l);
        }
        t
    })[letter as usize]
}

/// The elemental sum behind [`residue_mass`] (which caches it per letter).
fn residue_mass_from_formula(letter: u8) -> Option<f64> {
    FORMULAS.iter().find(|(l, _)| *l == letter).map(|(_, f)| {
        f[0] as f64 * M_C
            + f[1] as f64 * M_H
            + f[2] as f64 * M_N
            + f[3] as f64 * M_O
            + f[4] as f64 * M_S
    })
}

// ---- rounding conventions (all saturating; NaN -> 0, as Rust `as` casts do) -------------------

/// Score rounding: `floor(x + 0.5)` with the addition in f32 (halves go up).
#[inline]
pub(crate) fn round_score(x: f32) -> i32 {
    floor_to_i32(x + 0.5f32)
}

/// Nominal rounding: nearest, halves away from zero, on an f32 value.
#[inline]
pub(crate) fn round_nominal(x: f32) -> i32 {
    round_to_i32(x)
}

/// `y.floor() as i32` without a libm call (baseline x86-64 has no rounding instruction). Below
/// 2^23 in magnitude `y as i32` truncates exactly and one step down corrects negative
/// non-integers; at or above it every f32 is already an integer (or infinite), and NaN takes the
/// same saturating cast to 0.
#[inline(always)]
fn floor_to_i32(y: f32) -> i32 {
    if y.abs() < 8_388_608.0 {
        let t = y as i32;
        if (t as f32) > y {
            t - 1
        } else {
            t
        }
    } else {
        y as i32
    }
}

/// `x.round() as i32` (halves away from zero) without a libm call. Below 2^23 in magnitude the
/// fractional part `x - trunc(x)` is exact in f32.
#[inline(always)]
fn round_to_i32(x: f32) -> i32 {
    if x.abs() < 8_388_608.0 {
        let t = x as i32;
        let f = x - t as f32;
        if f >= 0.5 {
            t + 1
        } else if f <= -0.5 {
            t - 1
        } else {
            t
        }
    } else {
        x as i32
    }
}

/// Real mass -> nominal integer.
#[inline]
pub fn nominal(m: f64) -> i32 {
    round_nominal(m as f32 * NOMINAL_SCALE)
}

/// Nominal integer -> representative real mass (f32).
#[inline]
pub(crate) fn nominal_to_real(k: i32) -> f32 {
    k as f32 / NOMINAL_SCALE
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rounding_conventions() {
        assert_eq!(round_score(-2.5), -2);
        assert_eq!(round_score(2.5), 3);
        assert_eq!(round_score(f32::NAN), 0);
        assert_eq!(round_score(f32::NEG_INFINITY), i32::MIN);
        assert_eq!(round_nominal(-2.5), -3);
        assert_eq!(round_nominal(2.5), 3);
    }

    /// The libm-free floor and round agree with `f32::floor` / `f32::round` on every f32 class:
    /// a strided sweep over all bit patterns plus the boundary values.
    #[test]
    fn software_rounding_matches_std() {
        let mut check = |x: f32| {
            assert_eq!(
                floor_to_i32(x),
                x.floor() as i32,
                "floor {x:e} ({:#x})",
                x.to_bits()
            );
            assert_eq!(
                round_to_i32(x),
                x.round() as i32,
                "round {x:e} ({:#x})",
                x.to_bits()
            );
            assert_eq!(round_score(x), (x + 0.5f32).floor() as i32, "score {x:e}");
        };
        let mut b = 0u32;
        loop {
            check(f32::from_bits(b));
            match b.checked_add(9_973) {
                Some(n) => b = n,
                None => break,
            }
        }
        for k in -70_000i32..=70_000 {
            let x = k as f32 * 0.25;
            for y in [
                x,
                f32::from_bits(x.to_bits() + 1),
                f32::from_bits(x.to_bits().wrapping_sub(1)),
            ] {
                check(y);
            }
        }
        for x in [
            0.0f32,
            -0.0,
            0.5,
            -0.5,
            0.49999997,
            -0.49999997,
            8_388_607.5,
            -8_388_607.5,
            8_388_608.0,
            -8_388_608.0,
            2_147_483_520.0,
            2_147_483_648.0,
            -2_147_483_648.0,
            -2_147_483_904.0,
            f32::MAX,
            f32::MIN,
            f32::INFINITY,
            f32::NEG_INFINITY,
            f32::NAN,
            f32::MIN_POSITIVE,
            -f32::MIN_POSITIVE,
        ] {
            check(x);
        }
    }

    #[test]
    fn residue_nominal_masses() {
        let noms: Vec<i32> = RESIDUE_ORDER
            .iter()
            .map(|&l| nominal(residue_mass(l).unwrap()))
            .collect();
        assert_eq!(
            noms,
            [
                57, 71, 87, 97, 99, 101, 103, 113, 113, 114, 115, 128, 128, 129, 131, 137, 147,
                156, 163, 186
            ]
        );
        assert!(residue_mass(b'J').is_none());
    }

    /// Callers build candidates with `msgf_chem`'s masses and grid; the scorer must agree with
    /// them bit for bit, or a RawScore would depend on which crate computed a residue mass.
    #[test]
    fn agrees_with_msgf_chem() {
        for &l in RESIDUE_ORDER {
            let a = residue_mass(l).unwrap();
            let b = msgf_chem::residue_mass(l).unwrap();
            assert_eq!(a.to_bits(), b.to_bits(), "residue {}", l as char);
            assert_eq!(nominal(a), msgf_chem::scaling::nominal_bin(a as f32));
        }
        assert_eq!(WATER.to_bits(), msgf_chem::mass::WATER.to_bits());
        assert_eq!(PROTON.to_bits(), msgf_chem::mass::PROTON.to_bits());
        assert_eq!(NOMINAL_SCALE, msgf_chem::scaling::NOMINAL);
    }
}
