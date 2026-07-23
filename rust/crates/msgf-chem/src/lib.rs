//! msgf-chem — the mass/chemistry foundation for MSGF_Rust.
//!
//! Monoisotopic atomic masses, amino-acid residue masses (derived from fixed elemental
//! formulas), peptide masses, b/y fragment ions, mass tolerance, and the integer mass
//! discretization used by the generating function (nominal + high-precision scalers, mirroring
//! MS-GF+ `Constants.java`).
//!
//! Ground truth for this crate lives in `validation/golden/chemistry/` and is exercised by
//! `tests/golden_chemistry.rs`. The golden itself is guarded against published peptide calibrants.

pub mod peptide;

/// Monoisotopic masses of the most-abundant isotope (CODATA/AME-2020), plus derived constants.
pub mod mass {
    pub const H: f64 = 1.007_825_031_9;
    pub const C: f64 = 12.0;
    pub const N: f64 = 14.003_074_005_2;
    pub const O: f64 = 15.994_914_622_1;
    pub const S: f64 = 31.972_070_730_0;
    pub const P: f64 = 30.973_761_512_0;

    /// Mass of a proton (used for charging).
    pub const PROTON: f64 = 1.007_276_466_9;
    /// Mass of an electron.
    pub const ELECTRON: f64 = 0.000_548_579_9;
    /// Monoisotopic mass of H2O.
    pub const WATER: f64 = 2.0 * H + O;
}

/// Integer mass discretization for the generating-function DP (see MS-GF+ `Constants.java`).
///
/// `nominal` is the low-res grid (~1 bin/Da); `high_precision` is the high-res grid (~274
/// bins/Da, ~0.0036 Da) — the finer grid is why high-res MSGF is expensive and is the hot path
/// the Rust port optimizes.
// Constants mirror MS-GF+'s `float` literals (Constants.java); keep the full digits for
// traceability even though f32 can't distinguish them all.
#[allow(clippy::excessive_precision)]
pub mod scaling {
    /// `INTEGER_MASS_SCALER` — low-resolution / nominal mass. `f32` to match MS-GF+ exactly.
    pub const NOMINAL: f32 = 0.999_497;
    /// `INTEGER_MASS_SCALER_HIGH_PRECISION` — high-resolution mass.
    pub const HIGH_PRECISION: f32 = 274.335_215;

    /// Real mass → nominal integer, per `NominalMass.toNominalMass`: `round(mass * NOMINAL)`.
    #[inline]
    pub fn nominal_bin(mass: f32) -> i32 {
        (mass * NOMINAL).round() as i32
    }

    /// Nominal integer → representative real mass, per `NominalMass.getMass`: `nominal / NOMINAL`.
    #[inline]
    pub fn nominal_to_mass(nominal: i32) -> f32 {
        nominal as f32 / NOMINAL
    }

    /// Real mass → high-precision (high-res) integer bin.
    #[inline]
    pub fn high_res_bin(mass: f32) -> i32 {
        (mass * HIGH_PRECISION).round() as i32
    }
}

/// Elemental composition `[C, H, N, O, S]` of an amino-acid residue (free AA minus water),
/// or `None` for a non-standard character.
pub fn residue_formula(aa: u8) -> Option<[u16; 5]> {
    // [C, H, N, O, S]
    Some(match aa.to_ascii_uppercase() {
        b'G' => [2, 3, 1, 1, 0],
        b'A' => [3, 5, 1, 1, 0],
        b'S' => [3, 5, 1, 2, 0],
        b'P' => [5, 7, 1, 1, 0],
        b'V' => [5, 9, 1, 1, 0],
        b'T' => [4, 7, 1, 2, 0],
        b'C' => [3, 5, 1, 1, 1],
        b'L' => [6, 11, 1, 1, 0],
        b'I' => [6, 11, 1, 1, 0],
        b'N' => [4, 6, 2, 2, 0],
        b'D' => [4, 5, 1, 3, 0],
        b'Q' => [5, 8, 2, 2, 0],
        b'K' => [6, 12, 2, 1, 0],
        b'E' => [5, 7, 1, 3, 0],
        b'M' => [5, 9, 1, 1, 1],
        b'H' => [6, 7, 3, 1, 0],
        b'F' => [9, 9, 1, 1, 0],
        b'R' => [6, 12, 4, 1, 0],
        b'Y' => [9, 9, 1, 2, 0],
        b'W' => [11, 10, 2, 1, 0],
        _ => return None,
    })
}

/// Monoisotopic residue mass, or `None` for a non-standard character.
#[inline]
pub fn residue_mass(aa: u8) -> Option<f64> {
    residue_formula(aa).map(|[c, h, n, o, s]| {
        c as f64 * mass::C
            + h as f64 * mass::H
            + n as f64 * mass::N
            + o as f64 * mass::O
            + s as f64 * mass::S
    })
}

/// Neutral monoisotopic mass of a bare peptide sequence (Σ residues + H2O).
/// Returns `None` if any character is a non-standard residue.
pub fn peptide_neutral_mass(seq: &str) -> Option<f64> {
    let mut m = mass::WATER;
    for b in seq.bytes() {
        m += residue_mass(b)?;
    }
    Some(m)
}

/// m/z of a neutral mass at the given charge: `(neutral + z·proton) / z`.
#[inline]
pub fn mz(neutral: f64, charge: u32) -> f64 {
    debug_assert!(charge >= 1);
    (neutral + charge as f64 * mass::PROTON) / charge as f64
}

/// A singly/doubly-charged fragment ion at a cleavage position.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FragmentIon {
    /// Cleavage index (b_i / y_j; 1-based, counting residues from the fragment's own terminus).
    pub index: usize,
    /// m/z at charge 1.
    pub mz1: f64,
    /// m/z at charge 2.
    pub mz2: f64,
}

/// b-ion series (N-terminal). `b` neutral mass = Σ prefix residues.
pub fn b_ions(seq: &str) -> Vec<FragmentIon> {
    let bytes = seq.as_bytes();
    let mut out = Vec::with_capacity(bytes.len().saturating_sub(1));
    let mut prefix = 0.0f64;
    for i in 1..bytes.len() {
        prefix += residue_mass(bytes[i - 1]).expect("standard residue");
        out.push(FragmentIon {
            index: i,
            mz1: mz(prefix, 1),
            mz2: mz(prefix, 2),
        });
    }
    out
}

/// y-ion series (C-terminal). `y` neutral mass = Σ suffix residues + H2O.
pub fn y_ions(seq: &str) -> Vec<FragmentIon> {
    let bytes = seq.as_bytes();
    let n = bytes.len();
    let mut out = Vec::with_capacity(n.saturating_sub(1));
    let mut suffix = 0.0f64;
    for j in 1..n {
        suffix += residue_mass(bytes[n - j]).expect("standard residue");
        out.push(FragmentIon {
            index: j,
            mz1: mz(suffix + mass::WATER, 1),
            mz2: mz(suffix + mass::WATER, 2),
        });
    }
    out
}

/// Mass tolerance unit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Unit {
    /// Absolute Daltons.
    Da,
    /// Parts-per-million, relative to the theoretical mass.
    Ppm,
}

/// A symmetric mass tolerance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tolerance {
    pub value: f64,
    pub unit: Unit,
}

impl Tolerance {
    pub fn da(value: f64) -> Self {
        Self {
            value,
            unit: Unit::Da,
        }
    }
    pub fn ppm(value: f64) -> Self {
        Self {
            value,
            unit: Unit::Ppm,
        }
    }
    /// The tolerance window in Daltons around `theoretical`.
    #[inline]
    pub fn window_da(&self, theoretical: f64) -> f64 {
        match self.unit {
            Unit::Da => self.value,
            Unit::Ppm => theoretical * self.value * 1e-6,
        }
    }
    /// Whether `observed` is within tolerance of `theoretical`.
    #[inline]
    pub fn matches(&self, observed: f64, theoretical: f64) -> bool {
        (observed - theoretical).abs() <= self.window_da(theoretical)
    }
    /// `[low, high]` bounds around `theoretical`.
    #[inline]
    pub fn bounds(&self, theoretical: f64) -> (f64, f64) {
        let d = self.window_da(theoretical);
        (theoretical - d, theoretical + d)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() <= tol, "{a} vs {b} (Δ={})", (a - b).abs());
    }

    #[test]
    fn glycine_and_water() {
        approx(residue_mass(b'G').unwrap(), 57.021_463_7, 1e-6);
        approx(mass::WATER, 18.010_564_7, 1e-6);
    }

    #[test]
    fn mrfa_calibrant() {
        // Met-Arg-Phe-Ala, published monoisotopic [M+H]+ = 524.2649
        let m = peptide_neutral_mass("MRFA").unwrap();
        approx(mz(m, 1), 524.2649, 3e-3);
    }

    #[test]
    fn tolerance_ppm() {
        let t = Tolerance::ppm(10.0);
        approx(t.window_da(1000.0), 0.01, 1e-12);
        assert!(t.matches(1000.005, 1000.0));
        assert!(!t.matches(1000.02, 1000.0));
    }

    #[test]
    fn unknown_residue() {
        assert!(residue_mass(b'X').is_none());
        assert!(peptide_neutral_mass("PEPXIDE").is_none());
    }

    #[test]
    fn high_res_grid_is_finer() {
        // ~274 bins per Da vs ~1 bin per Da
        assert_eq!(scaling::nominal_bin(1000.0), 999);
        assert!(scaling::high_res_bin(1000.0) > 274_000);
    }

    #[test]
    fn nominal_mass_conversion() {
        // matches Java NominalMass: toNominalMass then getMass round-trips closely
        let nm = scaling::nominal_bin(1234.567);
        approx(scaling::nominal_to_mass(nm) as f64, 1234.567, 1.0);
        assert_eq!(scaling::nominal_bin(57.02146), 57); // glycine residue -> nominal 57
    }
}
