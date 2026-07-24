//! The candidate fragment-ion list the trainer measures frequencies for.
//!
//! Clean-room: these are textbook backbone fragments (b/y series plus the standard small neutral
//! losses and the first two ¹³C isotopologues), enumerated over fragment charge. Which of them
//! end up *in* a model is decided by counting — a candidate is kept only where its observed
//! frequency clears the threshold, per partition.
//!
//! Offsets follow the `.param` convention `mz = residue_mass / charge + offset`
//! (`msgf_scorer::FragOff::mz`), so for a fragment carrying `charge` protons with neutral
//! adjustment `Δ` on a series whose neutral base is `base`:
//!
//! ```text
//! offset = PROTON + (base + Δ) / charge
//! ```
//!
//! with `base = 0` for the prefix (b) series and `base = H₂O` for the suffix (y) series.

use msgf_chem::mass;

/// `¹³C − ¹²C`, the isotope spacing.
pub const ISOTOPE: f64 = 13.003_354_835 - 12.0;
/// Ammonia, `NH₃`.
pub const NH3: f64 = mass::N + 3.0 * mass::H;
/// Carbon monoxide, `CO`.
pub const CO: f64 = mass::C + mass::O;

/// One theoretical ion type under test.
#[derive(Debug, Clone, PartialEq)]
pub struct Candidate {
    pub is_prefix: bool,
    pub charge: i32,
    pub offset: f32,
    /// Human-readable chemistry label (`b`, `y-H2O`, `y+i`, …) — for reports, not the model.
    pub label: String,
    /// The name `read_param` will derive: `{P|S}_{charge}_{round(offset)}`.
    pub name: String,
}

/// Neutral adjustments applied to each series, as `(label, Δ)`.
fn neutral_losses() -> Vec<(&'static str, f64)> {
    vec![
        ("", 0.0),
        ("-H2O", -mass::WATER),
        ("-NH3", -NH3),
        ("-CO", -CO),
        ("+i", ISOTOPE),
        ("+2i", 2.0 * ISOTOPE),
    ]
}

/// Java `Math.round` on the offset, as `read_param` does when deriving the ion name.
fn round_name(x: f32) -> i64 {
    (x + 0.5).floor() as i64
}

/// Every candidate ion type for fragment charges `1..=max_charge`.
pub fn candidates(max_charge: i32) -> Vec<Candidate> {
    let mut out = Vec::new();
    for charge in 1..=max_charge {
        for (is_prefix, series, base) in [(true, "b", 0.0), (false, "y", mass::WATER)] {
            for (loss, delta) in neutral_losses() {
                let offset = (mass::PROTON + (base + delta) / charge as f64) as f32;
                let tag = if is_prefix { 'P' } else { 'S' };
                let z = if charge == 1 {
                    String::new()
                } else {
                    format!("^{charge}")
                };
                out.push(Candidate {
                    is_prefix,
                    charge,
                    offset,
                    label: format!("{series}{loss}{z}"),
                    name: format!("{tag}_{charge}_{}", round_name(offset)),
                });
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The singly-charged candidates must reproduce the ion *names* a `.param` reader derives,
    /// which is what makes a trained model interchangeable with any other in this format.
    #[test]
    fn singly_charged_names_match_the_format() {
        let c = candidates(1);
        let by_label = |l: &str| c.iter().find(|x| x.label == l).unwrap().clone();
        assert_eq!(by_label("b").name, "P_1_1"); // b
        assert_eq!(by_label("y").name, "S_1_19"); // y
        assert_eq!(by_label("b-CO").name, "P_1_-27"); // a
        assert_eq!(by_label("b-H2O").name, "P_1_-17");
        assert_eq!(by_label("b-NH3").name, "P_1_-16");
        assert_eq!(by_label("y-H2O").name, "S_1_1");
        assert_eq!(by_label("y-NH3").name, "S_1_2");
        assert_eq!(by_label("y+i").name, "S_1_20");
        assert_eq!(by_label("y+2i").name, "S_1_21");
        assert_eq!(by_label("b+i").name, "P_1_2");
    }

    #[test]
    fn y_offset_is_water_plus_proton() {
        let y = candidates(1).into_iter().find(|c| c.label == "y").unwrap();
        assert!((y.offset - 19.017_84).abs() < 1e-3);
        let y2 = candidates(2)
            .into_iter()
            .find(|c| c.label == "y^2")
            .unwrap();
        assert!((y2.offset - 10.012_8).abs() < 1e-3);
    }
}
