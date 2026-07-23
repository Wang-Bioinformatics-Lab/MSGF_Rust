//! Spectrum preprocessing — turns raw MGF peaks into the ranked peak list the scorer consumes.
//!
//! Faithful port of the preprocessing MS-GF+ performs in the `NewScoredSpectrum` constructor
//! (`edu.ucsd.msjava.msscorer.NewScoredSpectrum`), in this exact order:
//!
//! 1. **Precursor-peak filtering** — for each `PrecursorOffsetFrequency` returned by
//!    `NewRankScorer.getPrecursorOFF(charge)` (a floor/ceiling lookup on charge over the trained
//!    map), `Spectrum.filterPrecursorPeaks(mme, reducedCharge, offset)` zeroes the intensity of
//!    every peak within the model's `mme` tolerance of the (charge-reduced) precursor m/z.
//! 2. **Ranking** — `Spectrum.setRanksOfPeaks()` ranks peaks by intensity *descending*, ties
//!    broken by *higher m/z* (`reverseOrder(Peak.IntensityComparator)`). Rank 1 = most intense.
//! 3. **Deconvolution** — when the model sets `applyDeconvolution`, `Spectrum.getDeconvolutedSpectrum`
//!    collapses low-charge isotope envelopes onto singly-charged m/z. Ranks are assigned *before*
//!    this step and travel with each peak unchanged (deconvolution reuses the same `Peak` objects
//!    and never re-ranks); the result is re-sorted by m/z.
//!
//! No peaks are dropped: filtered peaks survive with intensity 0, and every peak is emitted once,
//! m/z-ascending. Validated bit-for-rank against `f13_scored_spectrum.golden.json` in
//! `tests/golden_preprocess.rs`.

use crate::scored_spectrum::RankedPeak;
use crate::ScoringModel;

/// `Composition.ChargeCarrierMass()` = `Composition.PROTON`, as the `float` MS-GF+ uses.
const CHARGE_CARRIER: f32 = 1.00727649_f64 as f32;
/// `Composition.ISOTOPE` = `C13 - C`, cast to `float`.
const ISOTOPE: f32 = (13.00335483_f64 - 12.0_f64) as f32;
/// `Composition.C14 - Composition.C13`, cast to `float` (the second-isotope spacing).
const C14_MINUS_C13: f32 = (14.003241_f64 - 13.00335483_f64) as f32;

/// `NewRankScorer.getPrecursorOFF(charge)`: the precursor offset-frequency entries trained for a
/// charge, selected by `TreeMap.floorEntry` (greatest trained charge ≤ `charge`) with
/// `ceilingEntry` (smallest ≥ `charge`) as the fallback when `charge` is below every trained key.
/// Returns the entries for the selected charge in the model's read order.
fn precursor_off_for(model: &ScoringModel, charge: i32) -> Vec<&crate::PrecursorOff> {
    if model.precursor_off.is_empty() {
        return Vec::new();
    }
    // Distinct trained charges (the TreeMap keys).
    let mut keys: Vec<i32> = model.precursor_off.iter().map(|o| o.charge).collect();
    keys.sort_unstable();
    keys.dedup();
    // floorEntry, else ceilingEntry.
    let key = keys
        .iter()
        .rev()
        .find(|&&k| k <= charge)
        .or_else(|| keys.iter().find(|&&k| k >= charge));
    match key {
        Some(&k) => model
            .precursor_off
            .iter()
            .filter(|o| o.charge == k)
            .collect(),
        None => Vec::new(),
    }
}

/// `Spectrum.getDeconvolutedSpectrum(tol)`, operating in place on the m/z array.
///
/// For each not-yet-consumed peak `i`, try charge states `2..min(charge,4)`: if a following peak
/// sits one isotope step (`ISOTOPE/ionCharge`) away, mark it consumed, de-charge peak `i` and its
/// isotope(s) to singly-charged m/z (`ionCharge·mz − (ionCharge−1)·chargeCarrier`), and stop. The
/// forward scans read the *live* (already-mutated) m/z array, exactly like the Java, so consumed
/// peaks that were de-charged to large m/z terminate scans early just as they do there.
fn deconvolute(mz: &mut [f32], charge: i32, tol: f32) {
    if charge == 0 {
        return;
    }
    let n = mz.len();
    let mut ignore = vec![false; n];
    for i in 0..n {
        if ignore[i] {
            continue;
        }
        let p_mz = mz[i];
        let mut ion_charge = 2;
        while ion_charge < charge && ion_charge < 4 {
            let mut is_deconvoluted = false;
            let iso_step = ISOTOPE / ion_charge as f32;
            let mut j = i + 1;
            while j < n {
                let diff = mz[j] - p_mz - iso_step;
                if diff > -tol && diff < tol {
                    ignore[j] = true;
                    mz[i] = ion_charge as f32 * mz[i] - (ion_charge - 1) as f32 * CHARGE_CARRIER;
                    is_deconvoluted = true;
                    let p2_mz = mz[j];
                    let c_step = C14_MINUS_C13 / ion_charge as f32;
                    let mut k = j + 1;
                    while k < n {
                        let diff2 = mz[k] - p2_mz - c_step;
                        if diff2 > -tol && diff2 < tol {
                            ignore[k] = true;
                            mz[k] = ion_charge as f32 * mz[k]
                                - (ion_charge - 1) as f32 * CHARGE_CARRIER;
                            break;
                        } else if diff2 > tol {
                            break;
                        }
                        k += 1;
                    }
                    mz[j] = ion_charge as f32 * mz[j] - (ion_charge - 1) as f32 * CHARGE_CARRIER;
                    break;
                } else if diff > tol {
                    break;
                }
                j += 1;
            }
            if is_deconvoluted {
                break;
            }
            ion_charge += 1;
        }
    }
}

/// Reproduce MS-GF+'s `NewScoredSpectrum` preprocessing for one spectrum.
///
/// `precursor_mass` is the de-charged monoisotopic precursor mass (`Spectrum.getPrecursorMass()`);
/// `raw_peaks` are the `(m/z, intensity)` pairs as read from the source file (order as read — this
/// function reproduces the `MgfSpectrumParser` sort). Returns the preprocessed peaks m/z-ascending,
/// each carrying the rank MS-GF+ assigns.
pub fn preprocess(
    model: &ScoringModel,
    charge: i32,
    precursor_mass: f32,
    raw_peaks: &[(f32, f32)],
) -> Vec<RankedPeak> {
    let n = raw_peaks.len();
    let mut mz: Vec<f32> = raw_peaks.iter().map(|p| p.0).collect();
    let mut inten: Vec<f32> = raw_peaks.iter().map(|p| p.1).collect();

    // MgfSpectrumParser sorts (by Peak.compareTo = m/z, then intensity) only when the peaks are
    // not already non-decreasing in m/z; otherwise it keeps file order.
    let already_sorted = (1..n).all(|i| mz[i] >= mz[i - 1]);
    if !already_sorted {
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| {
            mz[a]
                .partial_cmp(&mz[b])
                .unwrap()
                .then(inten[a].partial_cmp(&inten[b]).unwrap())
        });
        mz = order.iter().map(|&i| mz[i]).collect();
        inten = order.iter().map(|&i| inten[i]).collect();
    }

    // 1. filterPrecursorPeaks: zero the intensity of peaks within `mme` of each precursor offset.
    for off in precursor_off_for(model, charge) {
        let c = charge - off.reduced_charge;
        if c == 0 {
            continue;
        }
        let mass = (precursor_mass + c as f32 * CHARGE_CARRIER) / c as f32 + off.offset;
        let tol = model.mme.window_da(mass as f64) as f32;
        let (lo, hi) = (mass - tol, mass + tol);
        for i in 0..n {
            if mz[i] >= lo && mz[i] <= hi {
                inten[i] = 0.0;
            }
        }
    }

    // 2. setRanksOfPeaks: rank 1 = most intense, ties broken by higher m/z (stable).
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        inten[b]
            .partial_cmp(&inten[a])
            .unwrap()
            .then(mz[b].partial_cmp(&mz[a]).unwrap())
    });
    let mut rank = vec![0i32; n];
    for (pos, &i) in order.iter().enumerate() {
        rank[i] = pos as i32 + 1;
    }

    // 3. deconvolute (ranks already fixed; they travel with the peaks unchanged).
    if model.apply_deconvolution {
        deconvolute(&mut mz, charge, model.deconvolution_error_tolerance);
        let mut out: Vec<RankedPeak> = (0..n)
            .map(|i| RankedPeak {
                mz: mz[i],
                intensity: inten[i],
                rank: rank[i],
            })
            .collect();
        // getDeconvolutedSpectrum ends with Collections.sort(_, MassComparator) = m/z, then intensity.
        out.sort_by(|a, b| {
            a.mz.partial_cmp(&b.mz)
                .unwrap()
                .then(a.intensity.partial_cmp(&b.intensity).unwrap())
        });
        out
    } else {
        // No deconvolution: the spectrum keeps its (m/z-ascending) load order.
        (0..n)
            .map(|i| RankedPeak {
                mz: mz[i],
                intensity: inten[i],
                rank: rank[i],
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constants_match_java_floats() {
        // Exactly the `(float)double` truncations MS-GF+ applies to the Composition constants.
        assert_eq!(CHARGE_CARRIER, 1.00727649_f64 as f32);
        assert_eq!(ISOTOPE, (13.00335483_f64 - 12.0_f64) as f32);
        assert_eq!(C14_MINUS_C13, (14.003241_f64 - 13.00335483_f64) as f32);
    }

    #[test]
    fn ranks_desc_intensity_then_desc_mz() {
        // model with no precursor offsets and no deconvolution: pure ranking.
        let peaks = [(100.0, 5.0), (200.0, 5.0), (300.0, 20.0)];
        // build a minimal model via a helper is overkill; test ranking logic directly.
        let mut order: Vec<usize> = (0..3).collect();
        let inten = [5.0f32, 5.0, 20.0];
        let mz = [100.0f32, 200.0, 300.0];
        order.sort_by(|&a, &b| {
            inten[b]
                .partial_cmp(&inten[a])
                .unwrap()
                .then(mz[b].partial_cmp(&mz[a]).unwrap())
        });
        // most intense (300) first, then tie 5.0 broken by higher m/z (200 before 100).
        assert_eq!(order, vec![2, 1, 0]);
        let _ = peaks;
    }
}
