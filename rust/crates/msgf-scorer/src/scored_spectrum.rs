//! Per-node spectrum scoring — the bridge from a (preprocessed, ranked) peak list to the
//! `prefixScore[nm]` / `suffixScore[nm]` arrays MS-GF+ sums into RawScore.
//!
//! Mirrors `NewScoredSpectrum.getNodeScore(node, isPrefix)` (and `FastScorer`'s use of it): for a
//! nominal node mass, walk the model's ion types per segment, compute each theoretical ion m/z
//! (`FragOff::mz`), look up the matching peak (most intense within the model's mass tolerance),
//! and add the rank-derived `node_score` (or `missing_ion_score` when absent).
//!
//! Spectrum preprocessing (precursor-peak filtering, deconvolution, ranking) is not done here
//! yet — `from_ranked_peaks` takes peaks that already carry ranks, so this layer can be validated
//! against MS-GF+'s own preprocessed peaks independently of the preprocessing port.

use crate::ScoringModel;
use msgf_chem::scaling;

/// A peak with its intensity rank (1 = most intense), as MS-GF+ ranks them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RankedPeak {
    pub mz: f32,
    pub intensity: f32,
    pub rank: i32,
}

/// Assign intensity ranks (1 = highest) to `(mz, intensity)` pairs, per `Spectrum.setRanksOfPeaks`.
/// Ties keep input order (stable sort), matching Java's behaviour closely enough for real data.
pub fn rank_by_intensity(peaks: &[(f32, f32)]) -> Vec<RankedPeak> {
    let mut idx: Vec<usize> = (0..peaks.len()).collect();
    idx.sort_by(|&a, &b| peaks[b].1.partial_cmp(&peaks[a].1).unwrap());
    let mut out = vec![
        RankedPeak {
            mz: 0.0,
            intensity: 0.0,
            rank: 0
        };
        peaks.len()
    ];
    for (rank0, &i) in idx.iter().enumerate() {
        out[i] = RankedPeak {
            mz: peaks[i].0,
            intensity: peaks[i].1,
            rank: rank0 as i32 + 1,
        };
    }
    out
}

/// Most intense peak within `±tol_da` of `mz`, or `None`. Peaks must be sorted by `mz` ascending.
///
/// Mirrors `Spectrum.getPeakByMass` = `Collections.max(window, IntensityComparator)`, whose order
/// is (intensity, then m/z) — so among equal-intensity peaks the **highest m/z** wins. This matters
/// because MS-GF+'s preprocessed spectrum keeps filtered peaks at intensity 0, so windows often
/// contain several zero-intensity peaks and the tie-break decides which rank is scored.
pub fn peak_by_mass(peaks: &[RankedPeak], mz: f32, tol_da: f32) -> Option<&RankedPeak> {
    let (lo, hi) = (mz - tol_da, mz + tol_da);
    let start = peaks.partition_point(|p| p.mz < lo);
    let mut best: Option<&RankedPeak> = None;
    for p in &peaks[start..] {
        if p.mz > hi {
            break;
        }
        // iterating ascending m/z; update on >= so the highest-m/z max-intensity peak wins ties
        match best {
            Some(b) if p.intensity < b.intensity => {}
            _ => best = Some(p),
        }
    }
    best
}

/// A scored spectrum: a model plus a ranked peak list for one precursor.
pub struct ScoredSpectrum<'a> {
    model: &'a ScoringModel,
    parent_mass: f32,
    peaks: Vec<RankedPeak>, // sorted by mz ascending
    /// Partition index serving each segment. Constant per spectrum, so it is precomputed once
    /// rather than re-running the partition `floor` lookup for every nominal mass scored.
    seg_partition: Vec<Option<usize>>,
}

impl<'a> ScoredSpectrum<'a> {
    /// Build from peaks that already carry intensity ranks (e.g. MS-GF+'s preprocessed spectrum).
    pub fn from_ranked_peaks(
        model: &'a ScoringModel,
        charge: i32,
        parent_mass: f32,
        mut peaks: Vec<RankedPeak>,
    ) -> Self {
        peaks.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap());
        let seg_partition = (0..model.num_segments)
            .map(|seg| Self::partition_for(model, charge, parent_mass, seg))
            .collect();
        Self {
            model,
            parent_mass,
            peaks,
            seg_partition,
        }
    }

    fn tol_da(&self, mz: f32) -> f32 {
        self.model.mme.window_da(mz as f64) as f32
    }

    /// `getSegmentNum`: which mass segment a theoretical m/z falls in.
    fn segment_num(&self, mz: f32) -> i32 {
        let s = (mz / self.parent_mass * self.model.num_segments as f32) as i32;
        s.min(self.model.num_segments - 1)
    }

    /// TreeSet `floor` over partitions by (charge, seg, parent_mass): greatest partition ≤ key.
    fn floor(model: &ScoringModel, charge: i32, seg: i32, mass: f32) -> Option<usize> {
        let mut best = None;
        for (i, p) in model.partitions.iter().enumerate() {
            if (p.charge, p.seg, p.parent_mass) <= (charge, seg, mass) {
                best = Some(i);
            }
        }
        best
    }

    /// `getPartition(charge, parentMass, seg)` — the trained partition serving this segment.
    fn partition_for(
        model: &ScoringModel,
        charge: i32,
        parent_mass: f32,
        seg: i32,
    ) -> Option<usize> {
        match Self::floor(model, charge, seg, parent_mass) {
            None => {
                let first_charge = model.partitions.first()?.charge;
                Self::floor(model, first_charge, seg, parent_mass)
            }
            Some(i) => {
                let matched_charge = model.partitions[i].charge;
                if matched_charge == charge {
                    Some(i)
                } else {
                    Self::floor(model, matched_charge, seg, parent_mass)
                }
            }
        }
    }

    /// `NewScoredSpectrum.getNodeScore(node, isPrefix)` for a nominal node mass.
    pub fn node_score(&self, nominal_mass: i32, is_prefix: bool) -> f32 {
        let node_mass = scaling::nominal_to_mass(nominal_mass);
        let mut score = 0.0f32;
        for seg in 0..self.model.num_segments {
            let Some(part_idx) = self.seg_partition[seg as usize] else {
                continue;
            };
            for ion in &self.model.frag_off[part_idx] {
                if ion.is_prefix != is_prefix {
                    continue;
                }
                let theo = ion.mz(node_mass);
                if self.segment_num(theo) != seg {
                    continue;
                }
                score += match peak_by_mass(&self.peaks, theo, self.tol_da(theo)) {
                    Some(p) => self.model.node_score(part_idx, ion, p.rank),
                    None => self.model.missing_ion_score(part_idx, ion),
                };
            }
        }
        score
    }

    /// The `prefixScore[nm]` and `suffixScore[nm]` arrays for `nm` in `0..pep_mass_nominal`
    /// (index 0 = 0.0), as consumed by `FastScorer`.
    pub fn prefix_suffix_scores(&self, pep_mass_nominal: i32) -> (Vec<f32>, Vec<f32>) {
        let n = pep_mass_nominal.max(0) as usize;
        let mut prefix = vec![0.0f32; n];
        let mut suffix = vec![0.0f32; n];
        for nm in 1..n as i32 {
            prefix[nm as usize] = self.node_score(nm, true);
            suffix[nm as usize] = self.node_score(nm, false);
        }
        (prefix, suffix)
    }
}

/// FastScorer node-score summation: `Σ_cleavages round(prefix[pm] + suffix[pepMass − pm])`.
///
/// `nominal_prefix_masses` are the cumulative nominal residue masses at each cleavage, ending
/// with the full peptide nominal mass (last element), per `FastScorer.getScore`. Edge scores
/// (ion-existence + mass-error) are added on top for high-res models — that comes next.
pub fn raw_score_nodes(prefix: &[f32], suffix: &[f32], nominal_prefix_masses: &[i32]) -> i32 {
    let Some(&pep) = nominal_prefix_masses.last() else {
        return 0;
    };
    let mut score = 0i32;
    for &pm in &nominal_prefix_masses[..nominal_prefix_masses.len().saturating_sub(1)] {
        let sm = pep - pm;
        if pm >= 0 && sm >= 0 && (pm as usize) < prefix.len() && (sm as usize) < suffix.len() {
            score += (prefix[pm as usize] + suffix[sm as usize]).round() as i32;
        }
    }
    score
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn raw_score_sums_rounded_node_scores() {
        // pep nominal = 10; cleavages at 3 and 6; prefix/suffix arrays indexed by nominal mass
        let prefix = vec![0.0, 0.0, 0.0, 2.4, 0.0, 0.0, 1.6, 0.0, 0.0, 0.0, 0.0];
        let suffix = vec![0.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0];
        // cleavage 3: round(2.4 + suffix[7]=0.9)=round(3.3)=3; cleavage 6: round(1.6 + suffix[4]=1.5)=round(3.1)=3
        assert_eq!(raw_score_nodes(&prefix, &suffix, &[3, 6, 10]), 6);
    }

    #[test]
    fn ranks_by_intensity_desc() {
        let r = rank_by_intensity(&[(100.0, 5.0), (200.0, 50.0), (300.0, 20.0)]);
        assert_eq!(r[0].rank, 3); // intensity 5 -> lowest -> rank 3
        assert_eq!(r[1].rank, 1); // intensity 50 -> highest -> rank 1
        assert_eq!(r[2].rank, 2);
    }

    #[test]
    fn peak_lookup_picks_most_intense_in_window() {
        let mut peaks =
            rank_by_intensity(&[(99.9, 10.0), (100.05, 40.0), (100.2, 5.0), (500.0, 1.0)]);
        peaks.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap());
        // window ±0.15 around 100.0 -> 99.9, 100.05 in range; most intense is 100.05
        let got = peak_by_mass(&peaks, 100.0, 0.15).unwrap();
        assert!((got.mz - 100.05).abs() < 1e-4);
        // nothing near 300
        assert!(peak_by_mass(&peaks, 300.0, 0.15).is_none());
    }
}
