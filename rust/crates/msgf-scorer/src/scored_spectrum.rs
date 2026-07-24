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

use crate::{FragOff, ScoringModel};
use msgf_chem::{mass, scaling};
use std::collections::HashMap;

/// `FlexAminoAcidGraph.MODIFIED_EDGE_PENALTY` (0 in current MS-GF+).
const MODIFIED_EDGE_PENALTY: i32 = 0;

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

/// Build the integer-m/z bucket index: `bucket[b]` = index of the first peak with `mz >= b`.
/// Buckets past the last peak point one-past-the-end. Peaks must be sorted by `mz` ascending.
fn build_peak_bucket(peaks: &[RankedPeak]) -> Vec<u32> {
    let max_mz = peaks.last().map(|p| p.mz).unwrap_or(0.0);
    let nb = (max_mz.max(0.0).floor() as usize) + 2;
    let mut bucket = vec![peaks.len() as u32; nb];
    let mut bi = 0usize;
    for (i, p) in peaks.iter().enumerate() {
        let b = p.mz.max(0.0).floor() as usize;
        while bi <= b {
            bucket[bi] = i as u32;
            bi += 1;
        }
    }
    bucket
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
    /// `peak_bucket[b]` = index of the first peak with `mz >= b` (integer m/z). Lets
    /// `peak_by_mass_idx` start the window scan in O(1) instead of a per-lookup binary search —
    /// the scoring node scans thousands of theoretical masses per spectrum.
    peak_bucket: Vec<u32>,
    /// Partition index serving each segment. Constant per spectrum, so it is precomputed once
    /// rather than re-running the partition `floor` lookup for every nominal mass scored.
    seg_partition: Vec<Option<usize>>,
    /// The dominant ion type (`mainIon`) for this precursor — used for `getNodeMass`.
    main_ion: Option<FragOff>,
    /// `getIonExistenceScore` precomputed for the 4 (cur/prev-present) indices — a per-spectrum
    /// constant, so its `ln` is evaluated once here rather than on every edge. `None` when there
    /// is no edge partition (edge scoring would not be valid, matching the old `.expect`).
    ion_existence_cache: Option<[f32; 4]>,
    /// `getErrorScore` precomputed for every quantized error bin `ei in 0..2·esf+1`.
    error_score_cache: Vec<f32>,
}

/// Per-node quantities shared across a spectrum's candidate graphs (see
/// [`ScoredSpectrum::tables`]). Each vector is indexed by nominal mass `0..=max_nominal`.
#[derive(Debug, Clone)]
pub struct SpectrumTables {
    /// `getNodeMass(k)` — the main ion's matched-peak mass, or `-1` when absent.
    pub node_mass: Vec<f32>,
    /// `getNodeScore(k, isPrefix = true)`.
    pub prefix: Vec<f32>,
    /// `getNodeScore(k, isPrefix = false)`.
    pub suffix: Vec<f32>,
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
        let peak_bucket = build_peak_bucket(&peaks);
        let seg_partition: Vec<Option<usize>> = (0..model.num_segments)
            .map(|seg| Self::partition_for(model, charge, parent_mass, seg))
            .collect();
        let edge_partition = *seg_partition.last().unwrap_or(&None);
        let main_ion = Self::compute_main_ion(model, &seg_partition);
        // probPeak = |peaks| / max(peptideMass / (2·mme), 1), per NewScoredSpectrum
        let peptide_mass = parent_mass - mass::WATER as f32;
        let approx_bins = (peptide_mass / (model.mme.value as f32 * 2.0)).max(1.0);
        let prob_peak = (peaks.len().max(1) as f32) / approx_bins;

        // Precompute the two edge-score lookups. Both depend only on the (per-spectrum) edge
        // partition and probPeak, so the `ln` in each is evaluated once here — not per edge. The
        // expressions mirror `ion_existence_score` / `error_score` exactly, keeping values bit-equal.
        let (ion_existence_cache, error_score_cache) = match edge_partition {
            Some(part) if part < model.error_dist.len() => {
                let ed = &model.error_dist[part];
                let mut ie4 = [0.0f32; 4];
                for (index, slot) in ie4.iter_mut().enumerate() {
                    let noise = match index {
                        0 => (1.0 - prob_peak) * (1.0 - prob_peak),
                        3 => prob_peak * prob_peak,
                        _ => prob_peak * (1.0 - prob_peak),
                    };
                    let ip = if ed.ion_existence[index] == 0.0 {
                        0.01
                    } else {
                        ed.ion_existence[index]
                    };
                    *slot = (ip as f64 / noise as f64).ln() as f32;
                }
                let est: Vec<f32> = ed
                    .signal
                    .iter()
                    .zip(&ed.noise)
                    .map(|(&s, &n)| (s as f64 / n as f64).ln() as f32)
                    .collect();
                (Some(ie4), est)
            }
            _ => (None, Vec::new()),
        };

        Self {
            model,
            parent_mass,
            peaks,
            peak_bucket,
            seg_partition,
            main_ion,
            ion_existence_cache,
            error_score_cache,
        }
    }

    /// `determineIonTypes` main ion: sum fragment frequencies across the segment partitions of
    /// this precursor's (charge, parent_mass) group, keyed by (name, exact offset), and pick the max.
    fn compute_main_ion(model: &ScoringModel, seg_partition: &[Option<usize>]) -> Option<FragOff> {
        let mut acc: HashMap<(&str, u32), (f32, &FragOff)> = HashMap::new();
        for &pi in seg_partition.iter().flatten() {
            for fo in &model.frag_off[pi] {
                let e = acc
                    .entry((fo.name.as_str(), fo.offset.to_bits()))
                    .or_insert((0.0, fo));
                e.0 += fo.frequency;
            }
        }
        acc.into_values()
            .max_by(|a, b| a.0.partial_cmp(&b.0).unwrap())
            .map(|(_, fo)| fo.clone())
    }

    fn tol_da(&self, mz: f32) -> f32 {
        self.model.mme.window_da(mz as f64) as f32
    }

    /// [`peak_by_mass`] using the `peak_bucket` index for the window start — returns the identical
    /// peak (same `±tol_da` window, same intensity/highest-m/z tiebreak), just without the binary
    /// search. Bit-for-bit equivalent, so all scoring goldens are unaffected.
    fn peak_by_mass_idx(&self, mz: f32, tol_da: f32) -> Option<&RankedPeak> {
        let lo = mz - tol_da;
        let hi = mz + tol_da;
        let b = lo.max(0.0).floor() as usize;
        let mut start = self
            .peak_bucket
            .get(b)
            .copied()
            .unwrap_or(self.peaks.len() as u32) as usize;
        // The bucket gives the first peak with `mz >= b` (b ≤ lo); step to the first `mz >= lo`.
        while start < self.peaks.len() && self.peaks[start].mz < lo {
            start += 1;
        }
        let mut best: Option<&RankedPeak> = None;
        for p in &self.peaks[start..] {
            if p.mz > hi {
                break;
            }
            match best {
                Some(bp) if p.intensity < bp.intensity => {}
                _ => best = Some(p),
            }
        }
        best
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
                score += match self.peak_by_mass_idx(theo, self.tol_da(theo)) {
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

    // ---- edge scoring: high-res RawScore = node scores + edge scores (DBScanScorer) ----

    /// Whether the main ion is an N-terminal (prefix) ion; drives the edge summation direction.
    pub fn main_ion_is_prefix(&self) -> bool {
        self.main_ion.as_ref().map(|i| i.is_prefix).unwrap_or(false)
    }

    /// `getNodeMass` for every nominal mass `0..=max_nominal`, precomputed once. Edge scoring
    /// resolves each node's main-ion mass with a single peak lookup here instead of one per
    /// incident edge (~19 amino acids). Depends only on the spectrum — not on the candidate
    /// peptide mass — so the table is identical across the isotope-error candidate graphs.
    pub fn node_masses(&self, max_nominal: i32) -> Vec<f32> {
        (0..=max_nominal.max(0))
            .map(|k| self.node_mass(k))
            .collect()
    }

    /// Precompute every per-node quantity the graph builder needs, for all nominal masses
    /// `0..=max_nominal`: the main-ion mass and the prefix/suffix node scores. **None of these
    /// depend on the candidate peptide mass**, so a spectrum's isotope-error candidate graphs (e.g.
    /// `-ti 0,1`) share one instance instead of recomputing the peak lookups per graph. Build this
    /// once per spectrum with `max_nominal` = the largest candidate mass, then pass it to
    /// `build_reverse_graph`.
    pub fn tables(&self, max_nominal: i32) -> SpectrumTables {
        let m = max_nominal.max(0);
        let node_mass = self.node_masses(m);
        let mut prefix = vec![0.0f32; (m + 1) as usize];
        let mut suffix = vec![0.0f32; (m + 1) as usize];
        for k in 1..=m {
            prefix[k as usize] = self.node_score(k, true);
            suffix[k as usize] = self.node_score(k, false);
        }
        SpectrumTables {
            node_mass,
            prefix,
            suffix,
        }
    }

    /// `getNodeMass`: the main ion's matched-peak mass for a node, or -1 when absent.
    fn node_mass(&self, nominal: i32) -> f32 {
        if nominal == 0 {
            return 0.0;
        }
        let Some(main_ion) = &self.main_ion else {
            return -1.0;
        };
        let theo = main_ion.mz(scaling::nominal_to_mass(nominal));
        match self.peak_by_mass_idx(theo, self.tol_da(theo)) {
            Some(p) => main_ion.mass(p.mz),
            None => -1.0,
        }
    }

    /// `getIonExistenceScore(partition, index, probPeak)` — a per-spectrum table lookup (built in
    /// `from_ranked_peaks`).
    fn ion_existence_score(&self, index: usize) -> f32 {
        self.ion_existence_cache.expect("edge partition present")[index]
    }

    /// `getErrorScore(partition, error)` — quantize the mass error to its bin, then read the
    /// per-spectrum table (built in `from_ranked_peaks`).
    fn error_score(&self, error: f32) -> f32 {
        let esf = self.model.error_scaling_factor;
        let ei = msgf_chem::round_half_up(error * esf as f32).clamp(-esf, esf) + esf;
        self.error_score_cache[ei as usize]
    }

    /// `DBScanScorer.getEdgeScoreInt` for one edge (between two nominal node masses).
    pub fn edge_score(&self, cur: i32, prev: i32, theo_mass: f32, max_nominal: i32) -> i32 {
        if cur >= max_nominal || prev >= max_nominal || cur < 0 || prev < 0 {
            return 0;
        }
        self.edge_score_with(self.node_mass(cur), self.node_mass(prev), theo_mass)
    }

    /// `getEdgeScoreInt` given the two nodes' already-resolved main-ion masses (from
    /// [`node_masses`](Self::node_masses)). Splitting the mass resolution out lets the graph builder
    /// look each node mass up once instead of once per incident edge; the arithmetic is identical.
    pub fn edge_score_with(&self, cur_mass: f32, prev_mass: f32, theo_mass: f32) -> i32 {
        let mut index = 0usize;
        if cur_mass >= 0.0 {
            index += 1;
        }
        if prev_mass >= 0.0 {
            index += 2;
        }
        let mut edge = self.ion_existence_score(index);
        if index == 3 {
            edge += self.error_score(cur_mass - prev_mass - theo_mass);
        }
        msgf_chem::round_half_up(edge)
    }

    /// Full RawScore (node + edge) for a peptide, mirroring `DBScanScorer.getScore`. `nominal_prefix`
    /// / `accurate_prefix` are the cumulative prefix masses (no leading zero; last = full peptide).
    /// The trypsin cleavage credit `DBScanner` adds on top is applied by the caller.
    pub fn raw_score(&self, nominal_prefix: &[i32], accurate_prefix: &[f64], num_mods: i32) -> i32 {
        if nominal_prefix.len() < 2 {
            return 0;
        }
        // MS-GF+ uses a leading-zero convention with fromIndex=1; mirror it.
        let mut nominal = Vec::with_capacity(nominal_prefix.len() + 1);
        nominal.push(0);
        nominal.extend_from_slice(nominal_prefix);
        let mut accurate = Vec::with_capacity(accurate_prefix.len() + 1);
        accurate.push(0.0);
        accurate.extend_from_slice(accurate_prefix);
        let (from, to) = (1usize, nominal.len());
        let pep = nominal[to - 1];

        // node scores (FastScorer): Σ round(prefix[pm] + suffix[pep − pm])
        let mut score = 0i32;
        for &pm in &nominal[from..to - 1] {
            let sm = pep - pm;
            score +=
                msgf_chem::round_half_up(self.node_score(pm, true) + self.node_score(sm, false));
        }
        score += MODIFIED_EDGE_PENALTY * num_mods;

        // edge scores (DBScanScorer), direction set by the main ion
        let max_n = pep + 1;
        if !self.main_ion_is_prefix() {
            for i in (from..=to - 2).rev() {
                let theo = (accurate[i + 1] - accurate[i]) as f32;
                score += self.edge_score(pep - nominal[i], pep - nominal[i + 1], theo, max_n);
            }
        } else {
            for i in from..=to - 2 {
                let theo = (accurate[i] - accurate[i - 1]) as f32;
                score += self.edge_score(nominal[i], nominal[i - 1], theo, max_n);
            }
        }
        score
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
            score += msgf_chem::round_half_up(prefix[pm as usize] + suffix[sm as usize]);
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
