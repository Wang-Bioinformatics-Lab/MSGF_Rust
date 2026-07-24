//! The partition scheme — how the corpus is sliced before counting.
//!
//! A `.param` model conditions every trained table on a partition
//! *(precursor charge, parent-mass boundary, mass segment)*. Boundaries are **equal-count
//! quantiles of our own corpus** (target ≈ `min_psms_per_partition` PSMs per bin), not values
//! copied from any existing model — a partition scheme is a property of the training set.
//!
//! [`PartitionScheme::index_of`] reproduces the lookup `msgf-scorer`'s `ScoredSpectrum` performs
//! at scoring time (`TreeSet` floor over `(charge, seg, parent_mass)`, with the charge fallback),
//! so a PSM is counted into exactly the partition that will later score it.

use crate::corpus::TrainingPsm;
use msgf_scorer::Partition;

/// Partition list plus the per-charge boundary table it was built from.
#[derive(Debug, Clone)]
pub struct PartitionScheme {
    /// Sorted by `(charge, seg, parent_mass)` — the order `.param` requires.
    pub partitions: Vec<Partition>,
    pub num_segments: i32,
    /// `(charge, boundaries)` for reporting.
    pub boundaries: Vec<(i32, Vec<f32>)>,
    /// Number of training PSMs whose (charge, mass) fall in each `(charge, boundary)` group,
    /// parallel to `partitions`.
    pub group_spectra: Vec<u64>,
}

impl PartitionScheme {
    /// Build from the corpus: per charge, equal-count parent-mass quantiles.
    ///
    /// A charge with fewer than `min_psms` PSMs gets no partitions of its own; those PSMs are
    /// counted into the nearest trained charge, which is exactly where the scorer's floor lookup
    /// will send such spectra.
    pub fn build(
        psms: &[TrainingPsm],
        num_segments: i32,
        min_psms: usize,
        max_parts_per_charge: usize,
    ) -> Self {
        let mut by_charge: Vec<(i32, Vec<f32>)> = Vec::new();
        for p in psms {
            match by_charge.iter_mut().find(|(c, _)| *c == p.charge) {
                Some((_, v)) => v.push(p.parent_mass),
                None => by_charge.push((p.charge, vec![p.parent_mass])),
            }
        }
        by_charge.sort_by_key(|(c, _)| *c);

        let mut boundaries: Vec<(i32, Vec<f32>)> = Vec::new();
        for (charge, masses) in by_charge.iter_mut() {
            if masses.len() < min_psms {
                continue;
            }
            masses.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let n_parts = (masses.len() / min_psms).clamp(1, max_parts_per_charge);
            let mut b = Vec::with_capacity(n_parts);
            b.push(0.0f32); // the first partition always starts at 0
            for k in 1..n_parts {
                let idx = k * masses.len() / n_parts;
                let v = masses[idx];
                // Strictly increasing boundaries only (a mass tie would collapse a partition).
                if v > *b.last().unwrap() {
                    b.push(v);
                }
            }
            boundaries.push((*charge, b));
        }

        let mut partitions = Vec::new();
        for (charge, b) in &boundaries {
            for seg in 0..num_segments {
                for &m in b {
                    partitions.push(Partition {
                        charge: *charge,
                        parent_mass: m,
                        seg,
                    });
                }
            }
        }
        partitions.sort_by(|a, b| {
            a.charge
                .cmp(&b.charge)
                .then(a.seg.cmp(&b.seg))
                .then(a.parent_mass.partial_cmp(&b.parent_mass).unwrap())
        });

        let mut scheme = Self {
            partitions,
            num_segments,
            boundaries,
            group_spectra: Vec::new(),
        };
        // How many corpus spectra land in each partition (same count for every segment of a
        // (charge, boundary) group — each spectrum contributes sites to all its segments).
        let mut counts = vec![0u64; scheme.partitions.len()];
        for p in psms {
            for seg in 0..num_segments {
                if let Some(i) = scheme.index_of(p.charge, p.parent_mass, seg) {
                    counts[i] += 1;
                }
            }
        }
        scheme.group_spectra = counts;
        scheme
    }

    /// `TreeSet.floor` over `(charge, seg, parent_mass)`.
    fn floor(&self, charge: i32, seg: i32, mass: f32) -> Option<usize> {
        let mut best = None;
        for (i, p) in self.partitions.iter().enumerate() {
            if (p.charge, p.seg, p.parent_mass) <= (charge, seg, mass) {
                best = Some(i);
            }
        }
        best
    }

    /// The partition serving `(charge, parent_mass, seg)` — mirrors `ScoredSpectrum::partition_for`.
    pub fn index_of(&self, charge: i32, parent_mass: f32, seg: i32) -> Option<usize> {
        match self.floor(charge, seg, parent_mass) {
            None => {
                let first_charge = self.partitions.first()?.charge;
                self.floor(first_charge, seg, parent_mass)
            }
            Some(i) => {
                let matched = self.partitions[i].charge;
                if matched == charge {
                    Some(i)
                } else {
                    self.floor(matched, seg, parent_mass)
                }
            }
        }
    }

    /// `NewScoredSpectrum.getSegmentNum` — which mass segment a theoretical m/z falls in.
    #[inline]
    pub fn segment_num(&self, mz: f32, parent_mass: f32) -> i32 {
        let s = (mz / parent_mass * self.num_segments as f32) as i32;
        s.clamp(0, self.num_segments - 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use msgf_chem::peptide::Residue;

    fn psm(charge: i32, mass: f32) -> TrainingPsm {
        TrainingPsm {
            charge,
            parent_mass: mass,
            residues: vec![Residue {
                aa: b'K',
                mod_delta: 0.0,
            }],
            peaks: vec![],
        }
    }

    #[test]
    fn quantile_boundaries_and_lookup() {
        let psms: Vec<_> = (0..40).map(|i| psm(2, 800.0 + i as f32 * 10.0)).collect();
        let s = PartitionScheme::build(&psms, 2, 10, 30);
        // 40 PSMs / 10 per partition = 4 mass bins, × 2 segments
        assert_eq!(s.partitions.len(), 8);
        assert_eq!(s.boundaries[0].1.len(), 4);
        assert_eq!(s.boundaries[0].1[0], 0.0);
        // sorted by (charge, seg, mass)
        assert!(s
            .partitions
            .windows(2)
            .all(|w| (w[0].charge, w[0].seg, w[0].parent_mass)
                < (w[1].charge, w[1].seg, w[1].parent_mass)));
        // a charge with no partitions of its own folds onto the trained charge
        let i = s.index_of(5, 900.0, 1).unwrap();
        assert_eq!(s.partitions[i].charge, 2);
        assert_eq!(s.partitions[i].seg, 1);
    }
}
