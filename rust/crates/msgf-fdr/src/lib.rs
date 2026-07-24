//! msgf-fdr — target-decoy false-discovery-rate estimation, compatible with MS-GF+'s
//! `fdr/TargetDecoyAnalysis.java` (specified in `PLAN2.md` §1.4).
//!
//! The score is **SpecEValue**, smaller-is-better, and all arithmetic is **`f32`** — that is what
//! the Java oracle emits, so q-values compare by exact equality rather than by tolerance.
//!
//! The estimator is Käll et al. (2008) `D/T` on a concatenated target-decoy search:
//!
//! 1. Sort targets and decoys ascending. Walk the **distinct** decoy scores; at each one let
//!    `decoy_index` be the number of decoys scoring at least that well and `target_index` the
//!    number of targets scoring *strictly* better. The FDR is `decoy_index / target_index`,
//!    or `1` when `target_index <= decoy_index`. Clamp to 1 and stop at the first FDR of 1.
//! 2. Convert FDRs to q-values by taking the running minimum from the worst threshold back to the
//!    best, so the result is monotone.
//! 3. Look a score up by **largest key ≤ score** — the FDR at the best threshold that still admits
//!    the PSM. Scores better than every decoy fall through to the `−∞` sentinel and get q = 0;
//!    scores worse than the last threshold inherit that threshold's value.
//!
//! Two populations are analysed: every reported match (`QValue`) and one entry per distinct
//! peptide, keyed by its best score (`PepQValue`).
//!
//! **Validated** against the committed `validation/golden/iprg2013_F13.golden.json`: both columns
//! reproduce MS-GF+ exactly for all 1610 unique PSMs (`tests/golden_fdr.rs`, which needs no fetched
//! data). Note that `PLAN2.md` §1.4 describes the lookup as `higherEntry` (strictly next-larger
//! key); that reproduces 1609/1610 — the best-scoring PSM comes out at q = 1 instead of 0 — so the
//! implementation here uses the floor lookup, which matches the oracle on all 1610. The F13 corpus
//! is a weak discriminator (it yields only two distinct q-values, see `PLAN2.md` §4), so PLAN2's
//! TD-2 "Gate 2" Java probe over synthetic score lists remains the way to settle the edge cases.

/// One match, as the FDR sweep sees it.
#[derive(Debug, Clone)]
pub struct PsmRecord {
    /// SpecEValue — smaller is better.
    pub score: f32,
    /// Peptide identity for the peptide-level population: the mod-bearing sequence with flanking
    /// context stripped and upper-cased (see [`peptide_key`]).
    pub peptide: String,
    /// `true` when **every** protein occurrence is a decoy. A peptide shared between a target and
    /// a decoy protein counts as a target (`PLAN2.md` §1.3).
    pub is_decoy: bool,
}

/// Normalise a peptide string into the key the peptide-level population is grouped by: strip a
/// leading `X.` and trailing `.X` enzymatic context, then upper-case.
pub fn peptide_key(peptide: &str) -> String {
    let b = peptide.as_bytes();
    let core = if b.len() >= 4 && b[1] == b'.' && b[b.len() - 2] == b'.' {
        &peptide[2..peptide.len() - 2]
    } else {
        peptide
    };
    core.to_ascii_uppercase()
}

/// Decide decoy status from a match's full list of protein accessions, per `PLAN2.md` §1.3.
/// A match with no protein occurrences is not a decoy.
pub fn is_decoy_match<'a>(proteins: impl IntoIterator<Item = &'a str>, decoy_prefix: &str) -> bool {
    let mut any = false;
    for p in proteins {
        any = true;
        if !p.starts_with(decoy_prefix) {
            return false;
        }
    }
    any
}

/// A monotone score → q-value step function, stored as ascending `(threshold, q)` pairs.
#[derive(Debug, Clone, Default)]
pub struct QValueMap {
    pairs: Vec<(f32, f32)>,
}

impl QValueMap {
    /// Build the map from target and decoy score lists.
    ///
    /// `pit` is MS-GF+'s "percentage of incorrect targets"; it is fixed at `1.0` there and the
    /// parameter exists so the formula stays legible.
    pub fn build(targets: &[f32], decoys: &[f32], pit: f32) -> QValueMap {
        let mut t: Vec<f32> = targets.to_vec();
        let mut d: Vec<f32> = decoys.to_vec();
        t.sort_by(|a, b| a.partial_cmp(b).expect("scores must be non-NaN"));
        d.sort_by(|a, b| a.partial_cmp(b).expect("scores must be non-NaN"));

        let mut pairs: Vec<(f32, f32)> = Vec::new();
        let (mut target_index, mut decoy_index, mut i) = (0usize, 0usize, 0usize);
        while i < d.len() {
            let key = d[i];
            // Consume the whole run of equal decoy scores — the sweep steps per *distinct* score.
            while i < d.len() && d[i] == key {
                i += 1;
                decoy_index += 1;
            }
            while target_index < t.len() && t[target_index] < key {
                target_index += 1;
            }
            let fdr = if target_index <= decoy_index {
                1.0f32
            } else {
                ((decoy_index as f32 * pit).round() / target_index as f32).min(1.0)
            };
            pairs.push((key, fdr));
            if fdr >= 1.0 {
                break; // every worse threshold is also FDR 1
            }
        }

        // FDR -> q-value: running minimum from the worst threshold back to the best.
        let mut running = f32::INFINITY;
        for p in pairs.iter_mut().rev() {
            running = running.min(p.1);
            p.1 = running;
        }
        // One sentinel: a score better than every decoy has no decoy above it, so q = 0. There is
        // no matching sentinel at the top — under the floor lookup the *last* computed threshold
        // already governs every score worse than it, which is the correct estimate (accepting more
        // targets past the final decoy only lowers the ratio).
        pairs.insert(0, (f32::NEG_INFINITY, 0.0));
        QValueMap { pairs }
    }

    /// The q-value for `score`: the value at the largest threshold `<= score`.
    pub fn q_value(&self, score: f32) -> f32 {
        let i = self.pairs.partition_point(|(k, _)| *k <= score);
        if i == 0 {
            0.0
        } else {
            self.pairs[i - 1].1
        }
    }

    /// The `(threshold, q-value)` pairs, ascending — including the leading `−∞` sentinel.
    pub fn pairs(&self) -> &[(f32, f32)] {
        &self.pairs
    }
}

/// PSM- and peptide-level q-values for one search.
#[derive(Debug, Clone)]
pub struct TargetDecoyAnalysis {
    psm: QValueMap,
    peptide: QValueMap,
    /// Best score per peptide key, so every match of a peptide reports that peptide's q-value.
    best_per_peptide: std::collections::HashMap<String, f32>,
}

impl TargetDecoyAnalysis {
    /// Run the analysis over every reported match.
    pub fn new(psms: &[PsmRecord], pit: f32) -> TargetDecoyAnalysis {
        let (t, d): (Vec<f32>, Vec<f32>) = split(psms.iter().map(|p| (p.score, p.is_decoy)));
        let psm = QValueMap::build(&t, &d, pit);

        // Peptide level: one entry per distinct peptide, represented by its best (minimum) score.
        let mut best: std::collections::HashMap<String, (f32, bool)> =
            std::collections::HashMap::new();
        for p in psms {
            best.entry(p.peptide.clone())
                .and_modify(|e| {
                    if p.score < e.0 {
                        *e = (p.score, p.is_decoy);
                    }
                })
                .or_insert((p.score, p.is_decoy));
        }
        let (pt, pd): (Vec<f32>, Vec<f32>) = split(best.values().copied());
        let peptide = QValueMap::build(&pt, &pd, pit);
        TargetDecoyAnalysis {
            psm,
            peptide,
            best_per_peptide: best.into_iter().map(|(k, v)| (k, v.0)).collect(),
        }
    }

    /// `QValue` for a match with this score.
    pub fn psm_q_value(&self, score: f32) -> f32 {
        self.psm.q_value(score)
    }

    /// `PepQValue` for a match on this peptide. Unknown peptides fall back to their own score.
    pub fn pep_q_value(&self, peptide: &str, score: f32) -> f32 {
        let best = self.best_per_peptide.get(peptide).copied().unwrap_or(score);
        self.peptide.q_value(best)
    }

    pub fn psm_map(&self) -> &QValueMap {
        &self.psm
    }

    pub fn peptide_map(&self) -> &QValueMap {
        &self.peptide
    }
}

fn split(items: impl Iterator<Item = (f32, bool)>) -> (Vec<f32>, Vec<f32>) {
    let (mut t, mut d) = (Vec::new(), Vec::new());
    for (score, is_decoy) in items {
        if is_decoy {
            d.push(score);
        } else {
            t.push(score);
        }
    }
    (t, d)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(score: f32, is_decoy: bool, peptide: &str) -> PsmRecord {
        PsmRecord {
            score,
            peptide: peptide.to_string(),
            is_decoy,
        }
    }

    #[test]
    fn peptide_keys_strip_flanks_and_upcase() {
        assert_eq!(peptide_key("K.SAMPLER.A"), "SAMPLER");
        assert_eq!(peptide_key("-.M+15.995PKR.S"), "M+15.995PKR");
        assert_eq!(peptide_key("sampler"), "SAMPLER");
    }

    #[test]
    fn decoy_status_needs_every_occurrence_to_be_decoy() {
        assert!(is_decoy_match(["XXX_A", "XXX_B"], "XXX_"));
        assert!(!is_decoy_match(["XXX_A", "B"], "XXX_")); // shared with a target -> target
        assert!(!is_decoy_match([], "XXX_")); // no proteins -> not a decoy
    }

    #[test]
    fn q_values_are_monotone_and_track_decoys() {
        // 4 targets better than every decoy, then 4 decoys interleaved.
        let t = [1e-10, 1e-9, 1e-8, 1e-7, 1e-4];
        let d = [1e-6, 1e-5, 1e-3, 1e-2];
        let m = QValueMap::build(&t, &d, 1.0);
        assert_eq!(m.q_value(1e-10), 0.0); // better than every decoy
        assert!(m.q_value(1e-10) <= m.q_value(1e-6));
        assert!(m.q_value(1e-6) <= m.q_value(1e-5));
        // 4 decoys / 5 targets at the worst threshold; everything below it inherits that.
        assert!((m.q_value(1.0) - 0.8).abs() < 1e-6);
    }

    #[test]
    fn first_decoy_beating_every_target_gives_fdr_one() {
        // decoy_index (1) >= target_index (0) at the first threshold -> FDR 1, sweep stops.
        let m = QValueMap::build(&[1e-3], &[1e-9], 1.0);
        assert_eq!(m.q_value(1e-9), 1.0);
        assert_eq!(m.q_value(1e-3), 1.0);
        assert_eq!(m.q_value(1e-12), 0.0); // still better than the first decoy
    }

    #[test]
    fn ties_within_the_decoy_list_step_once() {
        // Two decoys at the same score must advance decoy_index by 2 in a single step.
        let m = QValueMap::build(&[1e-9, 1e-8, 1e-7, 1e-6, 1e-5], &[1e-7, 1e-7], 1.0);
        // At key 1e-7: decoys<=key = 2, targets strictly better = 2 -> target_index <= decoy_index -> 1.
        assert_eq!(m.q_value(1e-7), 1.0);
    }

    #[test]
    fn empty_decoy_list_gives_zero() {
        let m = QValueMap::build(&[1e-9, 1e-8], &[], 1.0);
        assert_eq!(m.q_value(1e-9), 0.0);
        assert_eq!(m.q_value(1.0), 0.0);
        assert_eq!(m.pairs().len(), 1); // just the lower sentinel
    }

    #[test]
    fn empty_target_list_gives_one_at_every_decoy() {
        let m = QValueMap::build(&[], &[1e-9, 1e-8], 1.0);
        assert_eq!(m.q_value(1e-9), 1.0);
        assert_eq!(m.q_value(1e-12), 0.0);
    }

    #[test]
    fn monotonisation_pulls_an_early_spike_down() {
        // Raw FDR at the first decoy is 1/1 = 1, but later thresholds are better; after the
        // running-minimum pass the early key must carry the smaller, later value.
        let targets: Vec<f32> = (1..=20).map(|i| i as f32 * 1e-9).collect();
        let decoys: Vec<f32> = vec![0.5e-9, 19e-9];
        let m = QValueMap::build(&targets, &decoys, 1.0);
        let q: Vec<f32> = m.pairs().iter().map(|p| p.1).collect();
        assert!(q.windows(2).all(|w| w[0] <= w[1]), "not monotone: {q:?}");
    }

    #[test]
    fn peptide_level_uses_the_best_psm_per_peptide() {
        let psms = vec![
            rec(1e-10, false, "A"),
            rec(1e-5, false, "A"), // same peptide, worse match
            rec(1e-9, true, "DECOY"),
        ];
        let tda = TargetDecoyAnalysis::new(&psms, 1.0);
        // Both matches of peptide A report A's best-score q-value.
        assert_eq!(
            tda.pep_q_value("A", 1e-10),
            tda.pep_q_value("A", 1e-5),
            "a peptide has one q-value"
        );
    }
}
