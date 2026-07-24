//! msgf-fdr — target-decoy false-discovery-rate estimation, compatible with MS-GF+'s
//! `fdr/TargetDecoyAnalysis.java` (specified in `PLAN2.md` §1.4).
//!
//! The score is **SpecEValue**, smaller-is-better, and all arithmetic is **`f32`** — that is what
//! the Java oracle emits, so q-values compare by exact equality rather than by tolerance.
//!
//! The estimator is Käll et al. (2008) `D/T` on a concatenated target-decoy search:
//!
//! 1. Seed the map with both sentinels: `−∞ → 0` and `+∞ → 1`.
//! 2. Sort targets and decoys ascending, then walk the decoy list by index, skipping repeats of a
//!    score already handled. At a decoy score `s` sitting at index `i`, `decoy_index = i` is the
//!    number of decoys scoring **strictly better** than `s` (a run of equal scores is charged at
//!    its first member only) and `target_index` is the number of targets scoring strictly better.
//!    When `target_index == 0` **no entry is written and the sweep continues**; otherwise the FDR
//!    is `1` if `target_index <= decoy_index`, else `round(decoy_index · pit) / target_index`,
//!    clamped to 1. The entry is written, and the sweep stops once an FDR reaches 1.
//!    With no decoys at all, the `+∞` sentinel is rewritten to 0.
//! 3. Convert FDRs to q-values by taking the running minimum from the worst threshold back to the
//!    best, so the result is monotone.
//! 4. Look a score up by the **least key strictly greater than it** (Java `TreeMap.higherEntry`).
//!    A score sitting exactly on a threshold therefore reports the *next* threshold's q-value.
//!
//! Two populations are analysed: every reported match (`QValue`) and one entry per distinct
//! peptide, keyed by its best score (`PepQValue`).
//!
//! **Validated** two ways. Both goldens are MS-GF+-derived, so they are generated locally rather
//! than committed (`validation/golden/README.md`), and each test skips when its golden is absent:
//!
//! - `tests/golden_fdr.rs` (TD-2 Gate 1) — both columns reproduce MS-GF+ exactly for all 1610
//!   unique PSMs of `validation/golden/iprg2013_F13.golden.json` (`reference/generate_golden.sh`).
//! - `tests/golden_fdrmap.rs` (TD-2 Gate 2) — 14 synthetic cases dumped straight out of
//!   `TargetDecoyAnalysis` by `validation/reference/java/DumpFdrMap.java`
//!   (`reference/make_fdr_golden.sh`), comparing every map
//!   entry and every lookup (including each threshold's immediate float neighbours). This is the
//!   gate that pins the three rules F13 cannot see, because it yields only two distinct q-values
//!   (`PLAN2.md` §4): the tie rule in step 2, the `target_index == 0` skip, and the step-4 lookup.

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

        // Both sentinels are seeded before the sweep, exactly as Java does: a score better than
        // every threshold has no decoy above it (q = 0), a score worse than every threshold is
        // indistinguishable from noise (q = 1).
        let mut swept: Vec<(f32, f32)> = Vec::new();
        let mut target_index = 0usize;
        let mut previous: Option<f32> = None;
        for (decoy_index, &key) in d.iter().enumerate() {
            // A run of equal decoy scores is charged once, at its first member — so `decoy_index`
            // is the count of decoys scoring *strictly* better than `key`.
            if previous == Some(key) {
                continue;
            }
            previous = Some(key);
            while target_index < t.len() && t[target_index] < key {
                target_index += 1;
            }
            if target_index == 0 {
                // No target is better than this decoy yet: MS-GF+ writes nothing and keeps going,
                // so a later threshold can still produce a usable ratio.
                continue;
            }
            let fdr = if target_index <= decoy_index {
                1.0f32
            } else {
                // Java `Math.round(float)` is floor(x + 0.5) — `msgf_chem::round_half_up`. The
                // arguments here are non-negative, so this also matches `f32::round`.
                let numerator = (decoy_index as f32 * pit + 0.5).floor();
                (numerator / target_index as f32).min(1.0)
            };
            swept.push((key, fdr));
            if fdr >= 1.0 {
                break; // every worse threshold is also FDR 1
            }
        }

        // With no decoys there is nothing to estimate from, so the upper sentinel drops to 0 and
        // every score comes out at q = 0.
        let top = if decoys.is_empty() { 0.0 } else { 1.0 };
        let mut pairs = Vec::with_capacity(swept.len() + 2);
        pairs.push((f32::NEG_INFINITY, 0.0));
        pairs.extend(swept);
        pairs.push((f32::INFINITY, top));

        // FDR -> q-value: running minimum from the worst threshold back to the best.
        let mut running = 1.0f32;
        for p in pairs.iter_mut().rev() {
            running = running.min(p.1);
            p.1 = running;
        }
        QValueMap { pairs }
    }

    /// The q-value for `score`: the value at the **least threshold strictly greater** than it
    /// (Java `TreeMap.higherEntry`), so a score sitting exactly on a threshold reports the next
    /// one's value.
    ///
    /// Java has no answer for `+∞` — `getPSMQValue` dereferences a null entry — so a SpecEValue of
    /// `+∞` (which cannot arise from the generating function) reports the upper sentinel here.
    pub fn q_value(&self, score: f32) -> f32 {
        let i = self.pairs.partition_point(|(k, _)| *k <= score);
        self.pairs[i.min(self.pairs.len() - 1)].1
    }

    /// The `(threshold, q-value)` pairs, ascending — including both `±∞` sentinels.
    pub fn pairs(&self) -> &[(f32, f32)] {
        &self.pairs
    }
}

/// PSM- and peptide-level q-values for one search.
#[derive(Debug, Clone)]
pub struct TargetDecoyAnalysis {
    psm: QValueMap,
    peptide: QValueMap,
    /// Best score per peptide key among **target** matches, then among **decoy** matches. MS-GF+
    /// keeps the two populations in separate tables and consults the target one first, so a
    /// peptide reported both ways is scored as a target.
    best_target_peptide: std::collections::HashMap<String, f32>,
    best_decoy_peptide: std::collections::HashMap<String, f32>,
}

impl TargetDecoyAnalysis {
    /// Run the analysis over every reported match.
    pub fn new(psms: &[PsmRecord], pit: f32) -> TargetDecoyAnalysis {
        let (t, d): (Vec<f32>, Vec<f32>) = split(psms.iter().map(|p| (p.score, p.is_decoy)));
        let psm = QValueMap::build(&t, &d, pit);

        // Peptide level: one entry per distinct peptide, represented by its best (minimum) score,
        // tabulated separately for the two populations.
        let mut best_target_peptide: std::collections::HashMap<String, f32> = Default::default();
        let mut best_decoy_peptide: std::collections::HashMap<String, f32> = Default::default();
        for p in psms {
            let table = if p.is_decoy {
                &mut best_decoy_peptide
            } else {
                &mut best_target_peptide
            };
            table
                .entry(p.peptide.clone())
                .and_modify(|e| *e = e.min(p.score))
                .or_insert(p.score);
        }
        let pt: Vec<f32> = best_target_peptide.values().copied().collect();
        let pd: Vec<f32> = best_decoy_peptide.values().copied().collect();
        let peptide = QValueMap::build(&pt, &pd, pit);
        TargetDecoyAnalysis {
            psm,
            peptide,
            best_target_peptide,
            best_decoy_peptide,
        }
    }

    /// `QValue` for a match with this score.
    pub fn psm_q_value(&self, score: f32) -> f32 {
        self.psm.q_value(score)
    }

    /// `PepQValue` for a match on this peptide: the q-value of that peptide's best score, taken
    /// from the target table first and the decoy table second. Peptides in neither table (MS-GF+
    /// returns null there) fall back to the score passed in.
    pub fn pep_q_value(&self, peptide: &str, score: f32) -> f32 {
        let best = self
            .best_target_peptide
            .get(peptide)
            .or_else(|| self.best_decoy_peptide.get(peptide))
            .copied()
            .unwrap_or(score);
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

    // The expectations below are hand-derived from the sweep in `TargetDecoyAnalysis.getFDRMap`;
    // `tests/golden_fdrmap.rs` checks the same rules against numbers dumped from the Java itself.

    #[test]
    fn q_values_are_monotone_and_track_decoys() {
        // 4 targets better than every decoy, then 4 decoys interleaved.
        let t = [1e-10, 1e-9, 1e-8, 1e-7, 1e-4];
        let d = [1e-6, 1e-5, 1e-3, 1e-2];
        let m = QValueMap::build(&t, &d, 1.0);
        // Thresholds: 1e-6 -> 0/4, 1e-5 -> 1/4, 1e-3 -> 2/5, 1e-2 -> 3/5.
        assert_eq!(m.q_value(1e-10), 0.0); // better than every decoy
        assert_eq!(m.q_value(1e-6), 0.25); // sits on a threshold -> the next one's value
        assert_eq!(m.q_value(1e-5), 0.4);
        assert_eq!(m.q_value(1.0), 1.0); // worse than every threshold -> upper sentinel
        let q: Vec<f32> = m.pairs().iter().map(|p| p.1).collect();
        assert!(q.windows(2).all(|w| w[0] <= w[1]), "not monotone: {q:?}");
    }

    #[test]
    fn no_target_better_than_any_decoy_leaves_only_sentinels() {
        // target_index stays 0 at the single decoy, so no threshold is ever written and every
        // finite score falls through to the upper sentinel.
        let m = QValueMap::build(&[1e-3], &[1e-9], 1.0);
        assert_eq!(m.pairs().len(), 2);
        assert_eq!(m.q_value(1e-12), 1.0);
        assert_eq!(m.q_value(1e-3), 1.0);
    }

    #[test]
    fn a_run_of_equal_decoys_is_charged_at_its_first_member() {
        // Two decoys at 1e-7 with two targets strictly better: decoy_index is 0 (the run's first
        // index), not 2, so the threshold's FDR is 0/2 and not 1.
        let m = QValueMap::build(&[1e-9, 1e-8, 1e-7, 1e-6, 1e-5], &[1e-7, 1e-7], 1.0);
        assert_eq!(
            m.pairs(),
            [(f32::NEG_INFINITY, 0.0), (1e-7, 0.0), (f32::INFINITY, 1.0)]
        );
        assert_eq!(m.q_value(9.9999994e-8), 0.0); // one ulp better than the threshold
        assert_eq!(m.q_value(1e-7), 1.0); // on the threshold -> the next one's value
    }

    #[test]
    fn a_skipped_threshold_lets_a_later_one_still_count() {
        // No target beats 1e-9, so MS-GF+ writes nothing there and carries on; at 1e-8 all three
        // targets are better, giving 1/3 — an implementation that stopped at the first decoy would
        // report 1 for everything.
        let m = QValueMap::build(&[5e-9, 6e-9, 7e-9], &[1e-9, 1e-8], 1.0);
        assert_eq!(m.pairs().len(), 3);
        assert_eq!(m.q_value(7e-9), 1.0 / 3.0);
        assert_eq!(m.q_value(1e-8), 1.0);
    }

    #[test]
    fn empty_decoy_list_gives_zero() {
        let m = QValueMap::build(&[1e-9, 1e-8], &[], 1.0);
        assert_eq!(m.q_value(1e-9), 0.0);
        assert_eq!(m.q_value(1.0), 0.0);
        // Both sentinels, and the upper one drops to 0 because there is nothing to estimate from.
        assert_eq!(m.pairs(), [(f32::NEG_INFINITY, 0.0), (f32::INFINITY, 0.0)]);
    }

    #[test]
    fn empty_target_list_gives_one_everywhere() {
        let m = QValueMap::build(&[], &[1e-9, 1e-8], 1.0);
        assert_eq!(m.q_value(1e-9), 1.0);
        assert_eq!(m.q_value(1e-12), 1.0);
    }

    #[test]
    fn monotonisation_pulls_an_early_spike_down() {
        // Raw FDR at 1.9e-8 is 1/18; the threshold before it is skipped (no target is better than
        // 0.5e-9), and the running-minimum pass must leave the map non-decreasing regardless.
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
