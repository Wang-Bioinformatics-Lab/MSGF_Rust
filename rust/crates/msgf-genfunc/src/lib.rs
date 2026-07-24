//! msgf-genfunc — MS-GF+'s generating function: the spectral probability (SpecEValue, the
//! p-value) of a peptide-spectrum match.
//!
//! The score distribution over *all* peptides of the precursor mass is computed by a dynamic
//! program over the de novo graph (nodes = nominal masses, edges = amino-acid transitions). Each
//! node carries a node score and each edge an edge score (both from the validated `msgf-scorer`
//! model) and an amino-acid probability. `ScoreDist` is the per-node distribution; the DP is a
//! shifted, probability-weighted convolution (mirroring `GeneratingFunction` + `ScoreDist` in
//! MS-GF+). The p-value is the upper tail of the final distribution at the observed RawScore.
//!
//! This module holds the graph-agnostic core (ScoreDist + the DP). The graph builder that turns a
//! `ScoredSpectrum` into nodes/edges (reverse graph, sink tolerance, cleavage credits) is layered
//! on top and validated against MS-GF+'s own DeNovoScore / SpecEValue.

pub mod graph;

/// A score distribution: `probs[i]` is the probability of score `min_score + i`. Mirrors
/// `edu.ucsd.msjava.msgf.ScoreDist`.
#[derive(Debug, Clone)]
pub struct ScoreDist {
    pub min_score: i32,
    pub probs: Vec<f64>,
}

impl ScoreDist {
    /// Empty distribution spanning scores `[min_score, max_score)`.
    pub fn new(min_score: i32, max_score: i32) -> Self {
        Self {
            min_score,
            probs: vec![0.0; (max_score - min_score).max(0) as usize],
        }
    }

    /// A point mass: probability `prob` at `score`.
    pub fn point(score: i32, prob: f64) -> Self {
        let mut d = Self::new(score, score + 1);
        if !d.probs.is_empty() {
            d.probs[0] = prob;
        }
        d
    }

    /// Exclusive upper score bound.
    #[inline]
    pub fn max_score(&self) -> i32 {
        self.min_score + self.probs.len() as i32
    }

    /// Add `other`, shifted by `score_diff` and scaled by `aa_prob`, into this distribution.
    /// Mirrors `ScoreDist.addProbDist`.
    pub fn add_prob_dist(&mut self, other: &ScoreDist, score_diff: i32, aa_prob: f64) {
        let lo = other.min_score.max(self.min_score - score_diff);
        for t in lo..other.max_score() {
            let src = (t - other.min_score) as usize;
            let dst = (t + score_diff - self.min_score) as usize;
            self.probs[dst] += other.probs[src] * aa_prob;
        }
    }

    /// Spectral probability = `P(score >= threshold)`, the upper tail. Mirrors
    /// `ScoreDist.getSpectralProbability(int)` (capped at 1).
    pub fn spectral_probability(&self, threshold: i32) -> f64 {
        let start = if threshold >= self.min_score {
            (threshold - self.min_score) as usize
        } else {
            0
        };
        let p: f64 = self.probs[start.min(self.probs.len())..].iter().sum();
        p.min(1.0)
    }
}

/// One incoming edge to a node: `prev` (predecessor index), its edge score, and the amino-acid
/// probability weighting it.
#[derive(Debug, Clone)]
pub struct Edge {
    pub prev: usize,
    pub edge_score: i32,
    pub aa_prob: f64,
}

/// A graph node with its node score and incoming edges. Nodes must be in topological (increasing
/// nominal mass) order; node 0 is the source.
#[derive(Debug, Clone, Default)]
pub struct Node {
    pub node_score: i32,
    pub edges: Vec<Edge>,
}

/// Neighboring-amino-acid cleavage weighting applied to the final distribution.
#[derive(Debug, Clone, Copy)]
pub struct Cleavage {
    pub credit: i32,
    pub penalty: i32,
    pub prob_cleavage_sites: f64,
}

/// The computed generating function.
pub struct GenFunc {
    pub dist: ScoreDist,
}

impl GenFunc {
    /// Spectral probability (the p-value) at the observed RawScore.
    pub fn spectral_probability(&self, raw_score: i32) -> f64 {
        self.dist.spectral_probability(raw_score)
    }
    /// Maximum achievable score (DeNovoScore) is `max_score() - 1`.
    pub fn max_score(&self) -> i32 {
        self.dist.max_score() - 1
    }
}

/// Merge a `GeneratingFunctionGroup`: sum the per-graph distributions (one graph per candidate
/// peptide mass in the isotope/precursor-tolerance range). Mirrors
/// `GeneratingFunctionGroup.computeGeneratingFunction`.
pub fn merge_group(gfs: &[GenFunc]) -> Option<GenFunc> {
    let min = gfs.iter().map(|g| g.dist.min_score).min()?;
    let max = gfs.iter().map(|g| g.dist.max_score()).max()?;
    if max <= min {
        return None;
    }
    let mut merged = ScoreDist::new(min, max);
    for g in gfs {
        merged.add_prob_dist(&g.dist, 0, 1.0);
    }
    Some(GenFunc { dist: merged })
}

/// Compute the generating function over `nodes` (topological order, node 0 = source), summing the
/// distributions of `sinks` and applying the neighboring-AA `cleavage` weighting. Mirrors
/// `GeneratingFunction.computeGeneratingFunction`. Returns `None` if the sinks are unreachable.
pub fn compute(nodes: &[Node], sinks: &[usize], cleavage: Option<Cleavage>) -> Option<GenFunc> {
    let mut fwd: Vec<Option<ScoreDist>> = vec![None; nodes.len()];
    fwd[0] = Some(ScoreDist::point(0, 1.0));

    for i in 1..nodes.len() {
        let node_score = nodes[i].node_score;
        let (mut cur_min, mut cur_max) = (i32::MAX, i32::MIN);
        for e in &nodes[i].edges {
            if let Some(prev) = &fwd[e.prev] {
                let combined = node_score + e.edge_score;
                cur_max = cur_max.max(prev.max_score() + combined);
                cur_min = cur_min.min(prev.min_score + combined);
            }
        }
        if cur_min >= cur_max {
            continue; // unreachable
        }
        let mut cur = ScoreDist::new(cur_min, cur_max);
        for e in &nodes[i].edges {
            if let Some(prev) = &fwd[e.prev] {
                cur.add_prob_dist(prev, node_score + e.edge_score, e.aa_prob);
            }
        }
        fwd[i] = Some(cur);
    }

    // merge sink distributions
    let (mut min, mut max) = (i32::MAX, i32::MIN);
    for &s in sinks {
        if let Some(d) = &fwd[s] {
            min = min.min(d.min_score);
            max = max.max(d.max_score());
        }
    }
    if max <= min {
        return None;
    }
    let mut merged = ScoreDist::new(min, max);
    for &s in sinks {
        if let Some(d) = &fwd[s] {
            merged.add_prob_dist(d, 0, 1.0);
        }
    }

    // neighboring amino-acid cleavage credit/penalty (probabilistic)
    let dist = match cleavage {
        Some(w) => {
            let mut f = ScoreDist::new(merged.min_score + w.penalty, merged.max_score() + w.credit);
            f.add_prob_dist(&merged, w.credit, w.prob_cleavage_sites);
            f.add_prob_dist(&merged, w.penalty, 1.0 - w.prob_cleavage_sites);
            f
        }
        None => merged,
    };
    Some(GenFunc { dist })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-12, "{a} vs {b}");
    }

    #[test]
    fn score_dist_shift_and_tail() {
        let mut d = ScoreDist::new(3, 6); // scores 3,4,5
        let other = ScoreDist::point(0, 1.0);
        d.add_prob_dist(&other, 4, 0.5); // put 0.5 at score 4
        approx(d.spectral_probability(4), 0.5);
        approx(d.spectral_probability(5), 0.0);
        approx(d.spectral_probability(0), 0.5); // below min → whole tail
    }

    #[test]
    fn tiny_generating_function() {
        // source(0) --aa(0.5), edgeScore 1--> node1(nodeScore 2) --aa(0.5), edgeScore 3--> sink2(nodeScore 0)
        let nodes = vec![
            Node::default(), // source
            Node {
                node_score: 2,
                edges: vec![Edge {
                    prev: 0,
                    edge_score: 1,
                    aa_prob: 0.5,
                }],
            },
            Node {
                node_score: 0,
                edges: vec![Edge {
                    prev: 1,
                    edge_score: 3,
                    aa_prob: 0.5,
                }],
            },
        ];
        let gf = compute(&nodes, &[2], None).unwrap();
        // path score = (2+1) + (0+3) = 6 with prob 0.5*0.5 = 0.25
        approx(gf.spectral_probability(6), 0.25);
        approx(gf.spectral_probability(7), 0.0);
        assert_eq!(gf.max_score(), 6); // max achievable score
    }

    #[test]
    fn neighboring_cleavage_splits_mass() {
        let nodes = vec![
            Node::default(),
            Node {
                node_score: 0,
                edges: vec![Edge {
                    prev: 0,
                    edge_score: 0,
                    aa_prob: 1.0,
                }],
            },
        ];
        // merged: score 0 prob 1. cleavage credit +2 (p=0.25), penalty -1 (p=0.75)
        let gf = compute(
            &nodes,
            &[1],
            Some(Cleavage {
                credit: 2,
                penalty: -1,
                prob_cleavage_sites: 0.25,
            }),
        )
        .unwrap();
        approx(gf.spectral_probability(2), 0.25); // score 2 with prob 0.25
        approx(gf.spectral_probability(-1), 1.0); // full mass
    }
}
