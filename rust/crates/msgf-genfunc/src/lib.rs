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
//! The graph is stored in flat CSR form (`Graph`) and the DP runs over a single reusable arena
//! (`DpScratch`) so a whole spectrum's work does no per-node allocation. The convolution
//! arithmetic (`axpy`) is unchanged from the direct MS-GF+ port and reproduces its numbers
//! bit-for-bit; only the data layout is optimized.

pub mod graph;

/// The shift-add convolution kernel shared by every distribution update: for each source score
/// `t`, `dst[t + score_diff] += src[t] * aa_prob`. This is the single source of truth for the
/// DP's floating-point arithmetic — keep it byte-identical to preserve MS-GF+ reproduction. Each
/// call writes a distinct `dst` index per `t`, so it vectorizes without changing results.
/// The convolution kernel `d[i] += s[i] * c` over equal-length slices. Selected at runtime so the
/// default `--release` build gets vectorized codegen without a `target-cpu` flag. **Every kernel
/// must produce bit-identical results** (packed mul then packed add — no FMA contraction — matches
/// the scalar rounding per element), so which one runs never changes the SpecEValue.
type AxpyKernel = unsafe fn(&mut [f64], &[f64], f64);

/// Portable fallback: a plain accumulate loop (LLVM may still auto-vectorize to SSE2).
unsafe fn axpy_scalar(d: &mut [f64], s: &[f64], c: f64) {
    for (di, si) in d.iter_mut().zip(s) {
        *di += *si * c;
    }
}

/// AVX kernel: 4 lanes of f64 per iteration (`vmulpd` + `vaddpd`, no `vfmadd`). Caller must ensure
/// the CPU has AVX (guaranteed by [`select_axpy`]).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx")]
unsafe fn axpy_avx(d: &mut [f64], s: &[f64], c: f64) {
    use std::arch::x86_64::*;
    let n = d.len().min(s.len());
    let cc = _mm256_set1_pd(c);
    let dp = d.as_mut_ptr();
    let sp = s.as_ptr();
    let mut k = 0usize;
    while k + 4 <= n {
        let sv = _mm256_loadu_pd(sp.add(k));
        let dv = _mm256_loadu_pd(dp.add(k));
        _mm256_storeu_pd(dp.add(k), _mm256_add_pd(dv, _mm256_mul_pd(sv, cc)));
        k += 4;
    }
    while k < n {
        *dp.add(k) += *sp.add(k) * c;
        k += 1;
    }
}

/// Pick the fastest available convolution kernel once (cheap; `is_x86_feature_detected!` caches).
#[inline]
fn select_axpy() -> AxpyKernel {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx") {
            return axpy_avx;
        }
    }
    axpy_scalar
}

/// Windowing shared by the DP and the sink/cleavage merges: `dst[t + score_diff] += src[t] * aa_prob`
/// over the overlapping score range. `kernel` is the (bit-identical) accumulate primitive.
#[inline(always)]
fn axpy_with(
    kernel: AxpyKernel,
    dst: &mut [f64],
    dst_min: i32,
    src: &[f64],
    src_min: i32,
    score_diff: i32,
    aa_prob: f64,
) {
    let lo = src_min.max(dst_min - score_diff);
    let hi = src_min + src.len() as i32;
    if lo >= hi {
        return;
    }
    let src_lo = (lo - src_min) as usize;
    let dst_lo = (lo + score_diff - dst_min) as usize;
    let len = (hi - lo) as usize;
    let s = &src[src_lo..src_lo + len];
    let d = &mut dst[dst_lo..dst_lo + len];
    // Safety: `d.len() == s.len() == len`; the kernel only touches those in-bounds lanes.
    unsafe { kernel(d, s, aa_prob) }
}

/// Fast path for the DP, whose destination range is constructed as the union of every shifted
/// predecessor range. Unlike [`axpy_with`], no overlap clipping is needed: `src` is known to fit
/// completely at `dst_lo`.
#[inline(always)]
fn axpy_full_with(
    kernel: AxpyKernel,
    dst: &mut [f64],
    dst_min: i32,
    src: &[f64],
    src_min: i32,
    score_diff: i32,
    aa_prob: f64,
) {
    let dst_lo = (src_min + score_diff - dst_min) as usize;
    debug_assert!(dst_lo + src.len() <= dst.len());
    // Safety: the range construction guarantees that the whole source fits at `dst_lo`.
    let d = unsafe { std::slice::from_raw_parts_mut(dst.as_mut_ptr().add(dst_lo), src.len()) };
    // Safety: `d` and `src` have the same length.
    unsafe { kernel(d, src, aa_prob) }
}

/// Scalar windowing convolution for cold paths (sink merge, cleavage, `ScoreDist::add_prob_dist`).
#[inline(always)]
fn axpy(dst: &mut [f64], dst_min: i32, src: &[f64], src_min: i32, score_diff: i32, aa_prob: f64) {
    axpy_with(axpy_scalar, dst, dst_min, src, src_min, score_diff, aa_prob);
}

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

    /// Add a raw distribution slice (spanning `[src_min, src_min + src.len())`), shifted by
    /// `score_diff` and scaled by `aa_prob`, into this distribution.
    #[inline]
    fn add_slice(&mut self, src_min: i32, src: &[f64], score_diff: i32, aa_prob: f64) {
        axpy(
            &mut self.probs,
            self.min_score,
            src,
            src_min,
            score_diff,
            aa_prob,
        );
    }

    /// Add `other`, shifted by `score_diff` and scaled by `aa_prob`, into this distribution.
    /// Mirrors `ScoreDist.addProbDist`.
    pub fn add_prob_dist(&mut self, other: &ScoreDist, score_diff: i32, aa_prob: f64) {
        self.add_slice(other.min_score, &other.probs, score_diff, aa_prob);
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

/// The de novo graph in flat CSR form: node scores plus a single edge list addressed by per-node
/// offsets. Node 0 is the source; nodes are in topological (increasing nominal mass) order. The
/// incoming edges of node `i` are the parallel slices `edge_prev/edge_score/edge_prob` over
/// `edge_start[i]..edge_start[i + 1]`, in amino-acid insertion order — the order the DP sums them,
/// preserved for bit-exact reproduction of MS-GF+. Replaces the old `Vec<Node>`/per-node
/// `Vec<Edge>` (which reallocated on every `push`) with five buffers allocated once per graph.
#[derive(Debug, Clone, Default)]
pub struct Graph {
    pub node_score: Vec<i32>,
    pub edge_start: Vec<u32>,
    pub edge_prev: Vec<u32>,
    pub edge_score: Vec<i32>,
    pub edge_prob: Vec<f64>,
}

/// One node for [`Graph::from_adj`]: its node score and incoming `(prev, edge_score, aa_prob)` edges.
pub type AdjNode = (i32, Vec<(usize, i32, f64)>);

impl Graph {
    /// Number of nodes (including the source).
    #[inline]
    pub fn n_nodes(&self) -> usize {
        self.node_score.len()
    }

    /// Number of edges.
    #[inline]
    pub fn n_edges(&self) -> usize {
        self.edge_prev.len()
    }

    /// Convenience builder from an adjacency list — `nodes[i] = (node_score, incoming edges)` where
    /// each edge is `(prev, edge_score, aa_prob)`. Utility/test helper; the hot-path builder is
    /// [`graph::build_reverse_graph`].
    pub fn from_adj(nodes: &[AdjNode]) -> Self {
        let n = nodes.len();
        let mut g = Graph {
            node_score: Vec::with_capacity(n),
            edge_start: Vec::with_capacity(n + 1),
            edge_prev: Vec::new(),
            edge_score: Vec::new(),
            edge_prob: Vec::new(),
        };
        g.edge_start.push(0);
        for (ns, edges) in nodes {
            g.node_score.push(*ns);
            for &(prev, es, prob) in edges {
                g.edge_prev.push(prev as u32);
                g.edge_score.push(es);
                g.edge_prob.push(prob);
            }
            g.edge_start.push(g.edge_prev.len() as u32);
        }
        g
    }
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

/// One intermediate node distribution, addressed as a slice `[start, start + len)` of the
/// [`DpScratch`] arena, spanning scores `[min_score, min_score + len)`. `len == 0` means the node
/// was never reached.
#[derive(Clone, Copy)]
struct NodeDist {
    min_score: i32,
    start: u32,
    len: u32,
}

impl NodeDist {
    const ABSENT: NodeDist = NodeDist {
        min_score: 0,
        start: 0,
        len: 0,
    };
}

/// Reusable scratch for the DP: an arena backing every intermediate node distribution, so a whole
/// spectrum's generating function does **no per-node allocation**. Create once and reuse across
/// candidate masses and spectra (one instance per thread); [`compute`] allocates a throwaway one.
#[derive(Default)]
pub struct DpScratch {
    arena: Vec<f64>,
    dists: Vec<NodeDist>,
}

impl DpScratch {
    /// Total score-distribution cells written for the last `compute_into` — i.e. the sum of all
    /// reachable nodes' support widths (the DP's convolution work). Profiling aid.
    pub fn arena_len(&self) -> usize {
        self.arena.len()
    }

    /// Number of reachable nodes in the last `compute_into`. Profiling aid.
    pub fn reachable(&self) -> usize {
        self.dists.iter().filter(|d| d.len > 0).count()
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

/// Compute the generating function over `graph` (topological order, node 0 = source), summing the
/// distributions of `sinks` and applying the neighboring-AA `cleavage` weighting. Allocates a
/// throwaway [`DpScratch`]; use [`compute_into`] with a reused scratch on hot paths. Mirrors
/// `GeneratingFunction.computeGeneratingFunction`. Returns `None` if the sinks are unreachable.
pub fn compute(graph: &Graph, sinks: &[usize], cleavage: Option<Cleavage>) -> Option<GenFunc> {
    let mut sc = DpScratch::default();
    compute_into(&mut sc, graph, sinks, cleavage)
}

/// Like [`compute`] but reuses `sc` for all intermediate distributions — no per-node allocation.
pub fn compute_into(
    sc: &mut DpScratch,
    graph: &Graph,
    sinks: &[usize],
    cleavage: Option<Cleavage>,
) -> Option<GenFunc> {
    let n_full = graph.n_nodes();
    if n_full == 0 {
        return None;
    }
    // Only visit nodes up to this candidate's largest sink. A graph built for the largest
    // isotope-error candidate serves the smaller ones by processing a prefix of it.
    let n = sinks
        .iter()
        .copied()
        .max()
        .map(|s| s + 1)
        .unwrap_or(n_full)
        .min(n_full);
    let kernel = select_axpy(); // chosen once; the per-edge convolution below dominates the DP
    sc.arena.clear();
    sc.dists.clear();
    sc.dists.resize(n, NodeDist::ABSENT);

    // Source: point mass (prob 1) at score 0.
    sc.arena.push(1.0);
    sc.dists[0] = NodeDist {
        min_score: 0,
        start: 0,
        len: 1,
    };

    for i in 1..n {
        let node_score = graph.node_score[i];
        let (e0, e1) = (
            graph.edge_start[i] as usize,
            graph.edge_start[i + 1] as usize,
        );
        // The sink's incoming edges carry errorScore 0 (setBackwardEdgesFromSink). The sink differs
        // per candidate, so it is zeroed here rather than baked into the (shared) edge array.
        let is_sink = sinks.contains(&i);

        // Score range of this node's distribution, from its reachable predecessors.
        let (mut cur_min, mut cur_max) = (i32::MAX, i32::MIN);
        for e in e0..e1 {
            let pd = sc.dists[graph.edge_prev[e] as usize];
            if pd.len == 0 {
                continue;
            }
            let combined = node_score + if is_sink { 0 } else { graph.edge_score[e] };
            cur_min = cur_min.min(pd.min_score + combined);
            cur_max = cur_max.max(pd.min_score + pd.len as i32 + combined);
        }
        if cur_min >= cur_max {
            continue; // unreachable
        }

        let len = (cur_max - cur_min) as usize;
        let start = sc.arena.len();
        sc.arena.resize(start + len, 0.0);

        // Predecessors live earlier in the arena (append order), so split at the new node's start:
        // `prev_part` is every earlier distribution (immutable), `cur` is this node (mutable).
        {
            let (prev_part, cur_part) = sc.arena.split_at_mut(start);
            let cur = &mut cur_part[..len];
            for e in e0..e1 {
                let pd = sc.dists[graph.edge_prev[e] as usize];
                if pd.len == 0 {
                    continue;
                }
                let src_start = pd.start as usize;
                let src_len = pd.len as usize;
                debug_assert!(src_start + src_len <= prev_part.len());
                // Safety: every `NodeDist` is recorded immediately after its arena range is
                // appended, and predecessor ranges are wholly before `start`.
                let src = unsafe {
                    std::slice::from_raw_parts(prev_part.as_ptr().add(src_start), src_len)
                };
                let score_diff = node_score + if is_sink { 0 } else { graph.edge_score[e] };
                axpy_full_with(
                    kernel,
                    cur,
                    cur_min,
                    src,
                    pd.min_score,
                    score_diff,
                    graph.edge_prob[e],
                );
            }
        }
        sc.dists[i] = NodeDist {
            min_score: cur_min,
            start: start as u32,
            len: len as u32,
        };
    }

    // Merge the sink distributions.
    let (mut min, mut max) = (i32::MAX, i32::MIN);
    for &s in sinks {
        let d = sc.dists[s];
        if d.len == 0 {
            continue;
        }
        min = min.min(d.min_score);
        max = max.max(d.min_score + d.len as i32);
    }
    if max <= min {
        return None;
    }
    let mut merged = ScoreDist::new(min, max);
    for &s in sinks {
        let d = sc.dists[s];
        if d.len == 0 {
            continue;
        }
        let src = &sc.arena[d.start as usize..d.start as usize + d.len as usize];
        merged.add_slice(d.min_score, src, 0, 1.0);
    }

    // Neighboring amino-acid cleavage credit/penalty (probabilistic).
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
        // source(0) --aa(0.5),es1--> node1(ns2) --aa(0.5),es3--> node2(ns0) --aa(1.0),es0--> sink3(ns0)
        // The sink's incoming edge is scored 0 (setBackwardEdgesFromSink), which compute enforces, so
        // the two scored interior edges determine the path.
        let g = Graph::from_adj(&[
            (0, vec![]),            // source
            (2, vec![(0, 1, 0.5)]), // node1
            (0, vec![(1, 3, 0.5)]), // node2
            (0, vec![(2, 0, 1.0)]), // sink3
        ]);
        let gf = compute(&g, &[3], None).unwrap();
        // path score = (2+1) + (0+3) + (0+0) = 6 with prob 0.5*0.5*1.0 = 0.25
        approx(gf.spectral_probability(6), 0.25);
        approx(gf.spectral_probability(7), 0.0);
        assert_eq!(gf.max_score(), 6); // max achievable score
    }

    #[test]
    fn compute_zeros_sink_edge_scores() {
        // A scored edge into the sink is ignored (its score is forced to 0), matching MS-GF+.
        let g = Graph::from_adj(&[
            (0, vec![]),            // source
            (0, vec![(0, 5, 1.0)]), // sink1 — edge score 5 must NOT count
        ]);
        let gf = compute(&g, &[1], None).unwrap();
        approx(gf.spectral_probability(0), 1.0); // score 0, not 5
        assert_eq!(gf.max_score(), 0);
    }

    #[test]
    fn neighboring_cleavage_splits_mass() {
        // merged: score 0 prob 1. cleavage credit +2 (p=0.25), penalty -1 (p=0.75)
        let g = Graph::from_adj(&[(0, vec![]), (0, vec![(0, 0, 1.0)])]);
        let gf = compute(
            &g,
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
