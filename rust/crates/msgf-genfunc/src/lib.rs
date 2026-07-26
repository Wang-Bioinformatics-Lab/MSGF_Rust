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
pub mod tilt;

/// The shift-add convolution kernel shared by every distribution update: for each source score
/// `t`, `dst[t + score_diff] += src[t] * aa_prob`. This is the single source of truth for the
/// DP's floating-point arithmetic — keep it byte-identical to preserve MS-GF+ reproduction. Each
/// call writes a distinct `dst` index per `t`, so it vectorizes without changing results.
/// The convolution kernel `d[i] += s[i] * c` over equal-length slices. Selected per DP call so the
/// default `--release` build gets vectorized codegen without a `target-cpu` flag. **Every kernel
/// must produce bit-identical results** (packed mul then packed add — no FMA contraction — matches
/// the scalar rounding per element), so which one runs never changes the SpecEValue.
///
/// Portable fallback: a plain accumulate loop (LLVM may still auto-vectorize to SSE2).
#[inline(always)]
unsafe fn axpy_scalar(d: &mut [f64], s: &[f64], c: f64) {
    for (di, si) in d.iter_mut().zip(s) {
        *di += *si * c;
    }
}

/// AVX kernel: 4 lanes of f64 per iteration (`vmulpd` + `vaddpd`, no `vfmadd`). Caller must ensure
/// the CPU has AVX (guaranteed by [`compute_inner`]'s once-per-call dispatch).
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

/// The DP's convolution kernel, selected at *compile* time rather than through a function
/// pointer: [`compute_inner`] monomorphizes the whole DP body on `AVX` and dispatches once per call,
/// so the kernel can inline into the edge loop.
///
/// Both instantiations must stay bit-identical (packed multiply then packed add, never an FMA);
/// `avx_matches_scalar_bitwise` pins that.
#[inline(always)]
unsafe fn axpy_lane<const AVX: bool>(d: &mut [f64], s: &[f64], c: f64) {
    #[cfg(target_arch = "x86_64")]
    if AVX {
        return axpy_avx(d, s, c);
    }
    axpy_scalar(d, s, c)
}

/// The DP's per-edge convolution. `CLIP` selects between two cases the caller knows statically:
///
/// - `false` — the destination range is the union of every shifted predecessor range, so `src` is
///   known to fit whole at `dst_lo` and no overlap arithmetic is needed;
/// - `true` — tail pruning raised the destination floor above that union, so the source has to be
///   clipped to the surviving window.
#[inline(always)]
fn axpy_edge<const AVX: bool, const CLIP: bool>(
    dst: &mut [f64],
    dst_min: i32,
    src: &[f64],
    src_min: i32,
    score_diff: i32,
    aa_prob: f64,
) {
    let (src_lo, dst_lo, len) = if CLIP {
        let lo = src_min.max(dst_min - score_diff);
        let hi = src_min + src.len() as i32;
        if lo >= hi {
            return;
        }
        (
            (lo - src_min) as usize,
            (lo + score_diff - dst_min) as usize,
            (hi - lo) as usize,
        )
    } else {
        (0, (src_min + score_diff - dst_min) as usize, src.len())
    };
    debug_assert!(src_lo + len <= src.len() && dst_lo + len <= dst.len());
    // Safety: both ranges are in bounds by the case analysis above, and `d.len() == s.len()`.
    unsafe {
        let d = std::slice::from_raw_parts_mut(dst.as_mut_ptr().add(dst_lo), len);
        let s = std::slice::from_raw_parts(src.as_ptr().add(src_lo), len);
        axpy_lane::<AVX>(d, s, aa_prob)
    }
}

/// Scalar windowing convolution for cold paths (sink merge, cleavage, `ScoreDist::add_prob_dist`).
#[inline(always)]
fn axpy(dst: &mut [f64], dst_min: i32, src: &[f64], src_min: i32, score_diff: i32, aa_prob: f64) {
    axpy_edge::<false, true>(dst, dst_min, src, src_min, score_diff, aa_prob);
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
/// incoming edges of node `i` are the parallel slices `edge_prev/edge_score/edge_aa` over
/// `edge_start[i]..edge_start[i + 1]`, in amino-acid insertion order — the order the DP sums them,
/// preserved for bit-exact reproduction of MS-GF+. Replaces the old `Vec<Node>`/per-node
/// `Vec<Edge>` (which reallocated on every `push`) with flat buffers allocated once per graph.
///
/// The per-edge amino-acid **probability** is *not* stored per edge: it takes one of `|alphabet|`
/// (~21) values, so the edge carries a two-byte index `edge_aa` into the `aa_prob` table and the DP
/// reads `aa_prob[edge_aa[e] as usize]`. That is the identical `f64` bit pattern the old
/// `edge_prob: Vec<f64>` held — the arithmetic is unchanged — at 1/4 the memory traffic. The index
/// is explicit rather than derived from the nominal mass because distinct amino acids share nominal
/// masses (I/L at 113, K/Q at 128), so two parallel edges can differ only by which one they are.
#[derive(Debug, Clone, Default)]
pub struct Graph {
    pub node_score: Vec<i32>,
    pub edge_start: Vec<u32>,
    pub edge_prev: Vec<u32>,
    pub edge_score: Vec<i32>,
    /// Per-edge index into [`Graph::aa_prob`], in amino-acid insertion order.
    pub edge_aa: Vec<u16>,
    /// Amino-acid background probabilities, indexed by [`Graph::edge_aa`].
    pub aa_prob: Vec<f64>,
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
            edge_aa: Vec::new(),
            aa_prob: Vec::new(),
        };
        g.edge_start.push(0);
        for (ns, edges) in nodes {
            g.node_score.push(*ns);
            for &(prev, es, prob) in edges {
                // Intern the probability by exact bit pattern — the DP must read back the same
                // `f64`, so equality is `to_bits`, never `==` (which conflates ±0.0 and traps NaN).
                let bits = prob.to_bits();
                let idx = match g.aa_prob.iter().position(|p| p.to_bits() == bits) {
                    Some(i) => i,
                    None => {
                        g.aa_prob.push(prob);
                        g.aa_prob.len() - 1
                    }
                };
                assert!(
                    idx <= u16::MAX as usize,
                    "Graph supports at most {} distinct edge probabilities",
                    u16::MAX
                );
                g.edge_prev.push(prev as u32);
                g.edge_score.push(es);
                g.edge_aa.push(idx as u16);
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
    /// Lowest score this distribution answers for. [`compute`] / [`compute_into`] produce the whole
    /// distribution and set `i32::MIN`; the pruning entry points set the threshold they were given.
    /// Querying below it is meaningless — the cells are absent by construction.
    pub valid_from: i32,
    /// Upper bound on how much probability mass the pruning discarded, i.e. the true tail at
    /// [`Self::valid_from`] lies in `[p, p + err_bound]`. Exactly `0.0` for every exact path
    /// ([`compute`], [`compute_into`], [`compute_tail_into`]); only
    /// [`compute_tail_with`] can make it nonzero.
    pub err_bound: f64,
    /// DeNovoScore, taken from the integer `max_remaining` sweep rather than from `dist` — a top
    /// [`Prune::cap`] removes the maximum-score cells, so the distribution can no longer supply it.
    /// `i32::MIN` on the unpruned paths, where `dist` is complete and does supply it.
    denovo: i32,
}

impl GenFunc {
    /// Spectral probability (the p-value) at the observed RawScore.
    /// Querying below [`Self::valid_from`] is a programming error, not a smaller number: the cells
    /// under the pruning threshold were never computed, so the sum would be silently too small.
    /// This is asserted in **release** too — on the pruned path a wrong SpecEValue is far worse
    /// than a panic, and the check is one predictable branch per PSM against a whole DP.
    pub fn spectral_probability(&self, raw_score: i32) -> f64 {
        assert!(
            self.valid_from == i32::MIN || raw_score >= self.valid_from,
            "spectral_probability({raw_score}) below the pruning threshold {}: those cells were \
             never computed",
            self.valid_from
        );
        self.dist.spectral_probability(raw_score)
    }
    /// Maximum achievable score (DeNovoScore). Exact on every path — under pruning it comes from
    /// the integer `max_remaining` sweep, which is unaffected by anything the DP discards.
    pub fn max_score(&self) -> i32 {
        if self.denovo != i32::MIN {
            self.denovo
        } else {
            self.dist.max_score() - 1
        }
    }
    /// Relative width of the certified interval around [`Self::spectral_probability`] at
    /// [`Self::valid_from`]: `0.0` when exact, else `err_bound / p`.
    pub fn relative_error(&self) -> f64 {
        if self.err_bound == 0.0 {
            return 0.0;
        }
        let p = self.dist.spectral_probability(self.valid_from);
        if p > 0.0 {
            self.err_bound / p
        } else {
            f64::INFINITY
        }
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

/// One incoming edge of the node currently being relaxed, fully resolved: the arena range of the
/// predecessor distribution, its score origin, the shift this edge applies, and the amino-acid
/// probability. See [`MAX_FUSED_IN`] for why these are cached.
#[derive(Clone, Copy)]
struct EdgeDesc {
    /// Arena offset of the predecessor's distribution (`NodeDist::start`).
    src_start: u32,
    /// Length of the predecessor's distribution (`NodeDist::len`, always non-zero here).
    src_len: u32,
    /// Score of the predecessor's first cell (`NodeDist::min_score`).
    src_min: i32,
    /// `node_score + edge_score` — what the convolution shifts by.
    score_diff: i32,
    /// The edge's amino-acid probability, resolved through `Graph::aa_prob[edge_aa[e]]`.
    prob: f64,
}

impl EdgeDesc {
    const ZERO: EdgeDesc = EdgeDesc {
        src_start: 0,
        src_len: 0,
        src_min: 0,
        score_diff: 0,
        prob: 0.0,
    };
}

/// How many incoming edges the DP will cache on the stack between its two passes over a node.
///
/// The DP visits each incoming edge twice — once to compute the destination's score range, once to
/// convolve — and both passes gather the same `NodeDist` (a scattered load, indexed by
/// `edge_prev`) and re-derive the same `node_score + edge_score`. Caching the resolved descriptors
/// during the range pass lets the convolution run straight off the stack, removing the second
/// gather and the re-derivation. This changes *which loads happen*, never the arithmetic: the same
/// edges are convolved in the same order with the same `score_diff` and `prob`.
///
/// A de novo graph node has one incoming edge per amino acid (20 standard + modified forms), so 32
/// covers every realistic model. A node with more falls back to the original two-pass code rather
/// than growing the array or allocating.
const MAX_FUSED_IN: usize = 32;

/// Reusable scratch for the DP: an arena backing every intermediate node distribution, so a whole
/// spectrum's generating function does **no per-node allocation**. Create once and reuse across
/// candidate masses and spectra (one instance per thread); [`compute`] allocates a throwaway one.
#[derive(Default)]
pub struct DpScratch {
    arena: Vec<f64>,
    dists: Vec<NodeDist>,
    /// `max_remaining` — best score still obtainable from each node to a sink.
    rem: Vec<i32>,
    /// `max_achievable` — best score reachable at each node from the source. Only filled when
    /// [`Prune::skip_dead_nodes`] is on.
    ach: Vec<i32>,
    /// Backward tilted sums, for bounded pruning.
    tilt: tilt::TiltScratch,
    /// `e^{θ·d}` for `d = score − cut`, indexed by `d − etab_lo`.
    etab: Vec<f64>,
    etab_lo: i32,
    /// Cells the last bounded run *retained*; the DP's real convolution work.
    cells: usize,
}

impl DpScratch {
    /// Total score-distribution cells written for the last `compute_into` — i.e. the sum of all
    /// reachable nodes' support widths. Profiling aid.
    pub fn arena_len(&self) -> usize {
        self.arena.len()
    }

    /// Cells *retained* — the ones successors actually read, so the quantity the DP's convolution
    /// work is proportional to. Equals [`Self::arena_len`] unless bounded pruning trimmed rows.
    pub fn cells(&self) -> usize {
        self.cells
    }

    /// Number of reachable nodes in the last `compute_into`. Profiling aid.
    pub fn reachable(&self) -> usize {
        self.dists.iter().filter(|d| d.len > 0).count()
    }

    /// Tilted sweeps consumed since the last [`Self::reset_sweeps`]. Profiling aid.
    pub fn sweeps(&self) -> u32 {
        self.tilt.sweeps
    }

    pub fn reset_sweeps(&mut self) {
        self.tilt.reset_sweeps();
    }

    /// The backward-tilted-sweep buffers, for callers driving [`tilt::solve_theta`] themselves
    /// before a [`Prune`] with a `tilt` policy. The sweep must be run on the same graph, sink set
    /// and `θ` the `Prune` names, or the DP ignores it.
    pub fn tilt_mut(&mut self) -> &mut tilt::TiltScratch {
        &mut self.tilt
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
    Some(GenFunc {
        dist: merged,
        // The merged tail is only as trustworthy as its least-complete member, and the discarded
        // masses add.
        valid_from: gfs.iter().map(|g| g.valid_from).max().unwrap_or(i32::MIN),
        err_bound: gfs.iter().map(|g| g.err_bound).sum(),
        denovo: gfs.iter().map(|g| g.denovo).max().unwrap_or(i32::MIN),
    })
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
    compute_inner::<false>(sc, graph, sinks, cleavage, None)
}

/// Upper-tail-only generating function: identical to [`compute_into`] for every score at or above
/// `threshold`, but it never materializes the score cells that provably cannot reach it.
///
/// **Why this is exact.** Let `max_rem[m]` be the best score any path can still earn between node
/// `m` and a sink. A cell `(m, s)` can only ever end at a final score `<= s + max_rem[m]`, so if
/// `s + max_rem[m] < threshold` every path through it lands strictly below `threshold` and it
/// cannot contribute to the tail. Dropping it is not an approximation: the retained cells are
/// reached by exactly the same multiply-adds, in exactly the same order, as in the unpruned DP, so
/// they are **bit-identical** — `spectral_probability(t)` for `t >= threshold` returns the same
/// `f64` the full computation would. Only the discarded low-score cells differ (they are absent).
///
/// The `threshold` is what a search already knows before it needs a p-value: the RawScore of the
/// worst PSM it will report for this spectrum. Cost falls as the threshold approaches the
/// DeNovoScore — i.e. the DP gets cheaper exactly as the match gets better. The bound is computed
/// by one integer sweep of the CSR (`max_remaining`), ~1% of the DP's own cost.
///
/// `cleavage` shifts the merged distribution by `credit`/`penalty`, so the threshold is lowered by
/// `credit` internally to keep the credited branch intact.
pub fn compute_tail_into(
    sc: &mut DpScratch,
    graph: &Graph,
    sinks: &[usize],
    cleavage: Option<Cleavage>,
    threshold: i32,
) -> Option<GenFunc> {
    compute_tail_with(sc, graph, sinks, cleavage, Prune::exact(threshold))
}

/// Run the tail DP under an explicit [`Prune`] policy — the general form of [`compute_tail_into`].
///
/// With [`Prune::exact`] the result is bit-identical to the full DP above the threshold. Add a
/// [`Prune::cap`] or a [`Prune::tilt`] and it becomes **one-sided and certified**: probability is
/// only ever removed, so the returned `p` is a lower bound and the true tail at the threshold lies
/// in `[p, p + err_bound]`, with [`GenFunc::err_bound`] the summed bound of everything discarded —
/// an accumulated quantity, not an estimate.
///
/// `DeNovoScore` ([`GenFunc::max_score`]) stays exact under every policy: it is read off the
/// integer `max_remaining` sweep, not off the (possibly capped) distribution.
pub fn compute_tail_with(
    sc: &mut DpScratch,
    graph: &Graph,
    sinks: &[usize],
    cleavage: Option<Cleavage>,
    prune: Prune,
) -> Option<GenFunc> {
    compute_inner::<true>(sc, graph, sinks, cleavage, Some(prune))
}

/// How the tail DP prunes: an exact floor always, plus two optional lossy trims that trade a
/// certified amount of probability for cells.
///
/// After the exact floor, a node's retained score window is `[cut − max_rem[m], cut + cap]`. The
/// two ends behave completely differently, and measurement (see `examples/prunelab.rs`) says so:
///
/// - **The top (`cap`) is where the slack is.** Cells far above the threshold are the rare,
///   high-scoring ones: `P(m,s)` decays by roughly `e^{−θ} ≈ 0.6` per score unit there, so the
///   whole discarded mass is bounded by the geometric remainder, and 30-odd units above the
///   threshold already puts it below `1e-6` of the answer. Nothing else in the DP prunes this end,
///   because the exact bound has nothing to say about it — those cells *can* reach the threshold,
///   they simply almost never do. Costs nothing to exploit: `Q_m(r) <= 1` for `r <= 0`, so the
///   discarded contribution is bounded by the discarded probability itself.
/// - **The bottom (`tilt`) is nearly exhausted.** The exact floor already sits close to where the
///   tilted (Chernoff) bound would put it, so budgeting `P(m,s) · e^{−θ(cut−s)} · B_θ(m)` buys only
///   a few percent more cells — and the backward tilted sweep it needs ([`tilt`]) costs ~13% of the
///   DP. Left in as an opt-in knob because it is the natural thing to reach for; it is not
///   recommended, and `Prune::capped` does not enable it.
#[derive(Clone, Copy)]
pub struct Prune {
    /// Lowest score whose tail probability the caller will read.
    pub threshold: i32,
    /// Score units to retain **above** the threshold, or `None` to retain all of them (exact).
    pub cap: Option<i32>,
    /// `(θ, budget)` for low-end trimming: discard from each node's low end while the summed
    /// `P · e^{−θ(cut−s)} · B_θ(m)` bound stays under `budget`. Requires [`tilt::solve_theta`] to
    /// have filled the scratch's backward sweep at the same `θ`.
    pub tilt: Option<(f64, f64)>,
    /// Add a forward `max_achievable` sweep so nodes that cannot lie on any path clearing the
    /// threshold are skipped without their edges being visited. Exact — it removes only nodes the
    /// score-range pass would have found empty anyway. Whether it wins depends on how many nodes
    /// die: it trades one integer sweep of every edge for skipping the (more expensive) range pass
    /// on dead ones.
    pub skip_dead_nodes: bool,
}

impl Prune {
    /// Threshold-only pruning: bit-identical to the full DP at and above `threshold`.
    pub fn exact(threshold: i32) -> Self {
        Prune {
            threshold,
            cap: None,
            tilt: None,
            skip_dead_nodes: false,
        }
    }

    /// Threshold pruning plus a top cap of `cap` score units. **Lossy and not recommended** — see
    /// [`Prune::cap`].
    pub fn capped(threshold: i32, cap: i32) -> Self {
        Prune {
            cap: Some(cap),
            ..Prune::exact(threshold)
        }
    }
}

/// Nodes the DP visits: a graph built for the largest isotope-error candidate serves the smaller
/// ones by processing a prefix of it.
#[inline]
fn visited_nodes(graph: &Graph, sinks: &[usize]) -> usize {
    let n_full = graph.n_nodes();
    sinks
        .iter()
        .copied()
        .max()
        .map(|s| s + 1)
        .unwrap_or(n_full)
        .min(n_full)
}

/// `ach[m]` = the largest score any source→`m` path can have collected, or `i32::MIN` when `m` is
/// unreachable. Paired with [`max_remaining`] it gives `best[m] = ach[m] + max_rem[m]`, the best
/// full-path score *through* `m` — and a node with `best[m] < cut` cannot appear on any path that
/// clears the threshold, so the DP can skip it **without visiting its edges at all**. That is the
/// only way to get under the DP's fixed per-edge cost; the score-range pass it replaces has to
/// gather a `NodeDist` per edge, where this sweep reads two `i32` streams.
fn max_achievable(ach: &mut Vec<i32>, graph: &Graph, sinks: &[usize], n: usize) {
    ach.clear();
    ach.resize(n, i32::MIN);
    ach[0] = 0;
    for i in 1..n {
        let is_sink = sinks.contains(&i);
        let ns = graph.node_score[i];
        let (e0, e1) = (
            graph.edge_start[i] as usize,
            graph.edge_start[i + 1] as usize,
        );
        let mut best = i32::MIN;
        for e in e0..e1 {
            let p = ach[graph.edge_prev[e] as usize];
            if p == i32::MIN {
                continue;
            }
            let cand = p + ns + if is_sink { 0 } else { graph.edge_score[e] };
            if cand > best {
                best = cand;
            }
        }
        ach[i] = best;
    }
}

/// `max_rem[m]` = the largest score still obtainable on any path from `m` to a sink, or `i32::MIN`
/// when no sink is reachable from `m`. One descending sweep of the same CSR: an edge's `prev` index
/// is always strictly below the node it enters, so every node is finalized before it is read.
/// `max_rem[0]` is the source's best full-path score — the DeNovoScore before cleavage weighting.
fn max_remaining(rem: &mut Vec<i32>, graph: &Graph, sinks: &[usize], n: usize) {
    rem.clear();
    rem.resize(n, i32::MIN);
    for &s in sinks {
        if s < n {
            rem[s] = 0;
        }
    }
    for i in (1..n).rev() {
        let ri = rem[i];
        if ri == i32::MIN {
            continue;
        }
        let is_sink = sinks.contains(&i);
        let base = ri + graph.node_score[i];
        let (e0, e1) = (
            graph.edge_start[i] as usize,
            graph.edge_start[i + 1] as usize,
        );
        for e in e0..e1 {
            let p = graph.edge_prev[e] as usize;
            let cand = base + if is_sink { 0 } else { graph.edge_score[e] };
            if cand > rem[p] {
                rem[p] = cand;
            }
        }
    }
}

/// Dispatch the DP's SIMD capability **once per call** rather than once per edge, and run the
/// matching monomorphization of [`dp`].
fn compute_inner<const PRUNED: bool>(
    sc: &mut DpScratch,
    graph: &Graph,
    sinks: &[usize],
    cleavage: Option<Cleavage>,
    prune: Option<Prune>,
) -> Option<GenFunc> {
    #[cfg(target_arch = "x86_64")]
    if is_x86_feature_detected!("avx") {
        // Safety: guarded by the feature detection above.
        return unsafe { dp_avx::<PRUNED, true>(sc, graph, sinks, cleavage, prune) };
    }
    dp::<false, PRUNED, true>(sc, graph, sinks, cleavage, prune)
}

/// The pre-fusion two-pass DP, kept callable so `fused_edge_pass_is_bit_exact` has something to
/// compare against. Test-only: release builds never instantiate `FUSE = false`.
#[cfg(test)]
fn compute_inner_unfused<const PRUNED: bool>(
    sc: &mut DpScratch,
    graph: &Graph,
    sinks: &[usize],
    cleavage: Option<Cleavage>,
    prune: Option<Prune>,
) -> Option<GenFunc> {
    #[cfg(target_arch = "x86_64")]
    if is_x86_feature_detected!("avx") {
        // Safety: guarded by the feature detection above.
        return unsafe { dp_avx::<PRUNED, false>(sc, graph, sinks, cleavage, prune) };
    }
    dp::<false, PRUNED, false>(sc, graph, sinks, cleavage, prune)
}

/// The AVX instantiation. `#[target_feature]` here is what lets the kernel inline into the edge
/// loop; the body is the same [`dp`] the portable path runs, so the two stay bit-identical.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx")]
unsafe fn dp_avx<const PRUNED: bool, const FUSE: bool>(
    sc: &mut DpScratch,
    graph: &Graph,
    sinks: &[usize],
    cleavage: Option<Cleavage>,
    prune: Option<Prune>,
) -> Option<GenFunc> {
    dp::<true, PRUNED, FUSE>(sc, graph, sinks, cleavage, prune)
}

#[inline(always)]
fn dp<const AVX: bool, const PRUNED: bool, const FUSE: bool>(
    sc: &mut DpScratch,
    graph: &Graph,
    sinks: &[usize],
    cleavage: Option<Cleavage>,
    prune: Option<Prune>,
) -> Option<GenFunc> {
    let n_full = graph.n_nodes();
    if n_full == 0 {
        return None;
    }
    // Only visit nodes up to this candidate's largest sink. A graph built for the largest
    // isotope-error candidate serves the smaller ones by processing a prefix of it.
    let n = visited_nodes(graph, sinks);
    // The ~21-entry probability table stays in L1; the per-edge stream is one byte instead of
    // eight, and the value handed to `axpy` is the identical `f64` the old `edge_prob[e]` held.
    let aa_prob = graph.aa_prob.as_slice();

    // Tail pruning: the per-node score floor `floor[i] = cut - max_rem[i]`. `cut` is clamped to the
    // best achievable score so the DeNovoScore cell always survives and `max_score()` stays exact
    // even when the caller asks for a threshold no peptide can reach.
    let mut cut = 0i32;
    let mut ceiling = i32::MAX;
    let mut theta = 0.0f64;
    let mut budget = 0.0f64;
    let mut denovo = i32::MIN;
    let mut skip_dead = false;
    let credit = credit_of(cleavage);
    if PRUNED {
        let p = prune.expect("PRUNED implies a Prune config");
        max_remaining(&mut sc.rem, graph, sinks, n);
        skip_dead = p.skip_dead_nodes;
        if skip_dead {
            max_achievable(&mut sc.ach, graph, sinks, n);
        }
        denovo = sc.rem[0];
        if denovo == i32::MIN {
            return None; // no sink is reachable from the source
        }
        cut = p.threshold.saturating_sub(credit).min(denovo);
        if let Some(cap) = p.cap {
            ceiling = cut.saturating_add(cap.max(0));
        }
        if let Some((t, b)) = p.tilt {
            if t > 0.0 && b > 0.0 && sc.tilt.b.len() >= n {
                theta = t;
                budget = b;
                build_etab(sc, theta, denovo);
            }
        }
    }
    let mut err_bound = 0.0f64;

    sc.arena.clear();
    sc.cells = 0;
    sc.dists.clear();
    sc.dists.resize(n, NodeDist::ABSENT);

    // Source: point mass (prob 1) at score 0.
    sc.arena.push(1.0);
    sc.cells = 1;
    sc.dists[0] = NodeDist {
        min_score: 0,
        start: 0,
        len: 1,
    };

    // Scratch for the fused edge pass, hoisted out of the node loop so no node pays to initialize
    // it. Written by the score-range pass, read by the convolution; see [`MAX_FUSED_IN`].
    let mut descs = [EdgeDesc::ZERO; MAX_FUSED_IN];

    for i in 1..n {
        // Cheapest possible rejection: this node is on no path that can clear the threshold, so
        // none of its ~21 edges need to be looked at.
        if PRUNED && skip_dead {
            let (a, r) = (sc.ach[i], sc.rem[i]);
            if a == i32::MIN || r == i32::MIN || a + r < cut {
                continue;
            }
        }
        let node_score = graph.node_score[i];
        let (e0, e1) = (
            graph.edge_start[i] as usize,
            graph.edge_start[i + 1] as usize,
        );
        // The sink's incoming edges carry errorScore 0 (setBackwardEdgesFromSink). The sink differs
        // per candidate, so it is zeroed here rather than baked into the (shared) edge array.
        let is_sink = sinks.contains(&i);

        // Score range of this node's distribution, from its reachable predecessors. When the node's
        // in-degree fits, this pass also records each live edge's resolved descriptor so the
        // convolution below never has to gather `sc.dists[edge_prev[e]]` a second time. The
        // descriptors are appended in `e0..e1` order and skip exactly the edges the convolution
        // skips (`pd.len == 0`), so the convolution's edge order — and therefore its summation
        // order — is unchanged.
        let fused = FUSE && e1 - e0 <= MAX_FUSED_IN;
        let mut n_desc = 0usize;
        let (mut cur_min, mut cur_max) = (i32::MAX, i32::MIN);
        if fused {
            for e in e0..e1 {
                let pd = sc.dists[graph.edge_prev[e] as usize];
                if pd.len == 0 {
                    continue;
                }
                let combined = node_score + if is_sink { 0 } else { graph.edge_score[e] };
                cur_min = cur_min.min(pd.min_score + combined);
                cur_max = cur_max.max(pd.min_score + pd.len as i32 + combined);
                // `n_desc <= e - e0 < MAX_FUSED_IN` by the `fused` guard.
                descs[n_desc] = EdgeDesc {
                    src_start: pd.start,
                    src_len: pd.len,
                    src_min: pd.min_score,
                    score_diff: combined,
                    prob: aa_prob[graph.edge_aa[e] as usize],
                };
                n_desc += 1;
            }
        } else {
            for e in e0..e1 {
                let pd = sc.dists[graph.edge_prev[e] as usize];
                if pd.len == 0 {
                    continue;
                }
                let combined = node_score + if is_sink { 0 } else { graph.edge_score[e] };
                cur_min = cur_min.min(pd.min_score + combined);
                cur_max = cur_max.max(pd.min_score + pd.len as i32 + combined);
            }
        }
        // Raise the floor to the lowest score that can still reach `cut` (see `compute_tail_into`).
        // `max_rem[i] == i32::MIN` means no sink lies beyond this node, so it is dropped entirely.
        // `clipped` records whether the floor actually moved: when it did not, the destination is
        // still the plain union of the shifted predecessors and the convolution can skip its
        // overlap arithmetic, exactly as in the unpruned DP.
        let mut clipped = false;
        if PRUNED {
            let rem = sc.rem[i];
            if rem == i32::MIN {
                continue;
            }
            let floor = cut - rem;
            if floor > cur_min {
                cur_min = floor;
                clipped = true;
            }
        }
        if cur_min >= cur_max {
            continue; // unreachable, or wholly below the tail threshold
        }

        let len = (cur_max - cur_min) as usize;
        let start = sc.arena.len();
        sc.arena.resize(start + len, 0.0);

        // Predecessors live earlier in the arena (append order), so split at the new node's start:
        // `prev_part` is every earlier distribution (immutable), `cur` is this node (mutable).
        {
            let (prev_part, cur_part) = sc.arena.split_at_mut(start);
            let cur = &mut cur_part[..len];
            // Two monomorphizations of one loop: `CLIP` is a compile-time constant inside it, so
            // the overlap arithmetic disappears entirely from the common (unclipped) case.
            //
            // `relax_fused` consumes the descriptors the range pass already resolved; `relax`
            // is the fallback for a node whose in-degree exceeded `MAX_FUSED_IN`. Both issue the
            // same `axpy_edge` calls, in the same order, with the same arguments.
            macro_rules! relax_fused {
                ($clip:expr) => {
                    for d in &descs[..n_desc] {
                        let src_start = d.src_start as usize;
                        let src_len = d.src_len as usize;
                        debug_assert!(src_start + src_len <= prev_part.len());
                        // Safety: every `NodeDist` is recorded immediately after its arena range is
                        // appended, and predecessor ranges are wholly before `start`. The
                        // descriptor was copied from such a `NodeDist` earlier in this iteration,
                        // and the arena has not been resized since the split.
                        let src = unsafe {
                            std::slice::from_raw_parts(prev_part.as_ptr().add(src_start), src_len)
                        };
                        axpy_edge::<AVX, $clip>(cur, cur_min, src, d.src_min, d.score_diff, d.prob);
                    }
                };
            }
            macro_rules! relax {
                ($clip:expr) => {
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
                        axpy_edge::<AVX, $clip>(
                            cur,
                            cur_min,
                            src,
                            pd.min_score,
                            score_diff,
                            aa_prob[graph.edge_aa[e] as usize],
                        );
                    }
                };
            }
            match (fused, clipped) {
                (true, true) => relax_fused!(true),
                (true, false) => relax_fused!(false),
                (false, true) => relax!(true),
                (false, false) => relax!(false),
            }
        }

        // Lossy trims. Both are applied *after* the row is convolved, not before: the row must
        // exist for its discarded mass to be measured, and trimming still pays for itself because
        // the DP's work is `Σ_nodes retained_width × out-degree` — a narrowed row makes every one
        // of this node's ~21 successors cheaper, and propagates on through `cur_min`/`cur_max`
        // without any further bookkeeping.
        let (mut min_score, mut start, mut len) = (cur_min, start, len);

        // Top cap. Above the threshold `Q_m(r) <= 1` — a path that has already cleared `cut` needs
        // nothing more from its suffix — so a discarded cell can add at most its own probability
        // to the tail, and `err_bound` accumulates exactly that.
        if PRUNED && ceiling < i32::MAX {
            let keep = (ceiling + 1 - min_score).max(1) as usize;
            if keep < len {
                let dropped: f64 = sc.arena[start + keep..start + len].iter().sum();
                err_bound += dropped;
                len = keep;
            }
        }

        // Low end: discard while the *contribution* the cells could still make stays in budget.
        if PRUNED && budget > 0.0 {
            let b = sc.tilt.b[i];
            if b > 0.0 {
                // Spread what is left of the budget over the nodes still to come, so an early node
                // with a fat low tail cannot spend it all.
                let allow = (budget / (n - i) as f64).min(budget);
                let mut acc = 0.0;
                let mut k = 0usize;
                while k + 1 < len {
                    let c = sc.arena[start + k] * b * etab(sc, min_score + k as i32 - cut, theta);
                    if acc + c > allow {
                        break;
                    }
                    acc += c;
                    k += 1;
                }
                if k > 0 {
                    err_bound += acc;
                    budget -= acc;
                    min_score += k as i32;
                    start += k;
                    len -= k;
                }
            }
        }
        sc.cells += len;
        sc.dists[i] = NodeDist {
            min_score,
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
    Some(GenFunc {
        dist,
        valid_from: if PRUNED { cut + credit } else { i32::MIN },
        err_bound,
        denovo: if PRUNED { denovo + credit } else { i32::MIN },
    })
}

#[inline]
fn credit_of(cleavage: Option<Cleavage>) -> i32 {
    cleavage.map_or(0, |c| c.credit.max(c.penalty))
}

/// `e^{θ·d}` for the score offsets `d = s − cut` the trim scan reads. Built once per graph over
/// `[−denovo − 1, 1]` — the trim only ever walks cells below `cut`, and `cut <= denovo` bounds how
/// far below they can sit. Anything outside falls back to a direct `exp`, which keeps the bound
/// valid at the cost of one call.
fn build_etab(sc: &mut DpScratch, theta: f64, denovo: i32) {
    const HI: i32 = 1;
    sc.etab.clear();
    let lo = -(denovo.max(0).saturating_add(1));
    if (HI - lo) as i64 > 1 << 16 {
        sc.etab_lo = 0; // implausible score range; fall back to direct `exp` per lookup
        return;
    }
    sc.etab_lo = lo;
    sc.etab.extend((lo..=HI).map(|d| (theta * d as f64).exp()));
}

#[inline]
fn etab(sc: &DpScratch, d: i32, theta: f64) -> f64 {
    let j = d - sc.etab_lo;
    if j >= 0 && (j as usize) < sc.etab.len() {
        sc.etab[j as usize]
    } else {
        (theta * d as f64).exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-12, "{a} vs {b}");
    }

    /// A deterministic pseudo-random layered DAG shaped like a de novo graph: nodes in topological
    /// order, each drawing edges from a window of predecessors, with signed node/edge scores (the
    /// signs matter — they are why partial scores are not monotone).
    fn random_graph(seed: u64, n: usize) -> Graph {
        let mut s = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        let mut next = move || {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (s >> 33) as i64
        };
        let nodes: Vec<AdjNode> = (0..n)
            .map(|i| {
                let ns = if i == 0 { 0 } else { (next() % 21) as i32 - 10 };
                let mut edges = Vec::new();
                if i > 0 {
                    for back in 1..=6usize.min(i) {
                        if next() % 3 == 0 {
                            continue;
                        }
                        edges.push((i - back, (next() % 17) as i32 - 8, 0.05));
                    }
                }
                (ns, edges)
            })
            .collect();
        Graph::from_adj(&nodes)
    }

    /// Like [`random_graph`] but with in-degrees straddling [`MAX_FUSED_IN`], so both the fused
    /// edge pass and its wide-node fallback run in the same graph.
    fn wide_random_graph(seed: u64, n: usize) -> Graph {
        let mut s = seed.wrapping_mul(6364136223846793005).wrapping_add(7);
        let mut next = move || {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (s >> 33) as i64
        };
        let nodes: Vec<AdjNode> = (0..n)
            .map(|i| {
                let ns = if i == 0 { 0 } else { (next() % 21) as i32 - 10 };
                let mut edges = Vec::new();
                if i > 0 {
                    for back in 1..=40usize.min(i) {
                        if next() % 8 == 0 {
                            continue;
                        }
                        edges.push((i - back, (next() % 17) as i32 - 8, 0.05));
                    }
                }
                (ns, edges)
            })
            .collect();
        Graph::from_adj(&nodes)
    }

    const CLEAVE: Cleavage = Cleavage {
        credit: 2,
        penalty: -11,
        prob_cleavage_sites: 0.1,
    };

    fn assert_same_gf(a: &GenFunc, b: &GenFunc, what: &str) {
        assert_eq!(a.dist.min_score, b.dist.min_score, "{what}: min_score");
        assert_eq!(a.dist.probs.len(), b.dist.probs.len(), "{what}: len");
        for (i, (p, q)) in a.dist.probs.iter().zip(&b.dist.probs).enumerate() {
            assert_eq!(p.to_bits(), q.to_bits(), "{what}: cell {i}");
        }
        assert_eq!(a.denovo, b.denovo, "{what}: denovo");
        assert_eq!(a.valid_from, b.valid_from, "{what}: valid_from");
        assert_eq!(a.err_bound.to_bits(), b.err_bound.to_bits(), "{what}: err");
    }

    /// Caching each edge's resolved descriptor during the score-range pass removes the convolution
    /// pass's second `NodeDist` gather. It must remove **only loads**: the same edges, in the same
    /// order, with the same `score_diff` and `prob`, hence bit-identical `f64` output. Compared
    /// against the pre-fusion two-pass body (`FUSE = false`), which is kept callable for exactly
    /// this. Covers unpruned and pruned, both in-degree regimes, both cleavage settings.
    #[test]
    fn fused_edge_pass_is_bit_exact() {
        let (mut a, mut b) = (DpScratch::default(), DpScratch::default());
        let mut saw_wide = false;
        let mut saw_narrow = false;
        for seed in 0..25u64 {
            for wide in [false, true] {
                let g = if wide {
                    wide_random_graph(seed, 90)
                } else {
                    random_graph(seed, 90)
                };
                let sinks = [89usize];
                for i in 1..g.n_nodes() {
                    let deg = g.edge_start[i + 1] - g.edge_start[i];
                    if deg as usize > MAX_FUSED_IN {
                        saw_wide = true;
                    } else if deg > 0 {
                        saw_narrow = true;
                    }
                }
                for cleavage in [None, Some(CLEAVE)] {
                    let x = compute_inner::<false>(&mut a, &g, &sinks, cleavage, None);
                    let y = compute_inner_unfused::<false>(&mut b, &g, &sinks, cleavage, None);
                    match (&x, &y) {
                        (Some(x), Some(y)) => assert_same_gf(x, y, &format!("full seed {seed}")),
                        (None, None) => continue,
                        _ => panic!("seed {seed}: reachability disagreed"),
                    }
                    let denovo = x.as_ref().unwrap().max_score();
                    for t in [denovo - 30, denovo - 12, denovo - 3, denovo + 1] {
                        let p = Some(Prune::exact(t));
                        let px = compute_inner::<true>(&mut a, &g, &sinks, cleavage, p);
                        let py = compute_inner_unfused::<true>(&mut b, &g, &sinks, cleavage, p);
                        match (&px, &py) {
                            (Some(px), Some(py)) => {
                                assert_same_gf(px, py, &format!("pruned seed {seed} t {t}"))
                            }
                            (None, None) => {}
                            _ => panic!("seed {seed} t {t}: reachability disagreed"),
                        }
                    }
                }
            }
        }
        assert!(saw_narrow, "no node exercised the fused path");
        assert!(saw_wide, "no node exercised the wide-node fallback");
    }

    /// Threshold pruning must reproduce the unpruned tail **bit for bit**, not approximately.
    #[test]
    fn tail_prune_is_bit_exact() {
        let mut sc = DpScratch::default();
        for seed in 0..40u64 {
            let g = random_graph(seed, 60);
            let sinks = [59usize];
            for cleavage in [None, Some(CLEAVE)] {
                let Some(full) = compute(&g, &sinks, cleavage) else {
                    continue;
                };
                let denovo = full.max_score();
                for t in (denovo - 30)..=(denovo + 2) {
                    let pruned =
                        compute_tail_into(&mut sc, &g, &sinks, cleavage, t).expect("reachable");
                    assert_eq!(
                        pruned.spectral_probability(t).to_bits(),
                        full.spectral_probability(t).to_bits(),
                        "seed {seed} threshold {t}"
                    );
                    assert_eq!(pruned.max_score(), denovo, "DeNovoScore, seed {seed}");
                    assert_eq!(pruned.err_bound, 0.0);
                }
            }
        }
    }

    /// A lossy [`Prune`] may only ever *remove* probability, and never by more than the bound it
    /// reports. DeNovoScore must survive even a cap that removes the top of every row.
    #[test]
    fn capped_prune_is_one_sided_and_certified() {
        let mut sc = DpScratch::default();
        for seed in 0..40u64 {
            let g = random_graph(seed, 60);
            let sinks = [59usize];
            let Some(full) = compute(&g, &sinks, Some(CLEAVE)) else {
                continue;
            };
            let denovo = full.max_score();
            for t in [denovo - 25, denovo - 10] {
                let exact = full.spectral_probability(t);
                for cap in [2, 8, 25] {
                    let gf =
                        compute_tail_with(&mut sc, &g, &sinks, Some(CLEAVE), Prune::capped(t, cap))
                            .expect("reachable");
                    let p = gf.spectral_probability(t);
                    assert!(
                        p <= exact * (1.0 + 1e-12),
                        "seed {seed} cap {cap}: {p} > {exact}"
                    );
                    assert!(
                        exact <= p + gf.err_bound * (1.0 + 1e-9) + 1e-300,
                        "seed {seed} cap {cap}: {exact} outside [{p}, {}]",
                        p + gf.err_bound
                    );
                    assert_eq!(gf.max_score(), denovo, "DeNovoScore survives the cap");
                }
            }
        }
    }

    /// The AVX and portable DP bodies are two spellings of one arithmetic, not two approximations.
    #[test]
    #[cfg(target_arch = "x86_64")]
    fn avx_matches_scalar_bitwise() {
        if !is_x86_feature_detected!("avx") {
            eprintln!("skip: no AVX");
            return;
        }
        let (mut a, mut b) = (DpScratch::default(), DpScratch::default());
        for seed in 0..25u64 {
            let g = random_graph(seed, 90);
            let sinks = [89usize];
            // Safety: AVX confirmed present above.
            let x = unsafe { dp_avx::<false, true>(&mut a, &g, &sinks, Some(CLEAVE), None) };
            let y = dp::<false, false, true>(&mut b, &g, &sinks, Some(CLEAVE), None);
            match (x, y) {
                (Some(x), Some(y)) => {
                    assert_eq!(x.dist.min_score, y.dist.min_score);
                    assert_eq!(x.dist.probs.len(), y.dist.probs.len());
                    for (i, (p, q)) in x.dist.probs.iter().zip(&y.dist.probs).enumerate() {
                        assert_eq!(p.to_bits(), q.to_bits(), "seed {seed} cell {i}");
                    }
                }
                (None, None) => {}
                _ => panic!("seed {seed}: one path found the sinks unreachable"),
            }
        }
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
