//! Exponentially tilted **backward** sweeps over the de novo graph — the machinery that lets the
//! DP prune by a cell's *contribution* to the tail rather than by whether it can reach the
//! threshold at all.
//!
//! # What a backward sweep computes
//!
//! For node `m`, let a *suffix* be any path `m → sink`, weighted by the product of its edge
//! probabilities and scored by the sum of the score increments it collects. Then
//!
//! ```text
//!   B_θ(m) = Σ_{suffixes of m} weight · e^{θ · score}
//! ```
//!
//! satisfies the same recursion the DP runs, at width 1 instead of the ~100-cell score support:
//! `B_θ(m) = Σ_{edges m→i} p_e · e^{θ·w_e} · B_θ(i)`, with `B_θ(sink) = 1` and
//! `w_e = nodeScore(i) + edgeScore(e)` (the sink's incoming edges carry `edgeScore = 0`, exactly as
//! in the DP). The CSR stores edges by *destination*, so one descending pass over nodes scatters
//! each node's contribution into its predecessors — the same traversal `max_remaining` uses.
//!
//! # Why it is worth a sweep
//!
//! Markov's inequality on `e^{θ·score}` bounds the suffix tail for **any** `θ ≥ 0`:
//!
//! ```text
//!   Q_m(r) = Σ_{suffixes scoring ≥ r} weight  ≤  e^{−θ·r} · B_θ(m)
//! ```
//!
//! so a DP cell `(m, s)` contributes at most `P(m,s) · e^{−θ(T−s)} · B_θ(m)` to the final tail at
//! `T`. That is the bound a [`crate::Prune::tilt`] policy prunes with. Correctness does not
//! depend on `θ` being any particular value — only the *tightness* does — so a stale or approximate
//! `θ` costs speed, never accuracy.
//!
//! `B_θ(0)` is the whole graph's moment generating function, so carrying `dB/dθ` and `d²B/dθ²`
//! through the same sweep yields `K = ln B_θ(0)` and its first two derivatives. [`solve_theta`]
//! Newton-solves `K'(θ) = T` for the saddlepoint — the `θ` that makes the bound tightest — and
//! returns the leading saddlepoint tail estimate, which is what sizes the pruning error budget.

use crate::{Cleavage, Graph};

/// Largest `θ` [`solve_theta`] will return. `e^{θ·w}` with `|w|` in the tens is already extreme at
/// `θ = 5`; the Newton iteration only wants this many when the threshold is unreachable.
const THETA_MAX: f64 = 5.0;

/// Reusable buffers for the tilted sweeps — held inside [`crate::DpScratch`], one per thread.
#[derive(Default)]
pub struct TiltScratch {
    /// `B_θ(m)` for the last backward sweep — read by the DP's per-node trim.
    pub b: Vec<f64>,
    b1: Vec<f64>,
    b2: Vec<f64>,
    /// `e^{θ·w}` over the graph's `w = nodeScore + edgeScore` range, indexed by `w − w_lo`.
    wtab: Vec<f64>,
    w_lo: i32,
    /// Number of sweeps performed since the last [`Self::reset_sweeps`] — a cost counter.
    pub sweeps: u32,
}

/// `K(θ)`, `K'(θ)`, `K''(θ)` of the source-to-sink score distribution, cleavage weighting included.
#[derive(Clone, Copy, Debug)]
pub struct Cumulants {
    pub k: f64,
    pub k1: f64,
    pub k2: f64,
}

/// What [`solve_theta`] found: a tilt to prune with, and the tail estimate that sizes the budget.
#[derive(Clone, Copy, Debug)]
pub struct Saddle {
    /// The tilt. Always `>= 0`, so the Chernoff bound direction is valid.
    pub theta: f64,
    /// Leading saddlepoint estimate of `P(score >= threshold)`. Used **only** to size an error
    /// budget — never reported as an answer — so its few-percent accuracy is ample.
    pub tail_est: f64,
    /// Sweeps this solve consumed.
    pub sweeps: u32,
}

impl TiltScratch {
    pub fn reset_sweeps(&mut self) {
        self.sweeps = 0;
    }

    /// Build the `e^{θ·w}` table for the node range `0..n`. Returns `false` if the score range is
    /// implausibly wide (then the caller should skip tilted pruning rather than allocate).
    fn plan(&mut self, graph: &Graph, n: usize, theta: f64) -> bool {
        let (mut ns_lo, mut ns_hi) = (0i32, 0i32);
        for &s in &graph.node_score[..n] {
            ns_lo = ns_lo.min(s);
            ns_hi = ns_hi.max(s);
        }
        let (mut es_lo, mut es_hi) = (0i32, 0i32);
        for &s in &graph.edge_score[..graph.edge_start[n] as usize] {
            es_lo = es_lo.min(s);
            es_hi = es_hi.max(s);
        }
        let (lo, hi) = (ns_lo + es_lo, ns_hi + es_hi);
        let span = (hi - lo) as usize + 1;
        if span > 1 << 16 {
            return false;
        }
        self.w_lo = lo;
        self.wtab.clear();
        self.wtab
            .extend((lo..=hi).map(|w| (theta * w as f64).exp()));
        true
    }

    /// One backward tilted sweep over nodes `0..n`. Fills [`Self::b`]; with `derivs` it also fills
    /// `dB/dθ` and `d²B/dθ²`, which costs about three times as much per edge.
    ///
    /// Returns `false` if the recursion left the representable range (then no bound is available).
    fn backward(
        &mut self,
        graph: &Graph,
        sinks: &[usize],
        n: usize,
        theta: f64,
        derivs: bool,
    ) -> bool {
        if !self.plan(graph, n, theta) {
            return false;
        }
        self.sweeps += 1;
        let (wtab, w_lo) = (&self.wtab[..], self.w_lo);

        self.b.clear();
        self.b.resize(n, 0.0);
        if derivs {
            self.b1.clear();
            self.b1.resize(n, 0.0);
            self.b2.clear();
            self.b2.resize(n, 0.0);
        }
        for &s in sinks {
            if s < n {
                self.b[s] = 1.0;
            }
        }

        // Descending: every edge out of `i` lands on a node `> i`, so `b[i]` is final when read.
        for i in (1..n).rev() {
            let bi = self.b[i];
            if bi == 0.0 {
                continue; // no sink downstream of this node
            }
            let (b1i, b2i) = if derivs {
                (self.b1[i], self.b2[i])
            } else {
                (0.0, 0.0)
            };
            let ns = graph.node_score[i];
            let is_sink = sinks.contains(&i);
            let (e0, e1) = (
                graph.edge_start[i] as usize,
                graph.edge_start[i + 1] as usize,
            );
            for e in e0..e1 {
                let p = graph.edge_prev[e] as usize;
                let w = ns + if is_sink { 0 } else { graph.edge_score[e] };
                let f = graph.edge_prob[e] * wtab[(w - w_lo) as usize];
                self.b[p] += f * bi;
                if derivs {
                    let wf = w as f64;
                    self.b1[p] += f * (wf * bi + b1i);
                    self.b2[p] += f * (wf * wf * bi + 2.0 * wf * b1i + b2i);
                }
            }
        }
        let b0 = self.b[0];
        b0.is_finite()
            && b0 > 0.0
            && (!derivs || (self.b1[0].is_finite() && self.b2[0].is_finite()))
    }

    /// Cumulants at the source after the last `backward(.., derivs = true)`.
    fn cumulants(&self, cleavage: Option<Cleavage>, theta: f64) -> Cumulants {
        let (b, b1, b2) = (self.b[0], self.b1[0], self.b2[0]);
        let m1 = b1 / b;
        let (mut k, mut k1, mut k2) = (b.ln(), m1, b2 / b - m1 * m1);
        // The cleavage weighting multiplies the distribution's MGF by a scalar, so cumulants add.
        if let Some(w) = cleavage {
            let (cr, pe) = (w.credit as f64, w.penalty as f64);
            let (a, c) = (
                w.prob_cleavage_sites * (theta * cr).exp(),
                (1.0 - w.prob_cleavage_sites) * (theta * pe).exp(),
            );
            let (c0, c1, c2) = (a + c, a * cr + c * pe, a * cr * cr + c * pe * pe);
            k += c0.ln();
            k1 += c1 / c0;
            k2 += c2 / c0 - (c1 / c0) * (c1 / c0);
        }
        Cumulants { k, k1, k2 }
    }
}

/// Newton-solve `K'(θ) = threshold` for the saddlepoint, and return the leading saddlepoint tail
/// estimate `e^{K − θT} / (θ √(2π K''))`.
///
/// `hint` warm-starts the iteration; `θ̂` moves smoothly between neighbouring candidate masses and
/// spectra, so a good hint typically halves the sweeps. `K'` is solved only to `TOL` score units —
/// the tail is stationary in `θ` at the saddlepoint, so a tighter solve buys nothing and the bound
/// this feeds is valid at any `θ ≥ 0` regardless.
///
/// Returns `None` when the sweep overflows or the estimate is not usable; the caller then falls
/// back to exact pruning.
pub fn solve_theta(
    tilt: &mut TiltScratch,
    graph: &Graph,
    sinks: &[usize],
    n: usize,
    threshold: i32,
    cleavage: Option<Cleavage>,
    hint: Option<f64>,
) -> Option<Saddle> {
    const TOL: f64 = 0.25; // score units
    const MAX_STEPS: u32 = 6;
    let sweeps_before = tilt.sweeps;
    let t = threshold as f64;
    let mut theta = hint.unwrap_or(0.4).clamp(0.0, THETA_MAX);
    let mut cum = None;

    for _ in 0..MAX_STEPS {
        if !tilt.backward(graph, sinks, n, theta, true) {
            return None;
        }
        let c = tilt.cumulants(cleavage, theta);
        cum = Some(c);
        if !c.k1.is_finite() || !c.k2.is_finite() || c.k2 <= 0.0 {
            break;
        }
        if (c.k1 - t).abs() <= TOL {
            break;
        }
        // K' is increasing in θ with slope K'' > 0, so plain Newton converges from either side.
        let next = (theta + (t - c.k1) / c.k2).clamp(0.0, THETA_MAX);
        if (next - theta).abs() < 1e-4 {
            theta = next;
            break;
        }
        theta = next;
    }

    let c = cum?;
    if theta <= 0.0 || !c.k2.is_finite() || c.k2 <= 0.0 {
        // Threshold at or below the mean: nothing to tilt towards, and the tail is O(1).
        return Some(Saddle {
            theta: 0.0,
            tail_est: 1.0,
            sweeps: tilt.sweeps - sweeps_before,
        });
    }
    // The `b` array left behind is the one for the *last* θ evaluated, which is the θ returned.
    let log_tail = c.k - theta * t - (theta * (std::f64::consts::TAU * c.k2).sqrt()).ln();
    let tail_est = log_tail.exp();
    if !tail_est.is_finite() || tail_est <= 0.0 {
        return None;
    }
    Some(Saddle {
        theta,
        tail_est: tail_est.min(1.0),
        sweeps: tilt.sweeps - sweeps_before,
    })
}

/// Refresh `B_θ` for a graph at a `θ` that is already known (no derivatives, no Newton) — about a
/// third the cost of a `solve_theta` sweep. This is how a caller reuses one spectrum's saddlepoint
/// across its isotope-error candidates.
pub fn refresh(
    tilt: &mut TiltScratch,
    graph: &Graph,
    sinks: &[usize],
    n: usize,
    theta: f64,
) -> bool {
    theta >= 0.0 && tilt.backward(graph, sinks, n, theta, false)
}
