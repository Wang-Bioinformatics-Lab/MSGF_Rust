//! **Approximate** SpecEValue by saddlepoint inversion — the score axis removed entirely.
//!
//! # What this is
//!
//! The exact DP ([`crate::compute_into`]) costs `nodes × edges × score-support-width`, and the
//! score axis exists only so that one number can be read off the end: the upper tail at the
//! observed RawScore. Evaluating the generating function at a single exponential tilt `θ` collapses
//! that axis to a scalar per node —
//!
//! ```text
//!   Z_θ(m) = Σ_aa e^{θ·(nodeScore(m) + edgeScore)} · p_aa · Z_θ(m − mass(aa))
//! ```
//!
//! — the same recursion at width 1. Carrying `(Z, dZ/dθ, d²Z/dθ²)` through it yields the cumulant
//! function `K = ln Z` and its first two derivatives, and the tail follows from the Lugannani–Rice
//! saddlepoint formula with the lattice continuity correction.
//!
//! # This is not bit-exact, by construction
//!
//! `CLAUDE.md` makes bit-exactness against MS-GF+ the contract, and this module does **not** meet
//! it — it is a different algorithm for the same quantity, not a faster spelling of the same
//! arithmetic. It is therefore **opt-in**: nothing in the search or CLI path calls it. Use
//! [`crate::compute_tail_into`] when the answer must be exact (it is bit-identical to the full DP
//! above its threshold and gets *cheaper* as the match gets better).
//!
//! # Measured accuracy
//!
//! Against the exact DP over the 1,406-spectrum F13 set, at MS-GF+'s own observed RawScores
//! (`cargo run -p msgf-genfunc --example saddlepoint --release`): median log10 ratio `−0.002`,
//! **96.1%** of spectra within `0.05`, **98.6%** within `0.10`, **100%** within `0.30` — at ~3.2×
//! the speed of the exact DP on the same graphs.
//!
//! The error is a smooth function of how far out the tail is, not noise. Measuring it per decade of
//! *relative depth* — `tail / Z(0)`, the tail as a fraction of the distribution's own total mass —
//! gives the method's real shape (`saddlepoint_tracks_the_exact_tail`):
//!
//! | relative depth | median `|log10|` | worst |
//! |---|---|---|
//! | 1e0 – 1e-6 | 0.002 – 0.004 | 0.021 |
//! | 1e-6 – 1e-9 | 0.009 | 0.027 |
//! | 1e-9 – 1e-12 | 0.010 | 0.096 |
//! | 1e-12 – 1e-16 | 0.040 | 0.323 |
//! | beyond 1e-16 | 0.229 | 1.013 |
//!
//! So it is dependable well past where SpecEValues actually live (~1e-6 relative depth) and decays
//! as the threshold approaches the maximum achievable score, where too few paths carry the tail for
//! a normal-based inversion to work. Note the scale: `Z(0)` is the probability of hitting the
//! precursor mass at all, so an *absolute* SpecEValue of 1e-9 is nowhere near 1e-9 relative depth.
//!
//! Use it as a *first tier* that decides which PSMs deserve the exact DP, not as a silent
//! replacement for one.

use crate::{Cleavage, Graph};

/// Reusable scratch for the tilted sweeps — one per thread, like [`crate::DpScratch`].
#[derive(Default)]
pub struct SaddleScratch {
    z: Vec<f64>,
    z1: Vec<f64>,
    z2: Vec<f64>,
    /// `exp(θ·es)` over the graph's edge-score range. The node score factors out of the edge loop
    /// (`e^{θ(ns+es)} = e^{θ·ns}·e^{θ·es}`), keeping this table small and L1-resident.
    tab: Vec<f64>,
    es_lo: i32,
    es_hi: i32,
    /// `exp(θ·ns)` per node.
    nw: Vec<f64>,
    /// Tilted sweeps performed by the last call — a cost counter for benchmarking.
    pub sweeps: u32,
}

/// `K(θ)`, `K'(θ)` and `K''(θ)` of the sink's score distribution, cleavage weighting included.
#[derive(Clone, Copy, Debug)]
pub struct Cumulants {
    pub k: f64,
    pub k1: f64,
    pub k2: f64,
}

/// The result of a saddlepoint tail estimate.
#[derive(Clone, Copy, Debug)]
pub struct SaddleTail {
    /// The estimated spectral probability `P(score >= threshold)`.
    pub p: f64,
    /// The saddlepoint `θ̂`. Feed it back as the next call's `theta_hint`: `θ̂` moves smoothly
    /// between neighbouring candidate masses and spectra, and a warm start halves the sweeps.
    pub theta: f64,
}

impl SaddleScratch {
    /// Per-graph setup, hoisted out of the θ iteration. Must be called before
    /// [`Self::cumulants`] / [`spectral_probability`] whenever the graph or `n` changes —
    /// re-deriving the edge-score range inside every sweep would double its memory traffic.
    pub fn plan(&mut self, graph: &Graph, n: usize) {
        let n = n.min(graph.n_nodes());
        let (mut lo, mut hi) = (0i32, 0i32);
        for &es in &graph.edge_score[..graph.edge_start[n] as usize] {
            lo = lo.min(es);
            hi = hi.max(es);
        }
        self.es_lo = lo;
        self.es_hi = hi;
        for v in [&mut self.z, &mut self.z1, &mut self.z2, &mut self.nw] {
            v.clear();
            v.resize(n, 0.0);
        }
        self.sweeps = 0;
    }

    /// One tilted sweep: the exact DP's traversal with three scalars per node instead of a score
    /// array, so the per-edge cost is constant rather than proportional to the support width.
    pub fn cumulants(
        &mut self,
        graph: &Graph,
        sinks: &[usize],
        n: usize,
        theta: f64,
        cleavage: Option<Cleavage>,
    ) -> Option<Cumulants> {
        self.sweeps += 1;
        self.tab.clear();
        self.tab
            .extend((self.es_lo..=self.es_hi).map(|d| (theta * d as f64).exp()));
        let (tab, es_lo) = (&self.tab[..], self.es_lo);

        self.z[..n].fill(0.0);
        self.z1[..n].fill(0.0);
        self.z2[..n].fill(0.0);
        self.z[0] = 1.0;
        for i in 0..n {
            self.nw[i] = (theta * graph.node_score[i] as f64).exp();
        }

        for i in 1..n {
            let ns = graph.node_score[i];
            // The sink's incoming edges carry errorScore 0, exactly as in the exact DP.
            let is_sink = sinks.contains(&i);
            let wn = self.nw[i];
            let (e0, e1) = (
                graph.edge_start[i] as usize,
                graph.edge_start[i + 1] as usize,
            );
            let (mut a, mut b, mut c) = (0.0f64, 0.0f64, 0.0f64);
            for e in e0..e1 {
                let p = graph.edge_prev[e] as usize;
                let es = if is_sink { 0 } else { graph.edge_score[e] };
                let w = wn * tab[(es - es_lo) as usize] * graph.edge_prob[e];
                let df = (ns + es) as f64;
                let (zp, z1p, z2p) = (self.z[p], self.z1[p], self.z2[p]);
                a += w * zp;
                b += w * (df * zp + z1p);
                c += w * (df * df * zp + 2.0 * df * z1p + z2p);
            }
            self.z[i] = a;
            self.z1[i] = b;
            self.z2[i] = c;
        }

        let (mut z, mut z1, mut z2) = (0.0f64, 0.0f64, 0.0f64);
        for &s in sinks {
            if s < n {
                z += self.z[s];
                z1 += self.z1[s];
                z2 += self.z2[s];
            }
        }
        if !z.is_finite() || !z2.is_finite() || z <= 0.0 {
            return None;
        }
        let (mut k, mut k1, mut k2) = (z.ln(), z1 / z, z2 / z - (z1 / z) * (z1 / z));
        // The cleavage weighting multiplies the distribution's MGF by a scalar, so its cumulants add.
        if let Some(w) = cleavage {
            let (cr, pe) = (w.credit as f64, w.penalty as f64);
            let (a, b) = (
                w.prob_cleavage_sites * (theta * cr).exp(),
                (1.0 - w.prob_cleavage_sites) * (theta * pe).exp(),
            );
            let (c0, c1, c2) = (a + b, a * cr + b * pe, a * cr * cr + b * pe * pe);
            k += c0.ln();
            k1 += c1 / c0;
            k2 += c2 / c0 - (c1 / c0) * (c1 / c0);
        }
        Some(Cumulants { k, k1, k2 })
    }
}

/// `Z(0)[m]` = the total probability of all peptides of nominal mass exactly `m`, given the
/// per-residue background probabilities `aa` (`(nominal mass, probability)` pairs).
///
/// At `θ = 0` every tilt weight is 1 and the recursion collapses to
/// `z[m] = Σ_aa p_aa · z[m − mass(aa)]` — **the spectrum drops out completely**, as does the
/// cleavage factor. So the normalizing constant every graph needs is one table built once per run,
/// not a sweep per graph. That is what makes [`spectral_probability`] affordable: without it the
/// mandatory `θ = 0` sweep would be a third of the cost.
pub fn total_probability_table(aa: &[(i32, f64)], max_mass: usize) -> Vec<f64> {
    let mut z = vec![0.0f64; max_mass + 1];
    z[0] = 1.0;
    for m in 1..=max_mass {
        let mut s = 0.0;
        for &(nominal, prob) in aa {
            let prev = m as i32 - nominal;
            if prev >= 0 {
                s += prob * z[prev as usize];
            }
        }
        z[m] = s;
    }
    z
}

/// `erfc`, Numerical-Recipes rational form (relative error < 1.2e-7). Only used for shallow tails,
/// where that is far more precision than the answer needs; the deep tail takes the asymptotic
/// branch in [`lugannani_rice`].
fn erfc(x: f64) -> f64 {
    const C: [f64; 10] = [
        -1.3026537197817094,
        6.419_697_923_564_902e-1,
        1.9476473204185836e-2,
        -9.561_514_786_808_63e-3,
        -9.46595344482036e-4,
        3.66839497852761e-4,
        4.2523324806907e-5,
        -2.0278578112534e-5,
        -1.624290004647e-6,
        1.303655835580e-6,
    ];
    let z = x.abs();
    let t = 2.0 / (2.0 + z);
    let ty = 4.0 * t - 2.0;
    let (mut d, mut dd) = (0.0f64, 0.0f64);
    for j in (1..10).rev() {
        let tmp = d;
        d = ty * d - dd + C[j];
        dd = tmp;
    }
    let ans = t * (-z * z + 0.5 * (C[0] + ty * d) - dd).exp();
    if x >= 0.0 {
        ans
    } else {
        2.0 - ans
    }
}

/// Lugannani–Rice upper tail of the *normalized* distribution, with the second continuity
/// correction for a unit lattice: `1 − Φ(ŵ) + φ(ŵ)(1/û − 1/ŵ)`, `û = 2 sinh(θ̂/2)√K''`.
///
/// `1 − Φ(ŵ)` and `−φ(ŵ)/ŵ` cancel catastrophically once `ŵ` is large, so beyond `ŵ = 3` their
/// difference is taken from the Mills-ratio expansion instead — which is the regime that matters,
/// since that is where SpecEValues live.
fn lugannani_rice(theta: f64, cu: &Cumulants, k0: f64, t: f64) -> f64 {
    if theta <= 0.0 {
        return 1.0; // at or below the mean — not a tail
    }
    let arg = 2.0 * (theta * t - (cu.k - k0));
    if arg <= 0.0 || cu.k2 <= 0.0 {
        return 1.0;
    }
    let w = arg.sqrt();
    let u = 2.0 * (theta / 2.0).sinh() * cu.k2.sqrt();
    if u <= 0.0 {
        return 1.0;
    }
    let phi = (-0.5 * w * w).exp() / (2.0 * std::f64::consts::PI).sqrt();
    if w > 3.0 {
        let w2 = w * w;
        let mills_minus = -1.0 / (w * w2) + 3.0 / (w * w2 * w2) - 15.0 / (w * w2 * w2 * w2);
        (phi * (mills_minus + 1.0 / u)).max(0.0)
    } else {
        (0.5 * erfc(w / std::f64::consts::SQRT_2) + phi * (1.0 / u - 1.0 / w)).max(0.0)
    }
}

/// How close `K'(θ)` must get to the target score. The tail is *stationary* in θ at the
/// saddlepoint, so a fraction of a score unit costs `O(Δθ²)` in the answer — solving tighter than
/// this only burns sweeps. Measured on F13: identical accuracy to a 1e-9 solve, ~20% faster.
const SOLVE_TOL: f64 = 0.25;

/// Solve `K'(θ) = target` (`K'` is increasing in θ) by Newton, bisecting whenever Newton would
/// leave the bracket.
fn solve(
    sc: &mut SaddleScratch,
    graph: &Graph,
    sinks: &[usize],
    n: usize,
    target: f64,
    cleavage: Option<Cleavage>,
    theta_hint: f64,
) -> Option<(f64, Cumulants)> {
    let mut theta = theta_hint.clamp(-2.0, 4.0);
    let (mut lo, mut hi) = (-4.0f64, 8.0f64);
    for _ in 0..40 {
        let cu = sc.cumulants(graph, sinks, n, theta, cleavage)?;
        let err = cu.k1 - target;
        if err.abs() < SOLVE_TOL {
            return Some((theta, cu));
        }
        if err > 0.0 {
            hi = theta;
        } else {
            lo = theta;
        }
        let next = if cu.k2 > 1e-12 {
            theta - err / cu.k2
        } else {
            f64::NAN
        };
        theta = if next.is_finite() && next > lo && next < hi {
            next
        } else {
            0.5 * (lo + hi)
        };
    }
    let cu = sc.cumulants(graph, sinks, n, theta, cleavage)?;
    Some((theta, cu))
}

/// Approximate `P(score >= threshold)` for `graph` by saddlepoint inversion — see the module docs
/// for the accuracy this does and does not give you.
///
/// `log_z0` is `ln Z(0)` for this sink mass: pass `Some(ln(table[mass]))` from
/// [`total_probability_table`] (the fast path — it is spectrum-independent), or `None` to have it
/// computed by an extra sweep. `theta_hint` seeds the Newton solve; pass the previous call's
/// [`SaddleTail::theta`], or `0.35` cold.
///
/// Call [`SaddleScratch::plan`] first for this `(graph, n)`.
#[allow(clippy::too_many_arguments)]
pub fn spectral_probability(
    sc: &mut SaddleScratch,
    graph: &Graph,
    sinks: &[usize],
    n: usize,
    threshold: i32,
    cleavage: Option<Cleavage>,
    log_z0: Option<f64>,
    theta_hint: f64,
) -> Option<SaddleTail> {
    let t = threshold as f64 - 0.5; // continuity correction for the integer lattice
    let k0 = match log_z0 {
        Some(v) => v,
        None => sc.cumulants(graph, sinks, n, 0.0, cleavage)?.k,
    };
    let (theta, cu) = solve(sc, graph, sinks, n, t, cleavage, theta_hint)?;
    let p = (k0.exp() * lugannani_rice(theta, &cu, k0, t)).min(1.0);
    Some(SaddleTail { p, theta })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{compute, AdjNode, DpScratch};

    /// Deterministic graph with the shape of a real de novo graph. The mass steps are spread like
    /// residue masses relative to the graph's extent, so a source→sink path is ~10–30 edges — the
    /// length of a tryptic peptide. (Unit mass steps would make paths hundreds of residues long,
    /// whose path probabilities underflow `f64` and exercise nothing real.)
    const STEPS: [usize; 6] = [9, 11, 14, 17, 22, 29];

    fn random_graph(n: usize, seed: u64) -> Vec<AdjNode> {
        let mut s = seed | 1;
        let mut rng = move || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            s
        };
        (0..n)
            .map(|i| {
                let ns = if i == 0 { 0 } else { (rng() % 13) as i32 - 6 };
                let edges = STEPS
                    .iter()
                    .filter(|d| i >= **d)
                    .map(|&d| (i - d, (rng() % 11) as i32 - 5, 0.05))
                    .collect();
                (ns, edges)
            })
            .collect()
    }

    /// How far the saddlepoint estimate may stray from the exact DP, per decade of tail depth.
    ///
    /// "Depth" is `tail / Z(0)` — the tail as a fraction of the distribution's own total mass, which
    /// is the only scale-free way to ask how far out we are. (SpecEValue itself is not: `Z(0)` is
    /// the probability of hitting the precursor mass at all, so an absolute `1e-9` can be the 70th
    /// percentile on one graph and the extreme tail on another.)
    ///
    /// The pattern these encode is the method's real character: comfortably inside the project's
    /// `0.05` SpecEValue bar until ~`1e-9`, then degrading as the threshold approaches the maximum
    /// achievable score, where the tail is carried by too few paths for a normal-based inversion.
    /// Real SpecEValues sit around `1e-6` relative depth, well inside the accurate regime.
    const DEPTH_BOUNDS: [(f64, f64); 6] = [
        (1e-2, 0.02),
        (1e-4, 0.03),
        (1e-6, 0.03),
        (1e-9, 0.05),
        (1e-12, 0.15),
        (1e-16, 0.40),
    ];

    /// The saddlepoint tail must track the exact DP across the upper tail, to the accuracy
    /// [`DEPTH_BOUNDS`] records. This is deliberately *not* the project's flat `|log10| <= 0.05`
    /// bar — that bar governs the exact implementation against MS-GF+, and this module is openly an
    /// approximation whose error grows with depth. Asserting the measured shape catches a
    /// regression without encoding a promise the method does not make.
    #[test]
    fn saddlepoint_tracks_the_exact_tail() {
        let cleave = Cleavage {
            credit: 2,
            penalty: -11,
            prob_cleavage_sites: 0.1,
        };
        // (relative depth, signed log10 ratio) for every probed threshold.
        let mut probed: Vec<(f64, f64)> = Vec::new();
        for seed in [1u64, 7, 12345, 999] {
            for cl in [None, Some(cleave)] {
                let g = Graph::from_adj(&random_graph(260, seed));
                let sinks = [259usize];
                let n = 260;
                let exact = compute(&g, &sinks, cl).expect("reachable");
                // The distribution is a *sub*-probability measure (its total is the probability of
                // reaching the sink mass at all), so "how far into the tail is this?" has to be
                // asked relative to that total.
                let total: f64 = exact.dist.probs.iter().sum();
                let mut sc = SaddleScratch::default();
                sc.plan(&g, n);
                let mut hint = 0.35;
                for t in exact.dist.min_score..=exact.max_score() {
                    let e = exact.spectral_probability(t);
                    let depth = e / total;
                    if e < 1e-280 || depth > 0.05 {
                        continue; // not yet in the upper tail, or below f64 range
                    }
                    let got =
                        spectral_probability(&mut sc, &g, &sinks, n, t, cl, None, hint).unwrap();
                    hint = got.theta;
                    probed.push((depth, (got.p / e).log10()));
                }
            }
        }
        assert!(
            probed.len() > 400,
            "expected a broad sweep, got {}",
            probed.len()
        );

        for (floor, bound) in DEPTH_BOUNDS {
            let band: Vec<f64> = probed
                .iter()
                .filter(|(d, _)| *d >= floor)
                .map(|(_, r)| r.abs())
                .collect();
            let worst = band.iter().fold(0.0f64, |m, r| m.max(*r));
            eprintln!(
                "  depth >= {floor:>7.0e}: n={:<4} worst |log10| {worst:.4}  (bound {bound})",
                band.len()
            );
            assert!(
                worst <= bound,
                "saddlepoint error {worst:.4} exceeds {bound} for tails at or above depth {floor:e}"
            );
        }

        let mut all: Vec<f64> = probed.iter().map(|(_, r)| *r).collect();
        all.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let median = all[all.len() / 2];
        assert!(
            median.abs() <= 0.02,
            "saddlepoint is systematically biased: median log10 ratio {median:+.4}"
        );
    }

    /// `Z(0)` is spectrum-independent, so the tabulated value must equal a `θ = 0` sweep of any
    /// graph built over the same alphabet. This is what lets the table replace a sweep.
    #[test]
    fn total_probability_table_matches_a_theta_zero_sweep() {
        // A graph whose edges are exactly "mass step d, probability 0.05" for d in 1..=6, i.e. the
        // alphabet below — but carrying arbitrary node/edge scores, which θ = 0 must ignore.
        let g = Graph::from_adj(&random_graph(120, 31));
        let alphabet: Vec<(i32, f64)> = STEPS.iter().map(|&d| (d as i32, 0.05)).collect();
        let table = total_probability_table(&alphabet, 119);
        let mut sc = SaddleScratch::default();
        sc.plan(&g, 120);
        for sink in [40usize, 77, 119] {
            let swept = sc.cumulants(&g, &[sink], sink + 1, 0.0, None).unwrap().k;
            let tabulated = table[sink].ln();
            assert!(
                (swept - tabulated).abs() < 1e-9,
                "sink {sink}: swept {swept} vs tabulated {tabulated}"
            );
        }
    }

    /// The estimator must not disagree with the exact DP about the total mass of the distribution.
    #[test]
    fn saddlepoint_and_exact_agree_on_total_probability() {
        let g = Graph::from_adj(&random_graph(150, 5));
        let sinks = [149usize];
        let exact = compute(&g, &sinks, None).unwrap();
        let total: f64 = exact.dist.probs.iter().sum();
        let mut sc = SaddleScratch::default();
        sc.plan(&g, 150);
        let k0 = sc.cumulants(&g, &sinks, 150, 0.0, None).unwrap().k;
        assert!(
            (k0.exp() / total - 1.0).abs() < 1e-9,
            "Z(0) {} vs summed exact distribution {total}",
            k0.exp()
        );
        let _ = DpScratch::default();
    }
}
