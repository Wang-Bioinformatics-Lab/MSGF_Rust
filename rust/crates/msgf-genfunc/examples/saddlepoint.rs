//! Experiment: **replace the score axis with a saddlepoint inversion of the cumulant function.**
//!
//! The DP's cost is `nodes x edges x score-support-width`. But the score axis only exists so we can
//! read one number off the end — the upper tail at the observed RawScore. The generating function
//! evaluated at a single exponential tilt collapses that axis to a *scalar* per node:
//!
//! ```text
//!   Z_θ(m) = Σ_aa e^{θ·(nodeScore(m) + edgeScore)} · p_aa · Z_θ(m − mass(aa))
//! ```
//!
//! which is the same recursion at `width = 1`. Propagating `(Z, dZ/dθ, d²Z/dθ²)` gives the
//! cumulant function `K = log Z` and its first two derivatives, and the tail follows from the
//! Lugannani–Rice saddlepoint formula with the lattice (continuity) correction.
//!
//! This example measures, against the exact DP over real F13 spectra: (a) how accurate that is on
//! the `|log10(approx/exact)| <= 0.05` scale the project holds SpecEValue to, and (b) how much
//! faster it is.
//!
//! Run: cargo run -p msgf-genfunc --example saddlepoint --release
//! Needs the gitignored validation/data/ + the F13 golden.

use std::collections::HashMap;
use std::path::PathBuf;
use std::time::{Duration, Instant};

use msgf_chem::{mass, scaling};
use msgf_genfunc::graph::{build_reverse_graph, standard_aa_nominal, Aa, PeptideCleavage};
use msgf_genfunc::{compute_into, merge_group, Cleavage, DpScratch, Graph};
use msgf_scorer::preprocess::preprocess;
use msgf_scorer::scored_spectrum::ScoredSpectrum;

fn repo(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

const CLEAVE: Cleavage = Cleavage {
    credit: 2,
    penalty: -11,
    prob_cleavage_sites: 0.10,
};

// ---------------------------------------------------------------------------------------------
// The tilted DP: one scalar (plus two θ-derivatives) per node instead of a score distribution.
// ---------------------------------------------------------------------------------------------

/// Scratch for the tilted recursion, reused across calls.
#[derive(Default)]
struct TiltScratch {
    z: Vec<f64>,
    z1: Vec<f64>,
    z2: Vec<f64>,
    /// `exp(θ·es)` for every integer *edge* score in `[es_lo, es_lo + tab.len())`. The node score
    /// factors out of the edge loop (`exp(θ(ns+es)) = exp(θ·ns)·exp(θ·es)`), so this table stays
    /// small and L1-resident.
    tab: Vec<f64>,
    es_lo: i32,
    es_hi: i32,
    /// `exp(θ·ns)` per node.
    nw: Vec<f64>,
}

impl TiltScratch {
    /// Per-graph setup, hoisted out of the θ iteration: the edge-score range (scanning it inside
    /// every sweep would double the sweep's memory traffic) and the buffer sizes.
    fn plan(&mut self, graph: &Graph, n: usize) {
        let (mut lo, mut hi) = (0i32, 0i32);
        for &es in &graph.edge_score[..graph.edge_start[n.min(graph.n_nodes())] as usize] {
            lo = lo.min(es);
            hi = hi.max(es);
        }
        self.es_lo = lo;
        self.es_hi = hi;
        self.z.clear();
        self.z.resize(n, 0.0);
        self.z1.clear();
        self.z1.resize(n, 0.0);
        self.z2.clear();
        self.z2.resize(n, 0.0);
        self.nw.clear();
        self.nw.resize(n, 0.0);
    }
}

/// `Z(0)[m]` = the total probability of all peptides whose nominal mass is exactly `m`.
///
/// At `θ = 0` every tilt weight `e^{θ·score}` is 1, so the recursion degenerates to
/// `z[m] = Σ_aa p_aa · z[m − mass(aa)]` — **the spectrum drops out entirely**. The normalizing
/// constant of every graph at every mass therefore comes from one table, built once for the whole
/// run, instead of a per-graph sweep. (The cleavage factor is also 1 at `θ = 0`.)
fn total_probability_table(aa: &[Aa], max_mass: usize) -> Vec<f64> {
    let mut z = vec![0.0f64; max_mass + 1];
    z[0] = 1.0;
    for m in 1..=max_mass {
        let mut s = 0.0;
        for a in aa {
            let prev = m as i32 - a.nominal;
            if prev >= 0 {
                s += a.prob * z[prev as usize];
            }
        }
        z[m] = s;
    }
    z
}

/// `K(θ)`, `K'(θ)`, `K''(θ)` of the sink's score distribution, including the cleavage factor.
#[derive(Clone, Copy, Debug)]
struct Cumulants {
    k: f64,
    k1: f64,
    k2: f64,
}

/// One tilted sweep. Same traversal as the exact DP, but each node carries three scalars rather
/// than a `[min_score, max_score)` array — so the per-edge cost is constant instead of
/// proportional to the score-support width.
fn cumulants(
    sc: &mut TiltScratch,
    graph: &Graph,
    sinks: &[usize],
    n: usize,
    theta: f64,
    cleavage: Option<Cleavage>,
) -> Option<Cumulants> {
    // exp(θ·es) for the (small) edge-score range, so the sweep itself calls no transcendental.
    sc.tab.clear();
    sc.tab
        .extend((sc.es_lo..=sc.es_hi).map(|d| (theta * d as f64).exp()));
    let (tab, es_lo) = (&sc.tab[..], sc.es_lo);

    sc.z[..n].fill(0.0);
    sc.z1[..n].fill(0.0);
    sc.z2[..n].fill(0.0);
    sc.z[0] = 1.0;
    for i in 0..n {
        sc.nw[i] = (theta * graph.node_score[i] as f64).exp();
    }

    for i in 1..n {
        let ns = graph.node_score[i];
        let is_sink = sinks.contains(&i);
        let wn = sc.nw[i];
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
            let (zp, z1p, z2p) = (sc.z[p], sc.z1[p], sc.z2[p]);
            a += w * zp;
            b += w * (df * zp + z1p);
            c += w * (df * df * zp + 2.0 * df * z1p + z2p);
        }
        sc.z[i] = a;
        sc.z1[i] = b;
        sc.z2[i] = c;
    }

    let (mut z, mut z1, mut z2) = (0.0f64, 0.0f64, 0.0f64);
    for &s in sinks {
        z += sc.z[s];
        z1 += sc.z1[s];
        z2 += sc.z2[s];
    }
    if !z.is_finite() || !z2.is_finite() || z <= 0.0 {
        return None;
    }
    // The cleavage weighting multiplies the distribution's MGF by a scalar factor, so its cumulants
    // simply add.
    let (mut k, mut k1, mut k2) = (z.ln(), z1 / z, z2 / z - (z1 / z) * (z1 / z));
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

/// Solve `K'(θ) = target` (K' is increasing in θ). Returns the saddlepoint and its cumulants.
#[allow(clippy::too_many_arguments)]
fn solve_saddlepoint(
    sc: &mut TiltScratch,
    graph: &Graph,
    sinks: &[usize],
    n: usize,
    target: f64,
    cleavage: Option<Cleavage>,
    theta0: f64,
    iters: &mut usize,
) -> Option<(f64, Cumulants)> {
    let mut theta = theta0.clamp(-2.0, 4.0);
    let (mut lo, mut hi) = (-4.0f64, 8.0f64);
    for _ in 0..40 {
        *iters += 1;
        let cu = cumulants(sc, graph, sinks, n, theta, cleavage)?;
        let err = cu.k1 - target;
        // The tail is stationary in θ at the saddlepoint, so an error of a fraction of a score unit
        // in K' costs only O(Δθ²) in the answer. Solving tighter than this just burns sweeps.
        if err.abs() < 0.25 {
            return Some((theta, cu));
        }
        if err > 0.0 {
            hi = theta;
        } else {
            lo = theta;
        }
        // Newton where it stays bracketed, bisection otherwise.
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
    let cu = cumulants(sc, graph, sinks, n, theta, cleavage)?;
    Some((theta, cu))
}

/// `erfc` (Numerical Recipes rational form, relative error < 1.2e-7) — ample for the regime where
/// the deep-tail asymptotic branch below is not used.
fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 2.0 / (2.0 + z);
    let ty = 4.0 * t - 2.0;
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

/// Lugannani–Rice tail with the second continuity correction for a unit lattice:
/// `P(S >= t) ≈ 1 - Φ(ŵ) + φ(ŵ)(1/û - 1/ŵ)`, `û = 2 sinh(θ̂/2)√K''`, evaluated at `t - 1/2`.
/// The `1 - Φ(ŵ)` and `-φ(ŵ)/ŵ` terms cancel catastrophically in the deep tail, so beyond `ŵ = 3`
/// the Mills-ratio expansion of their difference is used instead.
fn lugannani_rice(theta: f64, cu: &Cumulants, k0: f64, t: f64) -> f64 {
    if theta <= 0.0 {
        return 1.0; // at or below the mean: the tail is not a tail
    }
    let kn = cu.k - k0; // cumulants of the *normalized* distribution
    let arg = 2.0 * (theta * t - kn);
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
        // 1-Φ(w) - φ(w)/w = φ(w)·(-1/w³ + 3/w⁵ - 15/w⁷ + ...) — no cancellation.
        let w2 = w * w;
        let mills_minus = -1.0 / (w * w2) + 3.0 / (w * w2 * w2) - 15.0 / (w * w2 * w2 * w2);
        (phi * (mills_minus + 1.0 / u)).max(0.0)
    } else {
        (0.5 * erfc(w / std::f64::consts::SQRT_2) + phi * (1.0 / u - 1.0 / w)).max(0.0)
    }
}

/// SpecEValue by saddlepoint: `Z(0)` (the total probability of peptides at this mass) times the
/// normalized upper tail.
#[allow(clippy::too_many_arguments)]
fn saddlepoint_tail(
    sc: &mut TiltScratch,
    graph: &Graph,
    sinks: &[usize],
    n: usize,
    threshold: i32,
    cleavage: Option<Cleavage>,
    k0: f64,
    theta0: f64,
    iters: &mut usize,
) -> Option<(f64, f64)> {
    let t = threshold as f64 - 0.5; // continuity correction for the integer lattice
    let (theta, cu) = solve_saddlepoint(sc, graph, sinks, n, t, cleavage, theta0, iters)?;
    // `K(0)` comes from the global mass table — no sweep (see `total_probability_table`).
    let p = k0.exp() * lugannani_rice(theta, &cu, k0, t);
    Some((p, theta))
}

// ---------------------------------------------------------------------------------------------

struct Prepared {
    charge: i32,
    parent_mass: f32,
    pep_nominal: i32,
    raw: Vec<(f32, f32)>,
    threshold: Option<i32>,
}

fn golden_rawscore() -> HashMap<i32, i32> {
    let mut out = HashMap::new();
    let Ok(text) = std::fs::read_to_string(repo("validation/golden/iprg2013_F13.tsv")) else {
        return out;
    };
    let mut lines = text.lines();
    let hdr: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
    let (Some(ci), Some(cs)) = (
        hdr.iter().position(|h| *h == "ScanNum"),
        hdr.iter().position(|h| *h == "MSGFScore"),
    ) else {
        return out;
    };
    for l in lines {
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() <= ci.max(cs) {
            continue;
        }
        if let (Ok(scan), Ok(sc)) = (f[ci].parse::<i32>(), f[cs].parse::<i32>()) {
            let e = out.entry(scan).or_insert(i32::MIN);
            *e = (*e).max(sc);
        }
    }
    out
}

fn load() -> Option<(msgf_scorer::ScoringModel, Vec<Prepared>, Vec<Aa>)> {
    let param = repo("validation/data/models/HCD_HighRes_Tryp.param");
    let mgf = repo("validation/data/spectra/F13.mgf");
    if !param.exists() || !mgf.exists() {
        eprintln!("skipped: validation/data absent");
        return None;
    }
    let model = msgf_scorer::read_param_file(&param).unwrap();
    let golden = golden_rawscore();
    let spectra: Vec<Prepared> = msgf_io::read_mgf_file(&mgf)
        .unwrap()
        .into_iter()
        .filter_map(|s| {
            let charge = s.charge?;
            let mz = s.precursor_mz? as f32;
            let parent_mass = mz * charge as f32 - charge as f32 * mass::PROTON as f32;
            let pep_nominal = scaling::nominal_bin(parent_mass - mass::WATER as f32);
            if !(200..=6000).contains(&pep_nominal) {
                return None;
            }
            let scan: Option<i32> = s.scan.as_deref().and_then(|v| v.parse().ok());
            Some(Prepared {
                charge,
                parent_mass,
                pep_nominal,
                threshold: scan.and_then(|sc| golden.get(&sc).copied()),
                raw: s
                    .peaks
                    .iter()
                    .map(|p| (p.mz as f32, p.intensity as f32))
                    .collect(),
            })
        })
        .collect();
    let mut aa: Vec<Aa> = standard_aa_nominal()
        .into_iter()
        .map(|(r, n)| Aa {
            residue: r,
            nominal: n,
            accurate_mass: msgf_chem::residue_mass(r).unwrap() as f32,
            prob: 0.05,
        })
        .collect();
    let m_ox = msgf_chem::residue_mass(b'M').unwrap() as f32 + 15.994915;
    aa.push(Aa {
        residue: b'M',
        nominal: scaling::nominal_bin(m_ox),
        accurate_mass: m_ox,
        prob: 0.05,
    });
    Some((model, spectra, aa))
}

struct Case {
    graph: Graph,
    tables: msgf_scorer::scored_spectrum::SpectrumTables,
    candidates: Vec<i32>,
    threshold: i32,
}

fn main() {
    let Some((model, spectra, aa)) = load() else {
        return;
    };
    let mut cases: Vec<Case> = Vec::new();
    for s in &spectra {
        let Some(threshold) = s.threshold else {
            continue;
        };
        let peaks = preprocess(&model, s.charge, s.parent_mass, &s.raw);
        let scored = ScoredSpectrum::from_ranked_peaks(&model, s.charge, s.parent_mass, peaks);
        let tables = scored.tables(s.pep_nominal);
        let (graph, _) = build_reverse_graph(
            &scored,
            &tables,
            s.pep_nominal,
            &[s.pep_nominal],
            &aa,
            PeptideCleavage::TRYPSIN,
        );
        cases.push(Case {
            graph,
            tables,
            candidates: (s.pep_nominal - 1..=s.pep_nominal)
                .filter(|&p| p > 0)
                .collect(),
            threshold,
        });
    }
    println!("{} F13 spectra with a golden threshold", cases.len());
    let max_mass = cases
        .iter()
        .flat_map(|c| c.candidates.iter())
        .copied()
        .max()
        .unwrap_or(0) as usize;
    let z0 = total_probability_table(&aa, max_mass);

    // ---- exact reference ---------------------------------------------------------------------
    let mut dp = DpScratch::default();
    let mut exact: Vec<f64> = Vec::with_capacity(cases.len());
    for c in cases.iter_mut() {
        let mut gfs = Vec::new();
        for &p in &c.candidates {
            c.graph.recompute_node_scores(&c.tables, p, &[p]);
            if let Some(gf) = compute_into(&mut dp, &c.graph, &[p as usize], Some(CLEAVE)) {
                gfs.push(gf);
            }
        }
        exact.push(
            merge_group(&gfs)
                .map(|g| g.spectral_probability(c.threshold))
                .unwrap_or(0.0),
        );
    }

    // ---- saddlepoint approximation -----------------------------------------------------------
    let mut tilt = TiltScratch::default();
    let mut approx: Vec<f64> = Vec::with_capacity(cases.len());
    let mut total_iters = 0usize;
    let mut theta_warm = 0.35;
    for c in cases.iter_mut() {
        let mut sum = 0.0;
        for &p in &c.candidates {
            c.graph.recompute_node_scores(&c.tables, p, &[p]);
            let n = (p as usize + 1).min(c.graph.n_nodes());
            tilt.plan(&c.graph, n);
            if let Some((v, th)) = saddlepoint_tail(
                &mut tilt,
                &c.graph,
                &[p as usize],
                n,
                c.threshold,
                Some(CLEAVE),
                z0[p as usize].ln(),
                theta_warm,
                &mut total_iters,
            ) {
                sum += v;
                if th > 0.0 {
                    theta_warm = th; // warm-start the next solve; θ̂ moves smoothly
                }
            }
        }
        approx.push(sum);
    }

    // ---- accuracy ------------------------------------------------------------------------------
    let mut ratios: Vec<f64> = Vec::new();
    let mut zero_exact = 0usize;
    for (a, e) in approx.iter().zip(&exact) {
        if *e <= 0.0 || *a <= 0.0 {
            zero_exact += 1;
            continue;
        }
        ratios.push((a / e).log10());
    }
    ratios.sort_by(|x, y| x.partial_cmp(y).unwrap());
    let q = |f: f64| ratios[((ratios.len() as f64 * f) as usize).min(ratios.len() - 1)];
    let within =
        |b: f64| ratios.iter().filter(|r| r.abs() <= b).count() as f64 / ratios.len() as f64;
    println!(
        "\n=== saddlepoint vs exact DP: log10(approx/exact) over {} spectra ===",
        ratios.len()
    );
    println!("  ({zero_exact} skipped: exact or approx tail was 0)");
    println!(
        "  p01 {:+.3}  p10 {:+.3}  median {:+.3}  p90 {:+.3}  p99 {:+.3}   |max| {:.3}",
        q(0.01),
        q(0.10),
        q(0.50),
        q(0.90),
        q(0.99),
        ratios.iter().fold(0.0f64, |m, r| m.max(r.abs()))
    );
    println!(
        "  within the project's SpecEValue bar |log10| <= 0.05: {:.2}%   <= 0.10: {:.2}%   <= 0.30: {:.2}%",
        100.0 * within(0.05),
        100.0 * within(0.10),
        100.0 * within(0.30)
    );
    println!(
        "  saddlepoint solves: {:.1} tilted sweeps per spectrum",
        total_iters as f64 / cases.len() as f64
    );

    // ---- speed -----------------------------------------------------------------------------------
    let t_exact = {
        let mut best = Duration::MAX;
        for _ in 0..3 {
            let t0 = Instant::now();
            for c in cases.iter_mut() {
                let mut gfs = Vec::new();
                for &p in &c.candidates {
                    c.graph.recompute_node_scores(&c.tables, p, &[p]);
                    if let Some(gf) = compute_into(&mut dp, &c.graph, &[p as usize], Some(CLEAVE)) {
                        gfs.push(gf);
                    }
                }
                std::hint::black_box(
                    merge_group(&gfs).map(|g| g.spectral_probability(c.threshold)),
                );
            }
            best = best.min(t0.elapsed());
        }
        best
    };
    let t_saddle = {
        let mut best = Duration::MAX;
        for _ in 0..3 {
            let mut warm = 0.35;
            let mut it = 0;
            let t0 = Instant::now();
            for c in cases.iter_mut() {
                for &p in &c.candidates {
                    c.graph.recompute_node_scores(&c.tables, p, &[p]);
                    let n = (p as usize + 1).min(c.graph.n_nodes());
                    tilt.plan(&c.graph, n);
                    let r = saddlepoint_tail(
                        &mut tilt,
                        &c.graph,
                        &[p as usize],
                        n,
                        c.threshold,
                        Some(CLEAVE),
                        z0[p as usize].ln(),
                        warm,
                        &mut it,
                    );
                    if let Some((_, th)) = r {
                        if th > 0.0 {
                            warm = th;
                        }
                    }
                    std::hint::black_box(r);
                }
            }
            best = best.min(t0.elapsed());
        }
        best
    };
    println!("\n=== speed (DP only; graph + tables prebuilt) ===");
    println!(
        "  exact DP     {:>9.1}ms\n  saddlepoint  {:>9.1}ms   {:.2}x faster",
        t_exact.as_secs_f64() * 1e3,
        t_saddle.as_secs_f64() * 1e3,
        t_exact.as_secs_f64() / t_saddle.as_secs_f64()
    );
}
