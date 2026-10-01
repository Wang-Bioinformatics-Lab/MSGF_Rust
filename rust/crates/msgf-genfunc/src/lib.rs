//! # msgf-genfunc — the generating function: DeNovoScore and SpecEValue
//!
//! The null score distribution of one prepared spectrum: a generating-function dynamic program
//! over the amino-acid graph of each isotope sink, the neighbour-cleavage mixture, the merge over
//! sinks, DeNovoScore and the spectral probability (SpecEValue) read off the merged tail
//! (Kim, Gupta & Pevzner, J Proteome Res 2008; Kim & Pevzner, Nat Commun 2014).
//!
//! With `min_query = Some(θ)` ([`score_distribution`]) every DP cell that cannot reach a score of
//! θ is skipped; the tail is exact for every score `>= θ`, which is all a rescore or a search
//! ever queries (θ = the lowest RawScore it will look up).
//!
//! ## Provenance
//!
//! Clean-room implementation, written from a functional specification
//! (`docs/cleanroom/SPEC.md`) by an implementer who never saw MS-GF+ source or this repository's
//! prior generating-function code. Originally `rust/src/dda/specprob/null.rs` in
//! DIA_Proteomics_Rust (commit `b4485bb`, merged as `b930875`; MIT OR Apache-2.0, same author),
//! brought here unchanged apart from module paths and moving `Cleavage` to `msgf-scorer` (both
//! crates use it; re-exported here). See `docs/cleanroom/PROVENANCE.md`.

pub use msgf_scorer::Cleavage;
use msgf_scorer::PreparedSpectrum;
use msgf_scorer::{nominal, residue_mass, RESIDUE_ORDER};

/// One residue of the null alphabet (one edge label of the amino-acid graph).
#[derive(Clone, Copy, Debug)]
pub struct AlphabetEntry {
    /// Unmodified residue letter (used for the cleavage test).
    pub letter: u8,
    /// Real mass including any modification.
    pub mass: f64,
    /// Background probability.
    pub prob: f64,
}

/// Everything the null distribution depends on besides the spectrum.
#[derive(Clone, Debug)]
pub struct NullModel {
    /// Edge labels, in the order their contributions are summed.
    pub alphabet: Vec<AlphabetEntry>,
    pub cleavage: Cleavage,
    /// Isotope-error range `(lo, hi)`; sinks are `N0 - hi ..= N0 - lo`.
    pub isotope: (i32, i32),
}

impl NullModel {
    /// The 20 unmodified residues, probability 0.05 each, trypsin, isotope range (0, 1).
    pub fn uniform() -> NullModel {
        NullModel::from_probs(|_| 0.05, &[], Cleavage::trypsin(), (0, 1))
    }

    /// The 20 residues (with the given per-letter mass deltas, e.g. fixed mods) at probabilities
    /// `prob(letter)`, followed by one extra entry per `(letter, delta)` in `extra` (variable
    /// mods) carrying the base letter's probability.
    pub fn from_probs(
        prob: impl Fn(u8) -> f64,
        extra: &[(u8, f64)],
        cleavage: Cleavage,
        isotope: (i32, i32),
    ) -> NullModel {
        NullModel::from_probs_fixed(prob, &[], extra, cleavage, isotope)
    }

    /// As [`from_probs`](Self::from_probs) with fixed-modification deltas folded into the 20
    /// base entries.
    pub fn from_probs_fixed(
        prob: impl Fn(u8) -> f64,
        fixed: &[(u8, f64)],
        extra: &[(u8, f64)],
        cleavage: Cleavage,
        isotope: (i32, i32),
    ) -> NullModel {
        let base = |l: u8| {
            residue_mass(l).unwrap() + fixed.iter().filter(|f| f.0 == l).map(|f| f.1).sum::<f64>()
        };
        let mut alphabet: Vec<AlphabetEntry> = RESIDUE_ORDER
            .iter()
            .map(|&l| AlphabetEntry {
                letter: l,
                mass: base(l),
                prob: prob(l),
            })
            .collect();
        for &(l, d) in extra {
            alphabet.push(AlphabetEntry {
                letter: l,
                mass: base(l) + d,
                prob: prob(l),
            });
        }
        NullModel {
            alphabet,
            cleavage,
            isotope,
        }
    }

    /// Mixture weight: summed background probability of the cleavage letters (each letter once,
    /// from its first alphabet entry).
    pub fn mixture_pi(&self) -> f64 {
        let mut seen: Vec<u8> = Vec::new();
        let mut pi = 0.0;
        for &l in &self.cleavage.sites {
            if seen.contains(&l) {
                continue;
            }
            seen.push(l);
            if let Some(e) = self.alphabet.iter().find(|e| e.letter == l) {
                pi += e.prob;
            }
        }
        pi
    }
}

/// Per-sink detail (debug / conformance checks).
#[derive(Clone, Debug)]
pub struct SinkDetail {
    pub sink: i32,
    /// v(m) for m in 0..=sink.
    pub vertex_weights: Vec<i32>,
    /// Lowest score of the mixed distribution G.
    pub lowest: i32,
    /// G, cell i at score `lowest + i`.
    pub mass: Vec<f64>,
    /// Best path score to the sink (before the terminus mixture).
    pub best_path: i32,
}

/// The merged null distribution of one (spectrum, charge).
#[derive(Clone, Debug)]
pub struct NullTail {
    lowest: i32,
    mass: Vec<f64>,
    best: i32,
    min_query: Option<i32>,
    /// Filled only by [`score_distribution_detailed`].
    pub sinks: Vec<SinkDetail>,
}

impl NullTail {
    /// DeNovoScore: top of the structural support.
    pub fn best_possible(&self) -> i32 {
        self.best
    }
    /// Lowest stored score (the structural bottom when unpruned).
    pub fn lowest(&self) -> i32 {
        self.lowest
    }
    /// Stored cells, cell i at score `lowest() + i`.
    pub fn cells(&self) -> &[f64] {
        &self.mass
    }
    /// SpecEValue: `min(1, sum of mass at scores >= raw)`. 0 above the support. With pruning,
    /// `raw` must be at least the `min_query` the distribution was built for.
    pub fn tail_mass(&self, raw: i32) -> f64 {
        if let Some(q) = self.min_query {
            assert!(raw >= q, "tail_mass({raw}) below the pruning threshold {q}");
        }
        if raw > self.best {
            return 0.0;
        }
        let start = (raw - self.lowest).max(0) as usize;
        let mut s = 0.0f64;
        for &v in self.mass[start.min(self.mass.len())..].iter() {
            s += v;
        }
        s.min(1.0)
    }
}

/// Build the merged null distribution. With `min_query = Some(θ)` cells that cannot reach a
/// score of θ are skipped (exact for every score >= θ). `None` if no sink is reachable.
pub fn score_distribution(
    prep: &PreparedSpectrum,
    null: &NullModel,
    min_query: Option<i32>,
) -> Option<NullTail> {
    build(prep, null, min_query, false)
}

/// As [`score_distribution`], also keeping per-sink vertex weights and mixed distributions.
pub fn score_distribution_detailed(
    prep: &PreparedSpectrum,
    null: &NullModel,
    min_query: Option<i32>,
) -> Option<NullTail> {
    build(prep, null, min_query, true)
}

struct Edge {
    nom: usize,
    mass: f32,
    prob: f64,
    term: i32,
}

fn build(
    prep: &PreparedSpectrum,
    null: &NullModel,
    min_query: Option<i32>,
    detail: bool,
) -> Option<NullTail> {
    let edges: Vec<Edge> = null
        .alphabet
        .iter()
        .filter_map(|a| {
            let n = nominal(a.mass);
            (n >= 1).then(|| Edge {
                nom: n as usize,
                mass: a.mass as f32,
                prob: a.prob,
                term: null.cleavage.term(a.letter),
            })
        })
        .collect();
    let cl = &null.cleavage;
    let (kmin, kmax) = if cl.enabled {
        (cl.credit.min(cl.penalty), cl.credit.max(cl.penalty))
    } else {
        (0, 0)
    };
    let pi = null.mixture_pi();

    let n0 = prep.n0();
    let (ilo, ihi) = null.isotope;
    let mut per_sink: Vec<(i32, i32, Vec<f64>, i32, Vec<i32>)> = Vec::new(); // sink, G lowest, G, best, v
    for p in (n0 - ihi)..=(n0 - ilo) {
        if p <= 0 {
            continue;
        }
        let cut = min_query.map(|q| q - kmax);
        let Some(sd) = sink_dp(prep, &edges, p, cut) else {
            continue;
        };
        let (lo, g) = if cl.enabled {
            let glo = sd.low + kmin;
            let ghi = sd.best + kmax;
            let mut g = vec![0f64; (ghi - glo + 1) as usize];
            let dget = |s: i32| -> f64 {
                let i = s - sd.low;
                if i < 0 || i as usize >= sd.cells.len() {
                    0.0
                } else {
                    sd.cells[i as usize]
                }
            };
            for (i, c) in g.iter_mut().enumerate() {
                let s = glo + i as i32;
                *c = pi * dget(s - cl.credit) + (1.0 - pi) * dget(s - cl.penalty);
            }
            (glo, g)
        } else {
            (sd.low, sd.cells)
        };
        per_sink.push((
            p,
            lo,
            g,
            sd.best,
            if detail { sd.vertex } else { Vec::new() },
        ));
    }
    if per_sink.is_empty() {
        return None;
    }
    let lowest = per_sink.iter().map(|s| s.1).min().unwrap();
    let top = per_sink
        .iter()
        .map(|s| s.1 + s.2.len() as i32 - 1)
        .max()
        .unwrap();
    let mut mass = vec![0f64; (top - lowest + 1) as usize];
    for s in &per_sink {
        let off = (s.1 - lowest) as usize;
        for (i, v) in s.2.iter().enumerate() {
            mass[off + i] += v;
        }
    }
    let best = per_sink.iter().map(|s| s.3 + kmax).max().unwrap();
    let sinks = if detail {
        per_sink
            .into_iter()
            .map(
                |(sink, lowest, mass, best_path, vertex_weights)| SinkDetail {
                    sink,
                    vertex_weights,
                    lowest,
                    mass,
                    best_path,
                },
            )
            .collect()
    } else {
        Vec::new()
    };
    Some(NullTail {
        lowest,
        mass,
        best,
        min_query,
        sinks,
    })
}

struct SinkDist {
    low: i32,
    cells: Vec<f64>,
    best: i32,
    vertex: Vec<i32>,
}

/// Unmixed distribution D_p of one sink. `cut`: drop cells whose best completion is below it.
fn sink_dp(prep: &PreparedSpectrum, edges: &[Edge], p: i32, cut: Option<i32>) -> Option<SinkDist> {
    let pu = p as usize;
    let na = edges.len();
    // Vertex weights (index = suffix nominal mass) and anchor masses.
    let mut v = vec![0i32; pu + 1];
    for m in 1..pu {
        v[m] = prep.vertex(p, p - m as i32);
    }
    let alpha: Vec<f32> = (0..=p).map(|k| prep.alpha(k)).collect();
    // Edge scores e[m * na + a] (into m, label a); valid when m >= nom_a.
    let mut e = vec![0i32; (pu + 1) * na];
    for m in 1..pu {
        for (a, ed) in edges.iter().enumerate() {
            if m >= ed.nom {
                let mp = m - ed.nom;
                let mut s = prep.edge_between(alpha[m], alpha[mp], ed.mass);
                if mp == 0 {
                    s += ed.term;
                }
                e[m * na + a] = s;
            }
        }
    }
    // Forward: reachability and structural support [lo, hi].
    let mut reach = vec![false; pu + 1];
    let mut lo = vec![i32::MAX; pu + 1];
    let mut hi = vec![i32::MIN; pu + 1];
    reach[0] = true;
    lo[0] = 0;
    hi[0] = 0;
    for m in 1..=pu {
        for (a, ed) in edges.iter().enumerate() {
            if m >= ed.nom && reach[m - ed.nom] {
                let mp = m - ed.nom;
                let sh = v[m] + e[m * na + a];
                reach[m] = true;
                lo[m] = lo[m].min(lo[mp] + sh);
                hi[m] = hi[m].max(hi[mp] + sh);
            }
        }
    }
    if !reach[pu] {
        return None;
    }
    // Backward: best completion to the sink.
    let mut rem = vec![i32::MIN; pu + 1];
    rem[pu] = 0;
    for m in (0..pu).rev() {
        if !reach[m] {
            continue;
        }
        for (a, ed) in edges.iter().enumerate() {
            let to = m + ed.nom;
            if to <= pu && rem[to] != i32::MIN {
                let cand = v[to] + e[to * na + a] + rem[to];
                if cand > rem[m] {
                    rem[m] = cand;
                }
            }
        }
    }
    let best = hi[pu];
    let cut = cut.map(|c| c.min(best));
    // Stored range per vertex.
    let mut low = vec![0i32; pu + 1];
    let mut start = vec![0usize; pu + 2];
    let mut total = 0usize;
    for m in 0..=pu {
        start[m] = total;
        if reach[m] && rem[m] != i32::MIN {
            let l = match cut {
                Some(c) => lo[m].max(c - rem[m]),
                None => lo[m],
            };
            low[m] = l;
            if hi[m] >= l {
                total += (hi[m] - l + 1) as usize;
            }
        }
    }
    start[pu + 1] = total;
    let mut arena = vec![0f64; total];
    arena[0] = 1.0; // D_0 = {0: 1}; low[0] = 0 always (the cut never exceeds the best path).
    for m in 1..=pu {
        let (ms, me) = (start[m], start[m + 1]);
        if ms == me {
            continue;
        }
        let lm = low[m];
        for (a, ed) in edges.iter().enumerate() {
            if m < ed.nom {
                continue;
            }
            let mp = m - ed.nom;
            let (ps, pe) = (start[mp], start[mp + 1]);
            if ps == pe {
                continue;
            }
            let sh = v[m] + e[m * na + a];
            // Source cell i (score low[mp] + i) lands at target index low[mp] + i + sh - lm.
            let base = low[mp] + sh - lm;
            let skip = if base < 0 { (-base) as usize } else { 0 };
            let w = ed.prob;
            let (head, tail) = arena.split_at_mut(ms);
            let src = &head[ps..pe];
            let dst = &mut tail[..me - ms];
            if skip >= src.len() {
                continue;
            }
            let off = (base + skip as i32) as usize;
            for (d, s) in dst[off..].iter_mut().zip(&src[skip..]) {
                *d += w * s;
            }
        }
    }
    let cells = arena[start[pu]..start[pu + 1]].to_vec();
    Some(SinkDist {
        low: low[pu],
        cells,
        best,
        vertex: v,
    })
}
