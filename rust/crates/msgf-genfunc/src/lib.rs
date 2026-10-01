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
//! crates use it; re-exported here), then made faster on 2026-10-01 from the same spec and this
//! code alone (shared edge table, blocked integer passes, register-blocked convolution, AVX2
//! dispatch; output bit-identical). See `docs/cleanroom/PROVENANCE.md`.

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

/// Working buffers of [`build`], reused across calls on the same thread so a run allocates them
/// once rather than once per (spectrum, charge). Contents never carry over between calls: every
/// buffer is resized and (re)initialised before it is read.
#[derive(Default)]
struct Scratch {
    /// Anchor masses alpha(k), k in [0, pmax].
    alpha: Vec<f32>,
    /// Edge scores, label-major: `e[a * stride + m]` = score of the edge into `m` with label `a`
    /// (valid for `m >= nom_a`), shared by every sink. The sink's own column is patched to 0 while
    /// that sink runs (the edge into a sink carries no score).
    e: Vec<i32>,
    v: Vec<i32>,
    reach: Vec<i32>,
    lo: Vec<i32>,
    hi: Vec<i32>,
    rem: Vec<i32>,
    /// Per vertex: (arena start, stored length, lowest stored score).
    rows: Vec<(usize, u32, i32)>,
    /// Contributing edges of the row being convolved: (source offset, target range, weight).
    desc: Vec<(isize, isize, isize, f64)>,
    /// Holds the sink's edge-table column while it is patched to 0.
    saved: Vec<i32>,
    arena: Vec<f64>,
}

thread_local! {
    static SCRATCH: std::cell::RefCell<Scratch> = std::cell::RefCell::new(Scratch::default());
}

fn build(
    prep: &PreparedSpectrum,
    null: &NullModel,
    min_query: Option<i32>,
    detail: bool,
) -> Option<NullTail> {
    SCRATCH.with(|s| {
        let mut s = s.borrow_mut();
        #[cfg(target_arch = "x86_64")]
        {
            if std::is_x86_feature_detected!("avx2") {
                // SAFETY: the CPU supports AVX2 (checked above).
                return unsafe { build_avx2(prep, null, min_query, detail, &mut s) };
            }
        }
        build_impl(prep, null, min_query, detail, &mut s)
    })
}

/// [`build_impl`] compiled with AVX2 enabled (wider packed adds/multiplies, an inlined `floor`).
/// No FMA is enabled and Rust never contracts `a * b + c`, so every float result is the same
/// IEEE operation sequence as the baseline build: the output is bit-identical.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn build_avx2(
    prep: &PreparedSpectrum,
    null: &NullModel,
    min_query: Option<i32>,
    detail: bool,
    s: &mut Scratch,
) -> Option<NullTail> {
    build_impl(prep, null, min_query, detail, s)
}

#[inline(always)]
fn build_impl(
    prep: &PreparedSpectrum,
    null: &NullModel,
    min_query: Option<i32>,
    detail: bool,
    s: &mut Scratch,
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
    let pmax = n0 - ilo;
    if pmax >= 1 {
        fill_edge_table(prep, &edges, pmax as usize, s);
    }
    let mut per_sink: Vec<(i32, i32, Vec<f64>, i32, Vec<i32>)> = Vec::new(); // sink, G lowest, G, best, v
    for p in (n0 - ihi)..=(n0 - ilo) {
        if p <= 0 {
            continue;
        }
        let cut = min_query.map(|q| q - kmax);
        let Some(sd) = sink_dp(prep, &edges, p, pmax as usize + 1, cut, detail, s) else {
            continue;
        };
        let (lo, g) = if cl.enabled {
            let glo = sd.low + kmin;
            let ghi = sd.best + kmax;
            let mut g = vec![0f64; (ghi - glo + 1) as usize];
            let dget = |sc: i32| -> f64 {
                let i = sc - sd.low;
                if i < 0 || i as usize >= sd.cells.len() {
                    0.0
                } else {
                    sd.cells[i as usize]
                }
            };
            for (i, c) in g.iter_mut().enumerate() {
                let sc = glo + i as i32;
                *c = pi * dget(sc - cl.credit) + (1.0 - pi) * dget(sc - cl.penalty);
            }
            (glo, g)
        } else {
            (sd.low, sd.cells)
        };
        per_sink.push((p, lo, g, sd.best, sd.vertex));
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

/// Anchor masses up to `pmax` and the shared edge-score table: the score of the edge into vertex
/// `m` with label `a` depends only on `m`, `m - nom_a` and the label, never on the sink, so it is
/// computed once for the largest sink and read by every sink.
#[inline(always)]
fn fill_edge_table(prep: &PreparedSpectrum, edges: &[Edge], pmax: usize, s: &mut Scratch) {
    let stride = pmax + 1;
    s.alpha.clear();
    s.alpha.extend((0..=pmax as i32).map(|k| prep.alpha(k)));
    s.e.clear();
    s.e.resize(edges.len() * stride, 0);
    let alpha = &s.alpha[..];
    for (a, ed) in edges.iter().enumerate() {
        if ed.nom > pmax {
            continue;
        }
        let row = &mut s.e[a * stride..(a + 1) * stride];
        for m in ed.nom..=pmax {
            row[m] = prep.edge_between(alpha[m], alpha[m - ed.nom], ed.mass);
        }
        row[ed.nom] += ed.term; // the edge leaving the source (m - nom = 0)
    }
}

struct SinkDist {
    low: i32,
    cells: Vec<f64>,
    best: i32,
    vertex: Vec<i32>,
}

/// `dst[i] += w * src[i]` over the common length: one IEEE multiply and one add per cell, in
/// order, exactly as written (never contracted to an FMA).
#[inline(always)]
fn axpy(dst: &mut [f64], src: &[f64], w: f64) {
    let n = dst.len().min(src.len());
    let (dst, src) = (&mut dst[..n], &src[..n]);
    for i in 0..n {
        dst[i] += w * src[i];
    }
}

/// Cells per register block of the convolution; also the zero pad kept around every arena row.
const PAD: usize = 16;

/// One row of the convolution, register-blocked: for each block of PAD target cells the edges are
/// accumulated in alphabet order into a local block (`acc += w * src`, starting from +0.0) and
/// stored once. Per cell this is the same sequence of IEEE operations as zeroing the row and
/// running each edge's `axpy` in turn; lanes outside an edge's source range read pad zeros and
/// add +0.0 (see the caller's `blocked` condition).
#[inline(always)]
fn convolve_blocked(head: &[f64], dst: &mut [f64], desc: &[(isize, isize, isize, f64)]) {
    let ml = dst.len() as isize;
    let mut c = 0isize;
    while c < ml {
        let mut acc = [0f64; PAD];
        for &(d, ilo, ihi, w) in desc {
            if c + PAD as isize <= ilo || c >= ihi {
                continue;
            }
            let s0 = (d + c) as usize;
            let src: &[f64; PAD] = head[s0..s0 + PAD].try_into().unwrap();
            for t in 0..PAD {
                acc[t] += w * src[t];
            }
        }
        let c0 = c as usize;
        if c + PAD as isize <= ml {
            let out: &mut [f64; PAD] = (&mut dst[c0..c0 + PAD]).try_into().unwrap();
            *out = acc;
        } else {
            for (o, &x) in dst[c0..].iter_mut().zip(&acc) {
                *o = x;
            }
        }
        c += PAD as isize;
    }
}

/// Unmixed distribution D_p of one sink. `cut`: drop cells whose best completion is below it.
///
/// Reachability, the structural support `[lo, hi]` and the best completion `rem` are integer
/// min/max recurrences, so the order their candidates are combined in does not matter; they are
/// evaluated in blocks of `min nom` vertices (no edge spans less, so a block only reads finished
/// vertices), label by label over contiguous slices, which vectorises. The score convolution adds
/// into each cell in alphabet order exactly as the recurrence is written, so its f64 sums are
/// unchanged.
#[inline(always)]
fn sink_dp(
    prep: &PreparedSpectrum,
    edges: &[Edge],
    p: i32,
    stride: usize,
    cut: Option<i32>,
    detail: bool,
    s: &mut Scratch,
) -> Option<SinkDist> {
    let pu = p as usize;
    let na = edges.len();
    // The edge into the sink scores 0: patch the sink's column, restored before returning.
    s.saved.clear();
    for a in 0..na {
        let x = std::mem::replace(&mut s.e[a * stride + pu], 0);
        s.saved.push(x);
    }
    let out = sink_dp_inner(prep, edges, p, stride, cut, detail, s);
    for a in 0..na {
        s.e[a * stride + pu] = s.saved[a];
    }
    out
}

#[inline(always)]
fn sink_dp_inner(
    prep: &PreparedSpectrum,
    edges: &[Edge],
    p: i32,
    stride: usize,
    cut: Option<i32>,
    detail: bool,
    s: &mut Scratch,
) -> Option<SinkDist> {
    let pu = p as usize;
    let Scratch {
        e,
        v,
        reach,
        lo,
        hi,
        rem,
        rows,
        desc,
        arena,
        ..
    } = s;
    let e = &e[..];
    // Vertex weights (index = suffix nominal mass); v[0] = v[p] = 0.
    v.clear();
    v.resize(pu + 1, 0);
    for m in 1..pu {
        v[m] = prep.vertex(p, p - m as i32);
    }
    let v = &v[..];
    let bs = edges.iter().map(|ed| ed.nom).min().unwrap_or(1).max(1);

    // Forward: reachability (mask 0 / -1) and structural support [lo, hi].
    reach.clear();
    reach.resize(pu + 1, 0);
    lo.clear();
    lo.resize(pu + 1, i32::MAX);
    hi.clear();
    hi.resize(pu + 1, i32::MIN);
    reach[0] = -1;
    lo[0] = 0;
    hi[0] = 0;
    let mut b0 = 1usize;
    while b0 <= pu {
        let b1 = (b0 + bs).min(pu + 1);
        let (r_src, r_dst) = reach.split_at_mut(b0);
        let (l_src, l_dst) = lo.split_at_mut(b0);
        let (h_src, h_dst) = hi.split_at_mut(b0);
        for (a, ed) in edges.iter().enumerate() {
            let nom = ed.nom;
            let from = b0.max(nom);
            if from >= b1 {
                continue;
            }
            let n = b1 - from;
            let (s0, d0) = (from - nom, from - b0);
            let rs = &r_src[s0..s0 + n];
            let ls = &l_src[s0..s0 + n];
            let hs = &h_src[s0..s0 + n];
            let rd = &mut r_dst[d0..d0 + n];
            let ld = &mut l_dst[d0..d0 + n];
            let hd = &mut h_dst[d0..d0 + n];
            let vv = &v[from..b1];
            let ee = &e[a * stride + from..a * stride + b1];
            for i in 0..n {
                let r = rs[i];
                let sh = vv[i].wrapping_add(ee[i]);
                let nl = ld[i].min(ls[i].wrapping_add(sh));
                let nh = hd[i].max(hs[i].wrapping_add(sh));
                ld[i] = (nl & r) | (ld[i] & !r);
                hd[i] = (nh & r) | (hd[i] & !r);
                rd[i] |= r;
            }
        }
        b0 = b1;
    }
    if reach[pu] == 0 {
        return None;
    }
    // Backward: best completion to the sink (i32::MIN = cannot reach it).
    rem.clear();
    rem.resize(pu + 1, i32::MIN);
    rem[pu] = 0;
    let mut b1 = pu; // blocks [b0, b1), top down
    while b1 > 0 {
        let b0 = b1.saturating_sub(bs);
        let (r_dst, r_src) = rem.split_at_mut(b1);
        let r_dst = &mut r_dst[b0..];
        for (a, ed) in edges.iter().enumerate() {
            let nom = ed.nom;
            if nom > pu {
                continue;
            }
            let end = b1.min(pu - nom + 1);
            if b0 >= end {
                continue;
            }
            let n = end - b0;
            // m in [b0, end) reads to = m + nom in [b0 + nom, end + nom), all >= b1.
            let t0 = b0 + nom;
            let rs = &r_src[t0 - b1..t0 - b1 + n];
            let vv = &v[t0..t0 + n];
            let ee = &e[a * stride + t0..a * stride + t0 + n];
            let rd = &mut r_dst[..n];
            for i in 0..n {
                let rt = rs[i];
                let cand = vv[i].wrapping_add(ee[i]).wrapping_add(rt);
                let take = (rt != i32::MIN) & (cand > rd[i]);
                rd[i] = if take { cand } else { rd[i] };
            }
        }
        for (r, &ok) in r_dst.iter_mut().zip(&reach[b0..b1]) {
            if ok == 0 {
                *r = i32::MIN;
            }
        }
        b1 = b0;
    }
    let best = hi[pu];
    let cut = cut.map(|c| c.min(best));
    // Stored range per vertex: (arena start, length, lowest score). Rows are laid out in vertex
    // order with PAD zero cells before the first row and after every row.
    rows.clear();
    let mut total = PAD;
    for m in 0..=pu {
        let mut row = (total, 0u32, 0i32);
        if reach[m] != 0 && rem[m] != i32::MIN {
            let l = match cut {
                Some(c) => lo[m].max(c - rem[m]),
                None => lo[m],
            };
            row.2 = l;
            if hi[m] >= l {
                let len = (hi[m] - l + 1) as usize;
                row.1 = len as u32;
                total += len + PAD;
            }
        }
        rows.push(row);
    }
    // The arena only grows; the leading pad is zeroed here and each row's trailing pad when the
    // row is written, so every pad cell reads as +0.0.
    if arena.len() < total {
        arena.resize(total, 0.0);
    }
    arena[..PAD].fill(0.0);
    // D_0 = {0: 1}; low[0] = 0 always (the cut never exceeds the best path).
    let (s0, l0, _) = rows[0];
    debug_assert_eq!((s0, l0), (PAD, 1));
    arena[s0] = 1.0;
    arena[s0 + 1..s0 + 1 + PAD].fill(0.0);
    // Register blocking needs every weight finite and >= 0: then a pad cell contributes
    // w * (+0.0) = +0.0, and adding +0.0 to a sum of non-negative terms changes no bit, so
    // reading past a source row's ends is the same as skipping the missing cells.
    let blocked = edges.iter().all(|ed| ed.prob.is_finite() && ed.prob >= 0.0);
    for m in 1..=pu {
        let (ms, ml, lm) = rows[m];
        let ml = ml as usize;
        if ml == 0 {
            continue;
        }
        let (head, tail) = arena.split_at_mut(ms);
        let (dst, after) = tail.split_at_mut(ml);
        after[..PAD].fill(0.0);
        // Contributing edges in alphabet order: dst index i reads arena[d + i] for i in
        // [ilo, ihi) (the source row's cells), weight w.
        desc.clear();
        for (a, ed) in edges.iter().enumerate() {
            if m < ed.nom {
                continue;
            }
            let (ps, pl, lp) = rows[m - ed.nom];
            let pl = pl as usize;
            if pl == 0 {
                continue;
            }
            let sh = v[m] + e[a * stride + m];
            // Source cell j (score lp + j) lands at target index lp + j + sh - lm.
            let base = (lp + sh - lm) as isize;
            let (ilo, ihi) = (base.max(0), (base + pl as isize).min(ml as isize));
            if ilo >= ihi {
                continue;
            }
            desc.push((ps as isize - base, ilo, ihi, ed.prob));
        }
        if blocked {
            convolve_blocked(head, dst, desc);
        } else {
            dst.fill(0.0);
            for &(d, ilo, ihi, w) in desc.iter() {
                let (ilo, ihi) = (ilo as usize, ihi as usize);
                let s = (d + ilo as isize) as usize;
                axpy(&mut dst[ilo..ihi], &head[s..s + (ihi - ilo)], w);
            }
        }
    }
    let (ps, pl, low_p) = rows[pu];
    let cells = arena[ps..ps + pl as usize].to_vec();
    Some(SinkDist {
        low: low_p,
        cells,
        best,
        vertex: if detail { v.to_vec() } else { Vec::new() },
    })
}
