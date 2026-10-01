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
//!
//! ## Fixed N-terminal label
//!
//! [`NullModel::nterm_delta`] (default 0) puts a fixed peptide-N-terminal modification (a TMT /
//! iTRAQ label, …) on the edges into each sink, i.e. on the N-terminal residue of every null
//! string, as candidates carry it on their first residue. 0 is the unlabelled graph exactly. It
//! extends `SPEC.md` (step 5 of `docs/cleanroom/PROVENANCE.md`); `nterm_reference.rs` holds the
//! plain implementation the fast path is tested against.

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
    /// Fixed peptide-N-terminal modification (Da; 0 = none), e.g. a TMT / iTRAQ label on the free
    /// amine. The edges into the sink (the N-terminal residue of every null string) carry it, as
    /// candidates carry it on their first residue. 0 reproduces the unlabelled graph exactly.
    pub nterm_delta: f64,
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
            nterm_delta: 0.0,
        }
    }

    /// This null model with a fixed peptide-N-terminal modification of `delta` Da (see
    /// [`NullModel::nterm_delta`]).
    pub fn with_nterm_delta(mut self, delta: f64) -> NullModel {
        self.nterm_delta = delta;
        self
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
    /// Row layout as structure-of-arrays (start, length, lowest score), for the vectorised
    /// descriptor pass.
    rs: Vec<i64>,
    rl: Vec<i32>,
    rlo: Vec<i32>,
    /// Descriptors of one chunk of target rows, label-major (`a * CHUNK + i`): source offset and
    /// target range; an empty range (0, 0) marks an edge that contributes nothing.
    cd: Vec<i64>,
    clo: Vec<i32>,
    chi: Vec<i32>,
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
            // MSGF_GF_MAX_ISA (diagnostic: `avx2` or `sse2`) caps the path, so every path can be
            // checked for bit-identity on one machine.
            let cap = isa_cap();
            if cap >= 2 && std::is_x86_feature_detected!("avx512f") {
                // SAFETY: the CPU supports AVX-512F (checked above).
                return unsafe { build_avx512(prep, null, min_query, detail, &mut s) };
            }
            if cap >= 1 && std::is_x86_feature_detected!("avx2") {
                // SAFETY: the CPU supports AVX2 (checked above).
                return unsafe { build_avx2(prep, null, min_query, detail, &mut s) };
            }
        }
        build_impl::<Portable>(prep, null, min_query, detail, &mut s)
    })
}

/// The widest generating-function path allowed by `MSGF_GF_MAX_ISA` (2 = AVX-512, 1 = AVX2,
/// 0 = baseline SSE2); unset or unrecognised = no cap. Read once.
#[cfg(target_arch = "x86_64")]
fn isa_cap() -> u8 {
    static CAP: std::sync::OnceLock<u8> = std::sync::OnceLock::new();
    *CAP.get_or_init(|| match std::env::var("MSGF_GF_MAX_ISA").as_deref() {
        Ok("avx2") => 1,
        Ok("sse2") => 0,
        _ => 2,
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
    build_impl::<kernels::Avx2>(prep, null, min_query, detail, s)
}

/// [`build_impl`] compiled with AVX-512F (8-wide f64). Same argument as [`build_avx2`]: no FMA,
/// no contraction, so the IEEE operation sequence per cell is unchanged.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx2")]
unsafe fn build_avx512(
    prep: &PreparedSpectrum,
    null: &NullModel,
    min_query: Option<i32>,
    detail: bool,
    s: &mut Scratch,
) -> Option<NullTail> {
    build_impl::<kernels::Avx512>(prep, null, min_query, detail, s)
}

#[inline(always)]
fn build_impl<K: Kernel>(
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
    // Edges into the sink (the N-terminal residue) under a fixed N-terminal label: the same labels
    // in the same order as `edges`, each mass shifted by the label. `None` = no label, the regular
    // alphabet lands on the sink.
    let nedges: Option<Vec<Edge>> = (null.nterm_delta != 0.0).then(|| {
        null.alphabet
            .iter()
            .filter(|a| nominal(a.mass) >= 1)
            .map(|a| {
                let m = a.mass + null.nterm_delta;
                Edge {
                    nom: nominal(m).max(1) as usize,
                    mass: m as f32,
                    prob: a.prob,
                    term: null.cleavage.term(a.letter),
                }
            })
            .collect()
    });

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
        let Some(sd) = sink_dp::<K>(
            prep,
            &edges,
            nedges.as_deref(),
            p,
            pmax as usize + 1,
            cut,
            detail,
            s,
        ) else {
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
const PAD: usize = 32;

/// One row of the convolution, register-blocked: for each block of PAD target cells the edges are
/// accumulated in alphabet order into a local block (`acc += w * src`, starting from +0.0) and
/// stored once. Per cell this is the same sequence of IEEE operations as zeroing the row and
/// running each edge's `axpy` in turn; lanes outside an edge's source range read pad zeros and
/// add +0.0 (see the caller's `blocked` condition).
#[cfg_attr(target_arch = "x86_64", allow(dead_code))]
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

/// The register-blocked convolution of one row, per instruction set. Every implementation performs,
/// per target cell, the same IEEE sequence as [`convolve_blocked`]: `acc = +0.0`, then for each
/// contributing edge in alphabet order `acc = acc + (w * src)` (a separate multiply and add, never
/// an FMA), then one store. Only the lane width differs, so all are bit-identical.
trait Kernel {
    /// # Safety
    /// The CPU must support the kernel's target features; `head` must hold every block the
    /// descriptors address (`d + c .. d + c + PAD`), which the row layout guarantees.
    unsafe fn convolve(head: &[f64], dst: &mut [f64], desc: &[(isize, isize, isize, f64)]);
}

/// Scalar Rust (the compiler may or may not vectorise it).
struct Portable;
impl Kernel for Portable {
    #[inline(always)]
    unsafe fn convolve(head: &[f64], dst: &mut [f64], desc: &[(isize, isize, isize, f64)]) {
        #[cfg(target_arch = "x86_64")]
        {
            // SSE2 is part of the x86-64 baseline.
            kernels::sse2(head, dst, desc)
        }
        #[cfg(not(target_arch = "x86_64"))]
        convolve_blocked(head, dst, desc)
    }
}

#[cfg(target_arch = "x86_64")]
mod kernels {
    use super::{Kernel, PAD};
    use std::arch::x86_64::*;

    /// Store a finished block: whole when it fits, else only the cells inside the row.
    #[inline(always)]
    unsafe fn store_tail(dst: &mut [f64], c: usize, acc: &[f64; PAD]) {
        for (o, &x) in dst[c..].iter_mut().zip(acc) {
            *o = x;
        }
    }

    #[inline(always)]
    fn src_ptr(head: &[f64], d: isize, c: isize) -> *const f64 {
        let s0 = (d + c) as usize;
        assert!(s0 + PAD <= head.len());
        // SAFETY: in bounds (asserted).
        unsafe { head.as_ptr().add(s0) }
    }

    /// SSE2 has 16 vector registers, so a PAD block is computed as sub-blocks of 16 cells (8
    /// accumulators each); every sub-block walks the descriptors in alphabet order, so the per-cell
    /// sequence is unchanged.
    #[inline(always)]
    pub unsafe fn sse2(head: &[f64], dst: &mut [f64], desc: &[(isize, isize, isize, f64)]) {
        const SB: usize = 16;
        let ml = dst.len() as isize;
        let mut c = 0isize;
        while c < ml {
            let c0 = c as usize;
            let mut t = [0f64; PAD];
            for h in 0..PAD / SB {
                let ch = c + (h * SB) as isize;
                if ch >= ml {
                    break;
                }
                let mut a = [_mm_setzero_pd(); SB / 2];
                for &(d, ilo, ihi, w) in desc {
                    if ch + SB as isize <= ilo || ch >= ihi {
                        continue;
                    }
                    let p = src_ptr(head, d, c).add(h * SB);
                    let wv = _mm_set1_pd(w);
                    for (k, ak) in a.iter_mut().enumerate() {
                        *ak = _mm_add_pd(*ak, _mm_mul_pd(wv, _mm_loadu_pd(p.add(2 * k))));
                    }
                }
                for (k, ak) in a.iter().enumerate() {
                    _mm_storeu_pd(t.as_mut_ptr().add(h * SB + 2 * k), *ak);
                }
            }
            if c + PAD as isize <= ml {
                dst[c0..c0 + PAD].copy_from_slice(&t);
            } else {
                store_tail(dst, c0, &t);
            }
            c += PAD as isize;
        }
    }

    pub struct Avx2;
    impl Kernel for Avx2 {
        #[inline]
        #[target_feature(enable = "avx2")]
        unsafe fn convolve(head: &[f64], dst: &mut [f64], desc: &[(isize, isize, isize, f64)]) {
            let ml = dst.len() as isize;
            let mut c = 0isize;
            while c < ml {
                let mut a = [_mm256_setzero_pd(); PAD / 4];
                for &(d, ilo, ihi, w) in desc {
                    if c + PAD as isize <= ilo || c >= ihi {
                        continue;
                    }
                    let p = src_ptr(head, d, c);
                    let wv = _mm256_set1_pd(w);
                    for (k, ak) in a.iter_mut().enumerate() {
                        *ak = _mm256_add_pd(*ak, _mm256_mul_pd(wv, _mm256_loadu_pd(p.add(4 * k))));
                    }
                }
                let c0 = c as usize;
                if c + PAD as isize <= ml {
                    let q = dst.as_mut_ptr().add(c0);
                    for (k, ak) in a.iter().enumerate() {
                        _mm256_storeu_pd(q.add(4 * k), *ak);
                    }
                } else {
                    let mut t = [0f64; PAD];
                    for (k, ak) in a.iter().enumerate() {
                        _mm256_storeu_pd(t.as_mut_ptr().add(4 * k), *ak);
                    }
                    store_tail(dst, c0, &t);
                }
                c += PAD as isize;
            }
        }
    }

    pub struct Avx512;
    impl Kernel for Avx512 {
        #[inline]
        #[target_feature(enable = "avx512f,avx2")]
        unsafe fn convolve(head: &[f64], dst: &mut [f64], desc: &[(isize, isize, isize, f64)]) {
            let ml = dst.len() as isize;
            let mut c = 0isize;
            while c < ml {
                let mut a = [_mm512_setzero_pd(); PAD / 8];
                for &(d, ilo, ihi, w) in desc {
                    if c + PAD as isize <= ilo || c >= ihi {
                        continue;
                    }
                    let p = src_ptr(head, d, c);
                    let wv = _mm512_set1_pd(w);
                    for (k, ak) in a.iter_mut().enumerate() {
                        *ak = _mm512_add_pd(*ak, _mm512_mul_pd(wv, _mm512_loadu_pd(p.add(8 * k))));
                    }
                }
                let c0 = c as usize;
                if c + PAD as isize <= ml {
                    let q = dst.as_mut_ptr().add(c0);
                    for (k, ak) in a.iter().enumerate() {
                        _mm512_storeu_pd(q.add(8 * k), *ak);
                    }
                } else {
                    let mut t = [0f64; PAD];
                    for (k, ak) in a.iter().enumerate() {
                        _mm512_storeu_pd(t.as_mut_ptr().add(8 * k), *ak);
                    }
                    store_tail(dst, c0, &t);
                }
                c += PAD as isize;
            }
        }
    }
}

/// Target rows per chunk of the vectorised descriptor pass.
const CHUNK: usize = 64;

/// The blocked convolution of every row `1..=pu`. Descriptors are computed a chunk of target rows
/// at a time, label by label over contiguous slices (branch-free integer arithmetic, so it
/// vectorises), then each row's contributing edges are gathered in alphabet order — the same
/// `(d, ilo, ihi, w)` list, in the same order, as building it row by row — and convolved.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn convolve_rows<K: Kernel>(
    edges: &[Edge],
    v: &[i32],
    e: &[i32],
    stride: usize,
    pu: usize,
    rows: &[(usize, u32, i32)],
    desc: &mut Vec<(isize, isize, isize, f64)>,
    arena: &mut [f64],
    rs: &mut Vec<i64>,
    rl: &mut Vec<i32>,
    rlo: &mut Vec<i32>,
    cd: &mut Vec<i64>,
    clo: &mut Vec<i32>,
    chi: &mut Vec<i32>,
) {
    let na = edges.len();
    rs.clear();
    rl.clear();
    rlo.clear();
    for &(st, len, low) in rows {
        rs.push(st as i64);
        rl.push(len as i32);
        rlo.push(low);
    }
    cd.clear();
    cd.resize(na * CHUNK, 0);
    clo.clear();
    clo.resize(na * CHUNK, 0);
    chi.clear();
    chi.resize(na * CHUNK, 0);
    desc.clear();
    desc.resize(na, (0, 0, 0, 0.0));
    let mut m0 = 1usize;
    while m0 <= pu {
        let m1 = (m0 + CHUNK).min(pu + 1);
        for (a, ed) in edges.iter().enumerate() {
            let nom = ed.nom;
            let o = a * CHUNK;
            let from = m0.max(nom).min(m1);
            clo[o..o + (from - m0)].fill(0);
            chi[o..o + (from - m0)].fill(0);
            if from == m1 {
                continue;
            }
            let k = m1 - from;
            let (ss, sl, slo) = (
                &rs[from - nom..from - nom + k],
                &rl[from - nom..from - nom + k],
                &rlo[from - nom..from - nom + k],
            );
            let (tl, tlo) = (&rl[from..m1], &rlo[from..m1]);
            let (vv, ee) = (&v[from..m1], &e[a * stride + from..a * stride + m1]);
            let off = o + (from - m0);
            let (dd, dlo, dhi) = (
                &mut cd[off..off + k],
                &mut clo[off..off + k],
                &mut chi[off..off + k],
            );
            for j in 0..k {
                let sh = vv[j].wrapping_add(ee[j]);
                // Source cell j (score lp + j) lands at target index lp + j + sh - lm.
                let base = slo[j].wrapping_add(sh).wrapping_sub(tlo[j]);
                let ilo = base.max(0);
                let ihi = (base as i64 + sl[j] as i64).min(tl[j] as i64) as i32;
                let ok = (sl[j] > 0) & (ilo < ihi);
                dd[j] = ss[j] - base as i64;
                dlo[j] = if ok { ilo } else { 0 };
                dhi[j] = if ok { ihi } else { 0 };
            }
        }
        for m in m0..m1 {
            let (ms, ml, _) = rows[m];
            let ml = ml as usize;
            if ml == 0 {
                continue;
            }
            let (head, tail) = arena.split_at_mut(ms);
            let (dst, after) = tail.split_at_mut(ml);
            after[..PAD].fill(0.0);
            let i = m - m0;
            let mut n = 0usize;
            for (a, ed) in edges.iter().enumerate() {
                let x = a * CHUNK + i;
                let (l, h) = (clo[x], chi[x]);
                desc[n] = (cd[x] as isize, l as isize, h as isize, ed.prob);
                n += (l < h) as usize;
            }
            // SAFETY: K's target features are enabled by the dispatching caller.
            unsafe { K::convolve(head, dst, &desc[..n]) };
        }
        m0 = m1;
    }
}

/// Unmixed distribution D_p of one sink. `cut`: drop cells whose best completion is below it.
/// `nedges`: the edges into the sink when a fixed N-terminal label shifts them (same labels and
/// order as `edges`); `None` = the regular alphabet.
///
/// Reachability, the structural support `[lo, hi]` and the best completion `rem` are integer
/// min/max recurrences, so the order their candidates are combined in does not matter; they are
/// evaluated in blocks of `min nom` vertices (no edge spans less, so a block only reads finished
/// vertices), label by label over contiguous slices, which vectorises. The score convolution adds
/// into each cell in alphabet order exactly as the recurrence is written, so its f64 sums are
/// unchanged.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn sink_dp<K: Kernel>(
    prep: &PreparedSpectrum,
    edges: &[Edge],
    nedges: Option<&[Edge]>,
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
    let out = sink_dp_inner::<K>(prep, edges, nedges, p, stride, cut, detail, s);
    for a in 0..na {
        s.e[a * stride + pu] = s.saved[a];
    }
    out
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn sink_dp_inner<K: Kernel>(
    prep: &PreparedSpectrum,
    edges: &[Edge],
    nedges: Option<&[Edge]>,
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
        rs,
        rl,
        rlo,
        cd,
        clo,
        chi,
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
    // With an N-terminal label the regular edges stop below the sink; the shifted ones land on it
    // (after the blocks).
    let last = if nedges.is_some() { pu - 1 } else { pu };
    let mut b0 = 1usize;
    while b0 <= last {
        let b1 = (b0 + bs).min(last + 1);
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
    if let Some(ne) = nedges {
        for (a, ed) in ne.iter().enumerate() {
            if ed.nom > pu || reach[pu - ed.nom] == 0 {
                continue;
            }
            let mp = pu - ed.nom;
            let sh = v[pu] + e[a * stride + pu];
            reach[pu] = -1;
            lo[pu] = lo[pu].min(lo[mp] + sh);
            hi[pu] = hi[pu].max(hi[mp] + sh);
        }
    }
    if reach[pu] == 0 {
        return None;
    }
    // Backward: best completion to the sink (i32::MIN = cannot reach it).
    rem.clear();
    rem.resize(pu + 1, i32::MIN);
    rem[pu] = 0;
    // With an N-terminal label only the shifted edges reach the sink: seed their sources here
    // (a max, so the order against the blocks below does not matter; unreachable sources are
    // reset by the blocks' reach mask) and keep the regular edges off the sink.
    if let Some(ne) = nedges {
        for (a, ed) in ne.iter().enumerate() {
            if ed.nom > pu {
                continue;
            }
            let m = pu - ed.nom;
            let cand = v[pu] + e[a * stride + pu] + rem[pu];
            if cand > rem[m] {
                rem[m] = cand;
            }
        }
    }
    let to_sink = usize::from(nedges.is_none()); // 1: regular edges may land on the sink
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
            let end = b1.min(pu - nom + to_sink);
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
    if blocked {
        // With an N-terminal label the chunked pass stops below the sink, whose row (the only
        // one the shifted edges land on) is convolved on its own below.
        convolve_rows::<K>(
            edges, v, e, stride, last, rows, desc, arena, rs, rl, rlo, cd, clo, chi,
        );
        if let Some(ne) = nedges {
            let (ms, ml, lm) = rows[pu];
            let ml = ml as usize;
            if ml > 0 {
                let (head, tail) = arena.split_at_mut(ms);
                let (dst, after) = tail.split_at_mut(ml);
                after[..PAD].fill(0.0);
                desc.clear();
                for (a, ed) in ne.iter().enumerate() {
                    if pu < ed.nom {
                        continue;
                    }
                    let (ps, pl, lp) = rows[pu - ed.nom];
                    if pl == 0 {
                        continue;
                    }
                    let sh = v[pu] + e[a * stride + pu];
                    let base = (lp + sh - lm) as isize;
                    let (ilo, ihi) = (base.max(0), (base + pl as isize).min(ml as isize));
                    if ilo >= ihi {
                        continue;
                    }
                    desc.push((ps as isize - base, ilo, ihi, ed.prob));
                }
                // SAFETY: K's target features are enabled by the dispatching caller.
                unsafe { K::convolve(head, dst, desc) };
            }
        }
    } else {
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
            let into = match nedges {
                Some(ne) if m == pu => ne,
                _ => edges,
            };
            for (a, ed) in into.iter().enumerate() {
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

#[cfg(test)]
mod nterm_reference;

#[cfg(test)]
mod tests {
    use super::*;
    use msgf_scorer::{bundled, prepare_with_cache, Ms2, ScoreModel};

    const PROTON: f64 = 1.007276467;
    const WATER: f64 = 18.0105646863;

    /// Deterministic synthetic spectrum of `seq` with an N-terminal label `label`: b/y ladders
    /// (1+, plus 2+ for z >= 3), a few isotope peaks and uniform noise (64-bit LCG, `seed`).
    fn synthetic(seq: &[u8], label: f64, z: i32, seed: u64) -> (Vec<(f64, f64)>, f64) {
        let mut st = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let mut rnd = move || {
            st = st
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (st >> 11) as f64 / (1u64 << 53) as f64
        };
        let m: Vec<f64> = seq.iter().map(|&l| residue_mass(l).unwrap()).collect();
        let total: f64 = m.iter().sum::<f64>() + WATER + label;
        let mut peaks = Vec::new();
        let mut b = label;
        for x in &m[..m.len() - 1] {
            b += x;
            peaks.push((b + PROTON, 1e4 + 1e6 * rnd()));
            let y = total - b;
            peaks.push((y + PROTON, 1e4 + 1e6 * rnd()));
            peaks.push((y + PROTON + 1.00335, 1e3 + 1e5 * rnd()));
            if z >= 3 {
                peaks.push(((y + 2.0 * PROTON) / 2.0, 1e3 + 1e5 * rnd()));
            }
        }
        for _ in 0..60 {
            peaks.push((100.0 + (total - 100.0) * rnd(), 1e2 + 5e4 * rnd()));
        }
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        (peaks, (total + z as f64 * PROTON) / z as f64)
    }

    fn assert_same(a: &NullTail, b: &NullTail, what: &str) {
        assert_eq!(a.lowest(), b.lowest(), "{what}: lowest");
        assert_eq!(a.best_possible(), b.best_possible(), "{what}: best");
        let bits = |t: &NullTail| t.cells().iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(a), bits(b), "{what}: cells");
        assert_eq!(a.sinks.len(), b.sinks.len(), "{what}: sinks");
        for (x, y) in a.sinks.iter().zip(&b.sinks) {
            assert_eq!(
                (x.sink, x.lowest, x.best_path),
                (y.sink, y.lowest, y.best_path),
                "{what}: sink"
            );
            assert_eq!(x.vertex_weights, y.vertex_weights, "{what}: vertex weights");
            let mb = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
            assert_eq!(mb(&x.mass), mb(&y.mass), "{what}: sink mass");
        }
    }

    /// The fixed N-terminal label of the fast generating function reproduces, bit for bit, the
    /// plain reference implementation it was specified by (`nterm_reference.rs`), over labels
    /// (TMT, TMTpro, iTRAQ-4, acetyl, a negative delta, none), isotope ranges, pruning cuts,
    /// charges, uniform and composition alphabets, with and without a variable modification.
    #[test]
    fn nterm_label_matches_reference() {
        let model: ScoreModel = bundled::score_model().unwrap();
        let peptides: [&[u8]; 6] = [
            b"LVNELTEFAK",
            b"YLYEIARR",
            b"GDVTAQIALQPALK",
            b"AGFAGDDAPRAVFPSIVGR",
            b"HMTEVVR",
            b"SAMPLER",
        ];
        let labels = [
            229.162932, 304.207146, 144.102063, 42.010565, -17.026549, 0.0,
        ];
        let comp = |l: u8| match l {
            b'L' => 0.0996,
            b'K' | b'R' => 0.056,
            b'W' => 0.012,
            b'C' => 0.023,
            _ => 0.046,
        };
        let nulls = [
            NullModel::from_probs(|_| 0.05, &[], Cleavage::trypsin(), (0, 1)),
            NullModel::from_probs_fixed(
                comp,
                &[(b'C', 57.021464)],
                &[(b'M', 15.994915)],
                Cleavage::trypsin(),
                (0, 1),
            ),
            NullModel::from_probs(|_| 0.05, &[], Cleavage::disabled(), (0, 1)),
        ];
        let mut checked = 0;
        let mut changed = 0;
        for (pi, pep) in peptides.iter().enumerate() {
            for &label in &labels {
                for z in [2, 3] {
                    let (peaks, mz) = synthetic(pep, label, z, (pi * 31 + z as usize) as u64);
                    for iso in [(0, 1), (0, 0), (-1, 2)] {
                        let ms2 = Ms2 {
                            peaks: &peaks,
                            precursor_mz: mz,
                            charge: z,
                        };
                        let Some(prep) = prepare_with_cache(&model, &ms2, (-iso.0).max(0)) else {
                            continue;
                        };
                        for base in &nulls {
                            let mut null = base.clone().with_nterm_delta(label);
                            null.isotope = iso;
                            let full = nterm_reference::ref_build(&prep, &null, None, true);
                            let got = score_distribution_detailed(&prep, &null, None);
                            let what = format!("{} label {label} z {z} iso {iso:?}", pi);
                            match (&full, &got) {
                                (Some(a), Some(b)) => assert_same(a, b, &what),
                                (None, None) => continue,
                                _ => panic!("{what}: reachability differs"),
                            }
                            let best = full.as_ref().unwrap().best_possible();
                            for theta in [best - 3, best - 25, best - 80] {
                                let a = nterm_reference::ref_build(&prep, &null, Some(theta), true)
                                    .unwrap();
                                let b =
                                    score_distribution_detailed(&prep, &null, Some(theta)).unwrap();
                                assert_same(&a, &b, &format!("{what} theta {theta}"));
                                let c = score_distribution(&prep, &null, Some(theta)).unwrap();
                                assert_eq!(
                                    b.tail_mass(theta).to_bits(),
                                    c.tail_mass(theta).to_bits()
                                );
                            }
                            if label != 0.0 {
                                let mut plain = null.clone();
                                plain.nterm_delta = 0.0;
                                let p0 = score_distribution(&prep, &plain, None);
                                if p0.map(|t| t.cells().to_vec())
                                    != got.as_ref().map(|t| t.cells().to_vec())
                                {
                                    changed += 1;
                                }
                            }
                            checked += 1;
                        }
                    }
                }
            }
        }
        assert!(checked > 250, "only {checked} cases built");
        assert!(
            changed > 100,
            "the label changed only {changed} distributions"
        );
    }
}
