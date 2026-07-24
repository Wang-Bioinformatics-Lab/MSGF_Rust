//! Profiling harness for the SpecEValue pipeline over the full F13 high-res set.
//!
//! Answers, without a sampling profiler:
//!   - wall-time share of each stage (preprocess / scored spectrum / graph build / DP compute /
//!     merge / tail), single-thread, summed over all spectra;
//!   - allocations + bytes per stage (custom counting global allocator) → allocations per spectrum.
//!
//! Run:  cargo run -p msgf-genfunc --example profile --release
//!       cargo run -p msgf-genfunc --example profile --release -- hot   # tight loop for `perf stat`
//!
//! Needs the gitignored validation/data/. This file is a throwaway profiling aid.

use std::alloc::{GlobalAlloc, Layout, System};
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::time::{Duration, Instant};

use msgf_chem::{mass, scaling};
use msgf_genfunc::graph::{build_reverse_graph, standard_aa_nominal, Aa};
use msgf_genfunc::{compute_into, merge_group, Cleavage, DpScratch};
use msgf_scorer::preprocess::preprocess;
use msgf_scorer::scored_spectrum::ScoredSpectrum;

// ---- counting allocator (toggleable so timing passes stay clean) -----------------------------
static COUNT_ON: AtomicBool = AtomicBool::new(false);
static ALLOCS: AtomicU64 = AtomicU64::new(0);
static REALLOCS: AtomicU64 = AtomicU64::new(0);
static BYTES: AtomicU64 = AtomicU64::new(0);

struct Counting;
unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        if COUNT_ON.load(Ordering::Relaxed) {
            ALLOCS.fetch_add(1, Ordering::Relaxed);
            BYTES.fetch_add(l.size() as u64, Ordering::Relaxed);
        }
        System.alloc(l)
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        System.dealloc(p, l)
    }
    unsafe fn realloc(&self, p: *mut u8, l: Layout, new: usize) -> *mut u8 {
        if COUNT_ON.load(Ordering::Relaxed) {
            REALLOCS.fetch_add(1, Ordering::Relaxed);
            BYTES.fetch_add(new as u64, Ordering::Relaxed);
        }
        System.realloc(p, l, new)
    }
}
#[global_allocator]
static GLOBAL: Counting = Counting;

fn snap() -> (u64, u64, u64) {
    (
        ALLOCS.load(Ordering::Relaxed),
        REALLOCS.load(Ordering::Relaxed),
        BYTES.load(Ordering::Relaxed),
    )
}

fn repo(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

struct Prepared {
    charge: i32,
    parent_mass: f32,
    pep_nominal: i32,
    raw: Vec<(f32, f32)>,
}

fn load() -> Option<(msgf_scorer::ScoringModel, Vec<Prepared>, Vec<Aa>)> {
    let param = repo("validation/data/models/HCD_HighRes_Tryp.param");
    let mgf = repo("validation/data/spectra/F13.mgf");
    if !param.exists() || !mgf.exists() {
        eprintln!("skipped: validation/data absent (run validation/fetch_reference_data.sh)");
        return None;
    }
    let model = msgf_scorer::read_param_file(&param).unwrap();
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
            Some(Prepared {
                charge,
                parent_mass,
                pep_nominal,
                raw: s.peaks.iter().map(|p| (p.mz as f32, p.intensity as f32)).collect(),
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

#[derive(Default, Clone, Copy)]
struct Stats {
    pre: Duration,
    scored: Duration,
    tables: Duration,
    graph: Duration,
    compute: Duration,
    merge: Duration,
    tail: Duration,
    nodes: u64,
    edges: u64,
    graphs: u64,
    arena_cells: u64,
    reachable: u64,
}

const NSTAGE: usize = 7;

const CLEAVE: Cleavage = Cleavage { credit: 2, penalty: -11, prob_cleavage_sites: 0.10 };

/// One full pass over all spectra with per-stage timers. `census` accumulates alloc counts.
fn pass(
    model: &msgf_scorer::ScoringModel,
    spectra: &[Prepared],
    aa: &[Aa],
    st: &mut Stats,
    census: Option<&mut [(u64, u64, u64); NSTAGE]>,
) {
    let mut acc = [(0u64, 0u64, 0u64); NSTAGE]; // per-stage (allocs, reallocs, bytes) deltas
    let mut scratch = DpScratch::default(); // reused across all spectra — no per-node alloc
    macro_rules! timed {
        ($idx:expr, $field:ident, $body:block) => {{
            let a = snap();
            let t = Instant::now();
            let r = $body;
            st.$field += t.elapsed();
            let b = snap();
            acc[$idx].0 += b.0 - a.0;
            acc[$idx].1 += b.1 - a.1;
            acc[$idx].2 += b.2 - a.2;
            r
        }};
    }

    for s in spectra {
        let peaks = timed!(0, pre, { preprocess(model, s.charge, s.parent_mass, &s.raw) });
        let scored = timed!(1, scored, {
            ScoredSpectrum::from_ranked_peaks(model, s.charge, s.parent_mass, peaks)
        });
        // Shared per-spectrum tables, built once for the largest candidate mass and reused by both
        // isotope-error graphs.
        let tables = timed!(2, tables, { scored.tables(s.pep_nominal) });
        let mut gfs = Vec::new();
        for p in (s.pep_nominal - 1..=s.pep_nominal).filter(|&p| p > 0) {
            let (g, sinks) =
                timed!(3, graph, { build_reverse_graph(&scored, &tables, p, &[p], aa, 2, -11) });
            st.nodes += g.n_nodes() as u64;
            st.edges += g.n_edges() as u64;
            st.graphs += 1;
            let gf = timed!(4, compute, { compute_into(&mut scratch, &g, &sinks, Some(CLEAVE)) });
            st.arena_cells += scratch.arena_len() as u64;
            st.reachable += scratch.reachable() as u64;
            if let Some(gf) = gf {
                gfs.push(gf);
            }
        }
        let merged = timed!(5, merge, { merge_group(&gfs) });
        timed!(6, tail, {
            std::hint::black_box(merged.map(|g| g.spectral_probability(30)));
        });
    }

    if let Some(out) = census {
        *out = acc;
    }
}

fn main() {
    let Some((model, spectra, aa)) = load() else { return };
    let n = spectra.len();

    // "hot" mode: just hammer the single-thread throughput loop for `perf stat`.
    if std::env::args().nth(1).as_deref() == Some("hot") {
        for _ in 0..20 {
            let mut st = Stats::default();
            pass(&model, &spectra, &aa, &mut st, None);
        }
        return;
    }

    // Warm up.
    let mut warm = Stats::default();
    pass(&model, &spectra, &aa, &mut warm, None);

    // Timing pass (counting off) — take the best of 3 to reduce noise.
    COUNT_ON.store(false, Ordering::Relaxed);
    let mut best: Option<Stats> = None;
    for _ in 0..3 {
        let mut st = Stats::default();
        let t = Instant::now();
        pass(&model, &spectra, &aa, &mut st, None);
        let total = t.elapsed();
        if best.is_none_or(|b| total < b_total(&b)) {
            best = Some(st);
        }
    }
    let st = best.unwrap();

    // Allocation census pass (counting on).
    ALLOCS.store(0, Ordering::Relaxed);
    REALLOCS.store(0, Ordering::Relaxed);
    BYTES.store(0, Ordering::Relaxed);
    COUNT_ON.store(true, Ordering::Relaxed);
    let mut cst = Stats::default();
    let mut census = [(0u64, 0u64, 0u64); NSTAGE];
    pass(&model, &spectra, &aa, &mut cst, Some(&mut census));
    COUNT_ON.store(false, Ordering::Relaxed);

    let stages = [
        "preprocess",
        "scored-spec",
        "spec-tables",
        "graph-build",
        "DP compute",
        "merge-group",
        "final-tail",
    ];
    let times = [
        st.pre, st.scored, st.tables, st.graph, st.compute, st.merge, st.tail,
    ];
    let total: Duration = times.iter().sum();

    println!("\n=== SpecEValue profile — {n} F13 spectra, single-thread, nominal grid ===");
    println!("graphs built: {}  ({:.2} per spectrum)", st.graphs, st.graphs as f64 / n as f64);
    println!(
        "graph size:   {:.0} nodes, {:.0} edges per graph",
        st.nodes as f64 / st.graphs as f64,
        st.edges as f64 / st.graphs as f64
    );
    println!(
        "DP work:      {:.0} reachable nodes/graph, {:.0} dist-cells/graph, {:.1} avg support width\n",
        st.reachable as f64 / st.graphs as f64,
        st.arena_cells as f64 / st.graphs as f64,
        st.arena_cells as f64 / st.reachable.max(1) as f64,
    );

    println!(
        "{:<13} {:>9} {:>7}   {:>12} {:>10} {:>12}",
        "stage", "time", "%", "allocs/spec", "reallocs", "MB total"
    );
    println!("{}", "-".repeat(72));
    for i in 0..NSTAGE {
        let (al, re, by) = census[i];
        println!(
            "{:<13} {:>7.0}ms {:>6.1}%   {:>12.1} {:>10.1} {:>12.1}",
            stages[i],
            times[i].as_secs_f64() * 1e3,
            100.0 * times[i].as_secs_f64() / total.as_secs_f64(),
            al as f64 / n as f64,
            re as f64 / n as f64,
            by as f64 / 1e6,
        );
    }
    let total_allocs: u64 = census.iter().map(|c| c.0).sum();
    let total_reallocs: u64 = census.iter().map(|c| c.1).sum();
    let total_bytes: u64 = census.iter().map(|c| c.2).sum();
    println!("{}", "-".repeat(72));
    println!(
        "{:<13} {:>7.0}ms {:>6.1}%   {:>12.1} {:>10.1} {:>12.1}",
        "TOTAL",
        total.as_secs_f64() * 1e3,
        100.0,
        total_allocs as f64 / n as f64,
        total_reallocs as f64 / n as f64,
        total_bytes as f64 / 1e6,
    );
    println!(
        "\nthroughput: {:.0} spectra/s   |   {:.2} ms/spectrum   |   {:.1}M alloc calls total",
        n as f64 / total.as_secs_f64(),
        total.as_secs_f64() * 1e3 / n as f64,
        (total_allocs + total_reallocs) as f64 / 1e6,
    );
}

fn b_total(s: &Stats) -> Duration {
    s.pre + s.scored + s.tables + s.graph + s.compute + s.merge + s.tail
}
