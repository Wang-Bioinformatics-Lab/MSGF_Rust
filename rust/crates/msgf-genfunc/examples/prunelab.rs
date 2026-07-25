//! Measures the DP's pruning regimes against each other on the full F13 high-res set:
//!
//! 1. **full** — `compute_into`, the whole score distribution (today's shipping path).
//! 2. **exact** — `compute_tail_into`, cells that provably cannot reach the threshold dropped;
//!    bit-identical above the threshold.
//! 3. **cap +N** — plus a top cap N score units above the threshold. Lossy, one-sided, and
//!    **measured to be a bad trade** (see `ALGORITHMIDEAS.md`): partial scores are not monotone
//!    along a path, so clipping high intermediate scores removes real tail mass.
//! 4. **cap +N +tilt** — plus the Chernoff low-end trim. Buys a few percent of cells for a backward
//!    tilted sweep costing ~13% of the DP. Also a bad trade; kept so the evidence is reproducible.
//!
//! Reports, per regime: retained distribution cells, DP wall-time, tilted sweeps, the certified
//! relative error the run itself accumulated, and the *actual* log10 deviation from the exact tail.
//!
//! Thresholds come either from the frozen MS-GF+ F13 output (the observed top-hit RawScore — what a
//! search knows once it has scored its candidates) or from `DeNovoScore − k`, which sweeps the
//! quality of the match without depending on the golden.
//!
//! Run: cargo run -p msgf-genfunc --example prunelab --release [--] [k1 k2 ...]
//! Needs the gitignored validation/data/ (+ the F13 golden for the observed-RawScore row).

use std::collections::HashMap;
use std::path::PathBuf;
use std::time::{Duration, Instant};

use msgf_chem::{mass, scaling};
use msgf_genfunc::graph::{build_reverse_graph, standard_aa_nominal, Aa, PeptideCleavage};
use msgf_genfunc::{
    compute_into, compute_tail_into, compute_tail_with, merge_group, tilt, Cleavage, DpScratch,
    GenFunc, Graph, Prune,
};
use msgf_scorer::preprocess::preprocess;
use msgf_scorer::scored_spectrum::{ScoredSpectrum, SpectrumTables};

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

struct Prepared {
    charge: i32,
    parent_mass: f32,
    pep_nominal: i32,
    raw: Vec<(f32, f32)>,
    /// MS-GF+'s own top-hit RawScore for this scan, when the golden has one.
    golden_score: Option<i32>,
}

/// Everything for one spectrum that does not depend on the threshold — kept out of the timed loop.
struct Work {
    graph: Graph,
    tables: SpectrumTables,
    candidates: Vec<i32>,
    golden_score: Option<i32>,
    /// DeNovoScore from the reference (unpruned) pass.
    denovo: i32,
    /// Reference tail at each threshold this run measures, from the unpruned DP.
    reference: Vec<f64>,
}

fn load() -> Option<(msgf_scorer::ScoringModel, Vec<Prepared>, Vec<Aa>)> {
    let param = repo("validation/data/models/HCD_HighRes_Tryp.param");
    let mgf = repo("validation/data/spectra/F13.mgf");
    if !param.exists() || !mgf.exists() {
        eprintln!("skipped: validation/data absent (run validation/fetch_reference_data.sh)");
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
                golden_score: scan.and_then(|sc| golden.get(&sc).copied()),
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

/// Best observed RawScore (`MSGFScore`) per scan, from the frozen MS-GF+ F13 output.
fn golden_rawscore() -> HashMap<i32, i32> {
    let mut out = HashMap::new();
    let Ok(text) = std::fs::read_to_string(repo("validation/golden/iprg2013_F13.tsv")) else {
        eprintln!("note: F13 golden absent — the observed-RawScore row is skipped");
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

fn prepare(model: &msgf_scorer::ScoringModel, s: &Prepared, aa: &[Aa]) -> Work {
    let peaks = preprocess(model, s.charge, s.parent_mass, &s.raw);
    let scored = ScoredSpectrum::from_ranked_peaks(model, s.charge, s.parent_mass, peaks);
    let tables = scored.tables(s.pep_nominal);
    let (graph, _) = build_reverse_graph(
        &scored,
        &tables,
        s.pep_nominal,
        &[s.pep_nominal],
        aa,
        PeptideCleavage::TRYPSIN,
    );
    Work {
        graph,
        tables,
        candidates: (s.pep_nominal - 1..=s.pep_nominal)
            .filter(|&p| p > 0)
            .collect(),
        golden_score: s.golden_score,
        denovo: 0,
        reference: Vec::new(),
    }
}

/// Which threshold a run uses for a spectrum. `None` → the spectrum is skipped by that run.
#[derive(Clone, Copy)]
enum Source {
    /// MS-GF+'s own observed top-hit RawScore.
    Observed,
    /// `DeNovoScore − k` — sweeps match quality without needing the golden.
    Below(i32),
}

impl Source {
    fn threshold(&self, w: &Work) -> Option<i32> {
        match self {
            Source::Observed => w.golden_score,
            Source::Below(k) => Some(w.denovo - k),
        }
    }
    fn label(&self) -> String {
        match self {
            Source::Observed => "MS-GF+ RawScore".into(),
            Source::Below(k) => format!("DeNovoScore − {k}"),
        }
    }
}

#[derive(Clone, Copy)]
enum Mode {
    Full,
    Exact,
    /// Exact floor + the forward `max_achievable` sweep that skips dead nodes' edge lists.
    ExactSkip,
    /// Exact floor + a top cap of N score units above the threshold.
    Cap(i32),
    /// Top cap plus the tilted low-end trim at relative budget `eps` — the sweep-paying variant.
    CapTilt(i32, f64),
}

#[derive(Default, Clone, Copy)]
struct Run {
    time: Duration,
    cells: u64,
    written: u64,
    sweeps: u64,
    graphs: u64,
    /// Spectra where the run produced a usable tail.
    n: u64,
    /// Sum / max of |log10(p_run / p_exact)|.
    log_sum: f64,
    log_max: f64,
    /// Sum / max of the *certified* relative error the run accumulated.
    cert_sum: f64,
    cert_max: f64,
    /// Spectra whose certified error exceeded the requested eps.
    over_eps: u64,
    /// Spectra where the pruned tail was not bit-identical to the exact one.
    not_bitexact: u64,
    /// DeNovoScore mismatches against the reference pass.
    denovo_bad: u64,
}

/// One pass over every spectrum for one (mode, threshold source). `ti` indexes `Work::reference`.
fn pass(works: &mut [Work], mode: Mode, src: Source, ti: usize, sc: &mut DpScratch) -> Run {
    let mut r = Run::default();
    let mut gfs: Vec<GenFunc> = Vec::with_capacity(2);
    let mut theta_hint: Option<f64> = None;

    for w in works.iter_mut() {
        let Some(threshold) = src.threshold(w) else {
            continue;
        };
        gfs.clear();
        sc.reset_sweeps();
        let t0 = Instant::now();
        for &p in &w.candidates {
            w.graph.recompute_node_scores(&w.tables, p, &[p]);
            let sinks = [p as usize];
            let gf = match mode {
                Mode::Full => compute_into(sc, &w.graph, &sinks, Some(CLEAVE)),
                Mode::Exact => compute_tail_into(sc, &w.graph, &sinks, Some(CLEAVE), threshold),
                Mode::ExactSkip => {
                    let mut prune = Prune::exact(threshold);
                    prune.skip_dead_nodes = true;
                    compute_tail_with(sc, &w.graph, &sinks, Some(CLEAVE), prune)
                }
                Mode::Cap(cap) => compute_tail_with(
                    sc,
                    &w.graph,
                    &sinks,
                    Some(CLEAVE),
                    Prune::capped(threshold, cap),
                ),
                Mode::CapTilt(cap, eps) => {
                    let n = w.graph.n_nodes().min(p as usize + 1);
                    let cut = threshold - CLEAVE.credit.max(CLEAVE.penalty);
                    let mut prune = Prune::capped(threshold, cap);
                    if let Some(s) = tilt::solve_theta(
                        sc.tilt_mut(),
                        &w.graph,
                        &sinks,
                        n,
                        cut,
                        Some(CLEAVE),
                        theta_hint,
                    ) {
                        theta_hint = Some(s.theta);
                        prune.tilt = Some((s.theta, eps * s.tail_est));
                    }
                    compute_tail_with(sc, &w.graph, &sinks, Some(CLEAVE), prune)
                }
            };
            r.cells += sc.cells() as u64;
            r.written += sc.arena_len() as u64;
            r.graphs += 1;
            if let Some(gf) = gf {
                gfs.push(gf);
            }
        }
        let merged = merge_group(&gfs);
        let p = merged.as_ref().map(|g| g.spectral_probability(threshold));
        r.time += t0.elapsed();
        r.sweeps += sc.sweeps() as u64;

        let Some(gf) = merged else { continue };
        let p = p.unwrap();
        r.n += 1;
        if gf.max_score() != w.denovo {
            r.denovo_bad += 1;
        }

        match mode {
            Mode::Full => {
                // The reference pass: record the exact tail for every threshold this run measures.
                if w.reference.len() <= ti {
                    w.reference.resize(ti + 1, f64::NAN);
                }
                w.reference[ti] = p;
            }
            _ => {
                let exact = w.reference.get(ti).copied().unwrap_or(f64::NAN);
                if exact.is_finite() {
                    if p.to_bits() != exact.to_bits() {
                        r.not_bitexact += 1;
                    }
                    if p > 0.0 && exact > 0.0 {
                        let d = (p / exact).log10().abs();
                        r.log_sum += d;
                        r.log_max = r.log_max.max(d);
                    }
                }
                let cert = gf.relative_error();
                if cert.is_finite() {
                    r.cert_sum += cert;
                    r.cert_max = r.cert_max.max(cert);
                }
                if cert > 1e-3 {
                    r.over_eps += 1;
                }
            }
        }
    }
    r
}

fn main() {
    let Some((model, spectra, aa)) = load() else {
        return;
    };
    let ks: Vec<i32> = std::env::args()
        .skip(1)
        .filter_map(|a| a.parse().ok())
        .collect();
    let ks = if ks.is_empty() {
        vec![5, 10, 20, 40]
    } else {
        ks
    };

    let mut works: Vec<Work> = spectra.iter().map(|s| prepare(&model, s, &aa)).collect();
    let mut sc = DpScratch::default();

    // DeNovoScore for every spectrum (needed by the `DeNovoScore − k` sources).
    for w in works.iter_mut() {
        let mut gfs = Vec::new();
        for &p in &w.candidates {
            w.graph.recompute_node_scores(&w.tables, p, &[p]);
            if let Some(gf) = compute_into(&mut sc, &w.graph, &[p as usize], Some(CLEAVE)) {
                gfs.push(gf);
            }
        }
        w.denovo = merge_group(&gfs).map_or(i32::MIN, |g| g.max_score());
    }
    works.retain(|w| w.denovo > i32::MIN);

    let sources: Vec<Source> = std::iter::once(Source::Observed)
        .chain(ks.iter().map(|&k| Source::Below(k)))
        .collect();
    let mut modes: Vec<(String, Mode)> = vec![
        ("full".into(), Mode::Full),
        ("exact".into(), Mode::Exact),
        ("exact +skip".into(), Mode::ExactSkip),
    ];
    for cap in [40, 30, 20] {
        modes.push((format!("cap +{cap}"), Mode::Cap(cap)));
    }
    modes.push(("cap +30 +tilt".into(), Mode::CapTilt(30, 1e-3)));

    println!(
        "\n=== DP pruning lab — {} F13 spectra, single-thread, nominal grid ===",
        works.len()
    );

    for (ti, src) in sources.iter().enumerate() {
        // Warm + reference pass first: `Mode::Full` fills `Work::reference[ti]`.
        let mut base: Option<Run> = None;
        let mut rows: Vec<(String, Run)> = Vec::new();
        for (name, mode) in &modes {
            // Best of 2 for the timing, after one warm pass.
            let mut best: Option<Run> = None;
            for _ in 0..3 {
                let r = pass(&mut works, *mode, *src, ti, &mut sc);
                if best.is_none_or(|b: Run| r.time < b.time) {
                    best = Some(r);
                }
            }
            let r = best.unwrap();
            if matches!(mode, Mode::Full) {
                base = Some(r);
            }
            rows.push((name.clone(), r));
        }
        let base = base.unwrap();

        println!(
            "\n--- threshold = {} ({} spectra scored) ---",
            src.label(),
            base.n
        );
        println!(
            "{:<14} {:>9} {:>8} {:>10} {:>8} {:>10} {:>10} {:>9} {:>7}",
            "mode",
            "DP time",
            "speedup",
            "cells/gf",
            "vs full",
            "med |log10|",
            "max |log10|",
            "cert rel",
            "cert>1e-3"
        );
        println!("{}", "-".repeat(96));
        for (name, r) in &rows {
            let cells = r.cells as f64 / r.graphs.max(1) as f64;
            let bits = if r.not_bitexact > 0 {
                format!("{:>9}", r.not_bitexact)
            } else {
                format!("{:>9.1e}", r.cert_sum / r.n.max(1) as f64)
            };
            println!(
                "{:<14} {:>7.0}ms {:>7.2}× {:>10.0} {:>7.2}× {:>10.4} {:>10.4} {} {:>7}",
                name,
                r.time.as_secs_f64() * 1e3,
                base.time.as_secs_f64() / r.time.as_secs_f64(),
                cells,
                (base.cells as f64 / base.graphs.max(1) as f64) / cells,
                r.log_sum / r.n.max(1) as f64,
                r.log_max,
                bits,
                r.over_eps,
            );
        }
        let denovo_bad: u64 = rows.iter().map(|(_, r)| r.denovo_bad).sum();
        let tilt_row = rows.last().unwrap();
        println!(
            "DeNovoScore mismatches: {denovo_bad}   |   tilted sweeps ({}): {:.2} per graph   |   \
             exact-prune bit-exact: {}",
            tilt_row.0,
            tilt_row.1.sweeps as f64 / tilt_row.1.graphs.max(1) as f64,
            if rows[1].1.not_bitexact == 0 {
                "yes"
            } else {
                "NO"
            },
        );
    }
    println!(
        "\nnote: 'med |log10|' is the mean absolute log10 ratio against the unpruned DP; \
         'cert rel' is the run's own accumulated error bound (mean), which is a *guarantee*, \
         not a measurement."
    );
}
