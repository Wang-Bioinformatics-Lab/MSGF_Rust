//! Measures **tail-threshold pruning** of the generating-function DP on real data, and checks it
//! against the unpruned DP bit-for-bit.
//!
//! SpecEValue only needs the upper tail `P(score >= T)`. With `max_rem[m]` = the best score any
//! path can still earn from node `m` to the sink, a cell `(m, s)` with `s + max_rem[m] < T` can
//! only ever land below `T`, so it is irrelevant to the tail. `compute_tail_into` drops those.
//!
//! Thresholds come from the frozen MS-GF+ F13 output (the observed RawScore of its top hit) — the
//! same threshold a search knows once it has scored the candidates for a spectrum.
//!
//! Run: cargo run -p msgf-genfunc --example tailprune --release
//! Needs the gitignored validation/data/ + the F13 golden.

use std::collections::HashMap;
use std::path::PathBuf;
use std::time::{Duration, Instant};

use msgf_chem::{mass, scaling};
use msgf_genfunc::graph::{build_reverse_graph, standard_aa_nominal, Aa, PeptideCleavage};
use msgf_genfunc::{compute_into, compute_tail_into, merge_group, Cleavage, DpScratch, GenFunc};
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

struct Prepared {
    charge: i32,
    parent_mass: f32,
    pep_nominal: i32,
    raw: Vec<(f32, f32)>,
    /// MS-GF+'s own top-hit RawScore for this scan, when the golden has one.
    threshold: Option<i32>,
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

/// Best observed RawScore (`MSGFScore`) per scan, from the frozen MS-GF+ F13 output.
fn golden_rawscore() -> HashMap<i32, i32> {
    let mut out = HashMap::new();
    let Ok(text) = std::fs::read_to_string(repo("validation/golden/iprg2013_F13.tsv")) else {
        eprintln!("note: F13 golden absent — falling back to DeNovoScore-derived thresholds");
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

/// The per-spectrum work that does not depend on the threshold, kept out of the timed section.
struct SpectrumWork {
    graph: msgf_genfunc::Graph,
    tables: msgf_scorer::scored_spectrum::SpectrumTables,
    candidates: Vec<i32>,
}

fn prepare(model: &msgf_scorer::ScoringModel, s: &Prepared, aa: &[Aa]) -> SpectrumWork {
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
    SpectrumWork {
        graph,
        tables,
        candidates: (s.pep_nominal - 1..=s.pep_nominal)
            .filter(|&p| p > 0)
            .collect(),
    }
}

/// Run the isotope-error group, either fully or tail-pruned at `threshold`.
fn run(w: &mut SpectrumWork, sc: &mut DpScratch, threshold: Option<i32>) -> Option<GenFunc> {
    let mut gfs = Vec::with_capacity(w.candidates.len());
    for &p in &w.candidates {
        w.graph.recompute_node_scores(&w.tables, p, &[p]);
        let gf = match threshold {
            Some(t) => compute_tail_into(sc, &w.graph, &[p as usize], Some(CLEAVE), t),
            None => compute_into(sc, &w.graph, &[p as usize], Some(CLEAVE)),
        };
        if let Some(gf) = gf {
            gfs.push(gf);
        }
    }
    merge_group(&gfs)
}

fn main() {
    let Some((model, spectra, aa)) = load() else {
        return;
    };
    let mut sc = DpScratch::default();
    let mut prepared: Vec<(SpectrumWork, i32)> = Vec::new();

    // Threshold per spectrum: MS-GF+'s observed top-hit RawScore where the golden has one.
    let mut n_golden = 0;
    for s in &spectra {
        let w = prepare(&model, s, &aa);
        let t = match s.threshold {
            Some(t) => {
                n_golden += 1;
                t
            }
            None => continue,
        };
        prepared.push((w, t));
    }
    println!(
        "{} F13 spectra with a golden MS-GF+ RawScore (of {})",
        n_golden,
        spectra.len()
    );

    // --- correctness: the pruned tail must be the identical f64 -------------------------------
    let mut checked = 0usize;
    let mut cells_full = 0u64;
    let mut cells_tail = 0u64;
    for (w, t) in prepared.iter_mut() {
        let full = run(w, &mut sc, None);
        let tail = run(w, &mut sc, Some(*t));
        match (full, tail) {
            (Some(f), Some(p)) => {
                assert_eq!(
                    f.max_score(),
                    p.max_score(),
                    "DeNovoScore changed by pruning"
                );
                let (a, b) = (f.spectral_probability(*t), p.spectral_probability(*t));
                assert_eq!(
                    a.to_bits(),
                    b.to_bits(),
                    "SpecEValue differs at threshold {t}: full {a:e} vs pruned {b:e}"
                );
                checked += 1;
            }
            (None, None) => {}
            (f, p) => panic!("reachability disagreed: {} vs {}", f.is_some(), p.is_some()),
        }
    }
    println!("bit-identical SpecEValue on {checked} spectra (f64 to_bits)");

    // --- DP cell counts (the convolution work) -------------------------------------------------
    for (w, t) in prepared.iter_mut() {
        for &p in &w.candidates.clone() {
            w.graph.recompute_node_scores(&w.tables, p, &[p]);
            compute_into(&mut sc, &w.graph, &[p as usize], Some(CLEAVE));
            cells_full += sc.arena_len() as u64;
            w.graph.recompute_node_scores(&w.tables, p, &[p]);
            compute_tail_into(&mut sc, &w.graph, &[p as usize], Some(CLEAVE), *t);
            cells_tail += sc.arena_len() as u64;
        }
    }

    // --- timing ---------------------------------------------------------------------------------
    let time = |sc: &mut DpScratch, prep: &mut Vec<(SpectrumWork, i32)>, prune: bool| -> Duration {
        let mut best = Duration::MAX;
        for _ in 0..3 {
            let t0 = Instant::now();
            for (w, t) in prep.iter_mut() {
                std::hint::black_box(run(w, sc, if prune { Some(*t) } else { None }));
            }
            best = best.min(t0.elapsed());
        }
        best
    };
    let t_full = time(&mut sc, &mut prepared, false);
    let t_tail = time(&mut sc, &mut prepared, true);

    println!("\n=== DP only (graph + tables prebuilt, excluded from timing) ===");
    println!(
        "cells:  full {:>12}   pruned {:>12}   {:.2}x fewer",
        cells_full,
        cells_tail,
        cells_full as f64 / cells_tail.max(1) as f64
    );
    println!(
        "time:   full {:>9.1}ms   pruned {:>9.1}ms   {:.2}x faster   (F13's own MS-GF+ RawScores)",
        t_full.as_secs_f64() * 1e3,
        t_tail.as_secs_f64() * 1e3,
        t_full.as_secs_f64() / t_tail.as_secs_f64()
    );

    // --- how the payoff scales with match quality ------------------------------------------------
    // F13 identifies essentially nothing, so its DeNovoScore−RawScore gap (median 69) is a
    // worst case. Re-time at thresholds a fixed distance below each spectrum's own DeNovoScore to
    // show what a corpus with real identifications would see.
    let denovo: Vec<i32> = prepared
        .iter_mut()
        .map(|(w, _)| run(w, &mut sc, None).map(|g| g.max_score()).unwrap_or(0))
        .collect();
    println!("\n=== pruned-DP time vs. match quality (threshold = DeNovoScore − gap) ===");
    println!("  {:>5}   {:>10}   {:>9}", "gap", "time", "speedup");
    for gap in [5i32, 10, 20, 30, 40, 60, 80] {
        let mut best = Duration::MAX;
        for _ in 0..3 {
            let t0 = Instant::now();
            for ((w, _), dn) in prepared.iter_mut().zip(&denovo) {
                std::hint::black_box(run(w, &mut sc, Some(dn - gap)));
            }
            best = best.min(t0.elapsed());
        }
        println!(
            "  {:>5}   {:>8.1}ms   {:>8.2}x",
            gap,
            best.as_secs_f64() * 1e3,
            t_full.as_secs_f64() / best.as_secs_f64()
        );
    }
}
