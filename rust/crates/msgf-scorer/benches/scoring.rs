//! Benchmarks for the MSGF_Rust scoring pipeline as it stands: model parse, MGF read,
//! per-spectrum preprocessing, and per-spectrum node scoring (prefixScore/suffixScore) — the
//! work a database search does once per spectrum. Uses the real high-res F13 dataset.
//!
//! Run: `cargo bench -p msgf-scorer`. Needs the gitignored `validation/data/` present.

use std::path::PathBuf;
use std::time::Duration;

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use msgf_chem::{mass, scaling};
use msgf_scorer::preprocess::preprocess;
use msgf_scorer::scored_spectrum::ScoredSpectrum;

fn repo(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

/// One spectrum with everything the per-spectrum scoring needs precomputed.
struct Prepared {
    charge: i32,
    parent_mass: f32,
    pep_nominal: i32,
    raw: Vec<(f32, f32)>,
}

fn load() -> Option<(Vec<u8>, Vec<Prepared>)> {
    let param = repo("validation/data/models/HCD_QExactive_Tryp.param");
    let mgf = repo("validation/data/spectra/F13.mgf");
    if !param.exists() || !mgf.exists() {
        eprintln!("benches skipped: validation/data absent (run fetch_reference_data.sh)");
        return None;
    }
    let param_bytes = std::fs::read(&param).unwrap();
    let spectra = msgf_io::read_mgf_file(&mgf)
        .unwrap()
        .into_iter()
        .filter_map(|s| {
            let charge = s.charge?;
            let mz = s.precursor_mz? as f32;
            // de-charged neutral precursor mass; peptide mass = precursor − water
            let parent_mass = mz * charge as f32 - charge as f32 * mass::PROTON as f32;
            let pep_nominal = scaling::nominal_bin(parent_mass - mass::WATER as f32);
            if pep_nominal <= 0 {
                return None;
            }
            Some(Prepared {
                charge,
                parent_mass,
                pep_nominal,
                raw: s
                    .peaks
                    .iter()
                    .map(|p| (p.mz as f32, p.intensity as f32))
                    .collect(),
            })
        })
        .collect();
    Some((param_bytes, spectra))
}

fn benches(c: &mut Criterion) {
    let Some((param_bytes, spectra)) = load() else {
        return;
    };
    let model = msgf_scorer::read_param(&param_bytes).unwrap();
    // a representative spectrum: median peak count
    let mid = {
        let mut idx: Vec<usize> = (0..spectra.len()).collect();
        idx.sort_by_key(|&i| spectra[i].raw.len());
        idx[idx.len() / 2]
    };
    let s0 = &spectra[mid];

    c.bench_function("read_param_model", |b| {
        b.iter(|| msgf_scorer::read_param(black_box(&param_bytes)).unwrap())
    });

    c.bench_function("preprocess_one", |b| {
        b.iter(|| preprocess(&model, s0.charge, s0.parent_mass, black_box(&s0.raw)))
    });

    c.bench_function("score_one_spectrum", |b| {
        b.iter(|| {
            let peaks = preprocess(&model, s0.charge, s0.parent_mass, &s0.raw);
            let ss = ScoredSpectrum::from_ranked_peaks(&model, s0.charge, s0.parent_mass, peaks);
            black_box(ss.prefix_suffix_scores(s0.pep_nominal))
        })
    });

    let mut g = c.benchmark_group("throughput");
    g.throughput(Throughput::Elements(spectra.len() as u64));
    g.sample_size(20);
    g.measurement_time(Duration::from_secs(8));
    g.bench_function("preprocess_and_score_all", |b| {
        b.iter(|| {
            for s in &spectra {
                let peaks = preprocess(&model, s.charge, s.parent_mass, &s.raw);
                let ss = ScoredSpectrum::from_ranked_peaks(&model, s.charge, s.parent_mass, peaks);
                black_box(ss.prefix_suffix_scores(s.pep_nominal));
            }
        })
    });
    g.finish();
}

criterion_group!(benial, benches);
criterion_main!(benial);
