//! Benchmarks for the scoring half: model decode, and per-spectrum preparation (peak order,
//! precursor suppression, ranks, isotope reduction, node/edge tables) — the work a search does
//! once per (spectrum, charge). Uses the F13 spectra and the bundled model.
//!
//! Run: `cargo bench -p msgf-scorer`. Needs `validation/data/spectra/F13.mgf`
//! (`validation/fetch_reference_data.sh`).

use std::path::PathBuf;
use std::time::Duration;

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use msgf_scorer::{bundled, prepare, Ms2, ScoreModel};

fn repo(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

/// (charge, precursor m/z, peaks) per F13 spectrum.
type Raw = (i32, f64, Vec<(f64, f64)>);

fn load() -> Option<Vec<Raw>> {
    let mgf = repo("validation/data/spectra/F13.mgf");
    if !mgf.exists() {
        eprintln!("benches skipped: validation/data absent (run fetch_reference_data.sh)");
        return None;
    }
    Some(
        msgf_io::read_mgf_file(&mgf)
            .unwrap()
            .into_iter()
            .filter_map(|s| {
                let peaks = s.peaks.iter().map(|p| (p.mz, p.intensity)).collect();
                Some((s.charge?, s.precursor_mz?, peaks))
            })
            .collect(),
    )
}

fn benches(c: &mut Criterion) {
    let Some(spectra) = load() else {
        return;
    };
    let model = bundled::score_model().unwrap();
    fn ms2(s: &Raw) -> Ms2<'_> {
        Ms2 {
            peaks: &s.2,
            precursor_mz: s.1,
            charge: s.0,
        }
    }
    let mid = {
        let mut idx: Vec<usize> = (0..spectra.len()).collect();
        idx.sort_by_key(|&i| spectra[i].2.len());
        idx[idx.len() / 2]
    };

    c.bench_function("decode_bundled_model", |b| {
        b.iter(|| ScoreModel::from_bytes(black_box(bundled::PARAM)).unwrap())
    });
    c.bench_function("prepare_one", |b| {
        b.iter(|| prepare(&model, black_box(&ms2(&spectra[mid]))))
    });
    let mut g = c.benchmark_group("prepare_all");
    g.throughput(Throughput::Elements(spectra.len() as u64));
    g.measurement_time(Duration::from_secs(10));
    g.bench_function("prepare_all", |b| {
        b.iter(|| {
            for s in &spectra {
                black_box(prepare(&model, &ms2(s)));
            }
        })
    });
    g.finish();
}

criterion_group!(scoring, benches);
criterion_main!(scoring);
