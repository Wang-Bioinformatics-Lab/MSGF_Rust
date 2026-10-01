//! Benchmarks the significance stage per spectrum: prepare → null score distribution (DeNovoScore
//! and the full SpecEValue tail, unpruned) over the isotope range (0, 1). Uses the F13 spectra and
//! the bundled model. `cargo bench -p msgf-genfunc`. Needs `validation/data/spectra/F13.mgf`.

use std::path::PathBuf;
use std::time::Duration;

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use msgf_genfunc::{score_distribution, NullModel};
use msgf_scorer::{bundled, prepare, Ms2};

fn repo(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

fn benches(c: &mut Criterion) {
    let mgf = repo("validation/data/spectra/F13.mgf");
    if !mgf.exists() {
        eprintln!("benches skipped: validation/data absent");
        return;
    }
    let spectra: Vec<(i32, f64, Vec<(f64, f64)>)> = msgf_io::read_mgf_file(&mgf)
        .unwrap()
        .into_iter()
        .filter_map(|s| {
            let peaks = s.peaks.iter().map(|p| (p.mz, p.intensity)).collect();
            Some((s.charge?, s.precursor_mz?, peaks))
        })
        .collect();
    let model = bundled::score_model().unwrap();
    let null = NullModel::uniform();
    let mut g = c.benchmark_group("significance");
    g.throughput(Throughput::Elements(spectra.len() as u64));
    g.measurement_time(Duration::from_secs(20));
    g.sample_size(10);
    g.bench_function("prepare_and_distribution_all", |b| {
        b.iter(|| {
            for (z, mz, peaks) in &spectra {
                let ms2 = Ms2 {
                    peaks,
                    precursor_mz: *mz,
                    charge: *z,
                };
                if let Some(prep) = prepare(&model, &ms2) {
                    black_box(score_distribution(&prep, &null, None));
                }
            }
        })
    });
    g.finish();
}

criterion_group!(genfunc, benches);
criterion_main!(genfunc);
