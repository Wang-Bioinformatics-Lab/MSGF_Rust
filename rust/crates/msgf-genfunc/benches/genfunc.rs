//! Benchmarks the full MS-GF+ significance pipeline per spectrum: preprocess → scored spectrum →
//! generating function → SpecEValue distribution. This is the "MSGF scoring" a database search
//! runs once per spectrum. Uses the real high-res F13 dataset + HighRes model. `cargo bench -p
//! msgf-genfunc`. Needs the gitignored `validation/data/`.

use std::path::PathBuf;
use std::time::Duration;

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use msgf_chem::{mass, scaling};
use msgf_genfunc::graph::{build_reverse_graph, standard_aa_nominal, Aa, PeptideCleavage};
use msgf_genfunc::{compute_into, merge_group, Cleavage, DpScratch};
use msgf_scorer::preprocess::preprocess;
use msgf_scorer::scored_spectrum::ScoredSpectrum;
use rayon::prelude::*;

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
        eprintln!("benches skipped: validation/data absent");
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
                raw: s
                    .peaks
                    .iter()
                    .map(|p| (p.mz as f32, p.intensity as f32))
                    .collect(),
            })
        })
        .collect();
    // 21-aa de-novo set (20 standard + M-oxidation), uniform prob (timing is prob-independent)
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

/// Full per-spectrum SpecEValue: preprocess + scored spectrum + generating function over the
/// `-ti 0,1` isotope mass group (matching the MS-GF+ search: one graph per candidate mass, merged).
fn spec_evalue(
    model: &msgf_scorer::ScoringModel,
    s: &Prepared,
    aa: &[Aa],
    cleave: Cleavage,
    scratch: &mut DpScratch,
) {
    let peaks = preprocess(model, s.charge, s.parent_mass, &s.raw);
    let scored = ScoredSpectrum::from_ranked_peaks(model, s.charge, s.parent_mass, peaks);
    let tables = scored.tables(s.pep_nominal);
    // Build tables + edges once for the largest candidate; each candidate only recomputes node scores.
    let (mut g, _) = build_reverse_graph(
        &scored,
        &tables,
        s.pep_nominal,
        &[s.pep_nominal],
        aa,
        PeptideCleavage::TRYPSIN,
    );
    let mut gfs = Vec::new();
    for p in (s.pep_nominal - 1..=s.pep_nominal).filter(|&p| p > 0) {
        g.recompute_node_scores(&tables, p, &[p]);
        if let Some(gf) = compute_into(scratch, &g, &[p as usize], Some(cleave)) {
            gfs.push(gf);
        }
    }
    black_box(merge_group(&gfs).map(|g| g.spectral_probability(30)));
}

fn benches(c: &mut Criterion) {
    let Some((model, spectra, aa)) = load() else {
        return;
    };
    let cleave = Cleavage {
        credit: 2,
        penalty: -11,
        prob_cleavage_sites: 0.10,
    };
    let mid = {
        let mut idx: Vec<usize> = (0..spectra.len()).collect();
        idx.sort_by_key(|&i| spectra[i].pep_nominal);
        idx[idx.len() / 2]
    };

    c.bench_function("specevalue_one_spectrum", |b| {
        let mut scratch = DpScratch::default();
        b.iter(|| spec_evalue(&model, &spectra[mid], &aa, cleave, &mut scratch))
    });

    let mut g = c.benchmark_group("throughput");
    g.throughput(Throughput::Elements(spectra.len() as u64));
    g.sample_size(10);
    g.measurement_time(Duration::from_secs(15));
    g.bench_function("specevalue_all_spectra", |b| {
        let mut scratch = DpScratch::default();
        b.iter(|| {
            for s in &spectra {
                spec_evalue(&model, s, &aa, cleave, &mut scratch);
            }
        })
    });
    g.finish();

    // multi-core: the per-spectrum work is embarrassingly parallel
    let mut gp = c.benchmark_group("throughput_parallel");
    gp.throughput(Throughput::Elements(spectra.len() as u64));
    gp.sample_size(10);
    gp.measurement_time(Duration::from_secs(15));
    gp.bench_function("specevalue_all_spectra_rayon", |b| {
        b.iter(|| {
            spectra
                .par_iter()
                .for_each_init(DpScratch::default, |scratch, s| {
                    spec_evalue(&model, s, &aa, cleave, scratch)
                });
        })
    });
    gp.finish();

    eprintln!(
        "(benched {} F13 spectra; {} rayon threads)",
        spectra.len(),
        rayon::current_num_threads()
    );
}

criterion_group!(benial, benches);
criterion_main!(benial);
