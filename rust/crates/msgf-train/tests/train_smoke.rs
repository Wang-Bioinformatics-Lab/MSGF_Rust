//! End-to-end trainer test on a **synthetic** corpus — no fetched bytes, no UC artifacts.
//!
//! This is the trainer's half of the clean-room guarantee that `author_a_model_from_scratch`
//! makes for the encoder: a model is produced, written, re-read and scored without any licensed
//! input. It also pins the two properties a trainer must have — a model that round-trips through
//! the format, and scores with the right sign (a top-ranked real ion scores positive, a missing
//! ion negative) — plus reproducibility (counting, so the same corpus gives the same bytes).

use msgf_chem::{mass, peptide};
use msgf_train::corpus::{self, CorpusFilter, CorpusStats};
use msgf_train::{counts, TrainConfig};
use std::path::PathBuf;

/// Deterministic pseudo-random peptides ending in K/R, with b/y peaks and background noise.
fn synthetic_mgf(n: usize) -> String {
    const AAS: &[u8] = b"AGSPVTLNDQKEMHFRYW";
    let mut out = String::new();
    let mut seed: u64 = 0x2026_0724;
    let mut next = move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed
    };
    for i in 0..n {
        let len = 8 + (next() % 8) as usize;
        let mut seq: String = (0..len - 1)
            .map(|_| AAS[(next() % AAS.len() as u64) as usize] as char)
            .collect();
        seq.push(if next() % 2 == 0 { 'K' } else { 'R' });

        let residues = peptide::parse(&seq).unwrap();
        let acc = peptide::accurate_prefix_masses(&residues);
        let pep_mass = acc[acc.len() - 1] + mass::WATER;
        let charge = 2 + (next() % 2) as i32;

        // b and y ions for every cleavage, most intense first; y is the dominant series.
        let mut peaks: Vec<(f64, f64)> = Vec::new();
        for (k, &prefix) in acc[..acc.len() - 1].iter().enumerate() {
            let suffix = acc[acc.len() - 1] - prefix;
            let inten = 1000.0 - k as f64 * 10.0;
            peaks.push((suffix + mass::WATER + mass::PROTON, inten * 2.0)); // y
            if next() % 4 != 0 {
                peaks.push((prefix + mass::PROTON, inten)); // b, sometimes missing
            }
        }
        // background peaks, uniformly spread
        for j in 0..40 {
            peaks.push((150.0 + j as f64 * 23.7 + (next() % 100) as f64 * 0.01, 20.0));
        }
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        out.push_str("BEGIN IONS\n");
        out.push_str(&format!(
            "PEPMASS={}\nCHARGE={}\nSCANS={}\nSEQ={}\n",
            (pep_mass + charge as f64 * mass::PROTON) / charge as f64,
            charge,
            i + 1,
            seq
        ));
        for (mz, it) in peaks {
            out.push_str(&format!("{mz:.5}\t{it:.1}\n"));
        }
        out.push_str("END IONS\n");
    }
    out
}

fn small_config() -> TrainConfig {
    TrainConfig {
        min_psms_per_partition: 60,
        max_partitions_per_charge: 4,
        ..TrainConfig::high_res_hcd_tryptic()
    }
}

fn load(path: &PathBuf) -> Vec<corpus::TrainingPsm> {
    let mut psms = Vec::new();
    let mut stats = CorpusStats::default();
    corpus::read_annotated_mgf(path, &CorpusFilter::default(), &mut psms, &mut stats).unwrap();
    assert!(
        stats.kept > 0 && stats.kept == psms.len(),
        "corpus reader kept nothing: {stats:?}"
    );
    psms
}

#[test]
fn trains_a_scoring_model_from_scratch() {
    let tmp = PathBuf::from(env!("CARGO_TARGET_TMPDIR"));
    let mgf = tmp.join("synthetic_corpus.mgf");
    std::fs::write(&mgf, synthetic_mgf(400)).unwrap();
    let psms = load(&mgf);

    let cfg = small_config();
    let (model, report, _scheme, spectra) = counts::train(&psms, &cfg);
    assert!(spectra > 0, "no spectra counted");
    assert!(!model.partitions.is_empty());
    assert_eq!(
        model.error_dist.len(),
        model.partitions.len(),
        "§7 covers every partition"
    );

    // The trainer must find the series it was fed: y everywhere, b in most partitions.
    let scored: Vec<&msgf_train::counts::PartitionReport> =
        report.iter().filter(|r| !r.ions.is_empty()).collect();
    assert!(!scored.is_empty(), "no partition got a scored ion type");
    assert!(
        scored
            .iter()
            .all(|r| r.ions.iter().any(|(label, ..)| label == "y")),
        "y should be scored in every populated partition"
    );

    // Round-trip through the binary format, then score with the re-read model.
    let bytes = msgf_scorer::write_param(&model);
    let reread = msgf_scorer::read_param(&bytes).expect("re-read our own model");
    assert_eq!(reread, model, "write(read(m)) == m");

    // A y ion at an intensity rank the corpus actually produced must score positive, and its
    // absence must cost, in every partition that scores it. (Which rank a partition sees depends
    // on its m/z segment — the biggest y ions live in the upper segment — so scan the top ranks
    // rather than assuming rank 1.)
    let mut checked = 0;
    for rd in &reread.rank_dist {
        let pi = rd.partition_index;
        let Some(y) = reread.frag_off[pi].iter().find(|f| f.name == "S_1_19") else {
            continue;
        };
        let best = (1..=20)
            .map(|r| reread.node_score(pi, y, r))
            .fold(f32::NEG_INFINITY, f32::max);
        let missing = reread.missing_ion_score(pi, y);
        assert!(
            best > 0.0,
            "a top-ranked y ion must score positive, got {best}"
        );
        assert!(missing < 0.0, "a missing y ion must cost, got {missing}");
        assert!(best > missing);
        checked += 1;
    }
    assert!(checked > 0, "no partition scored the y series");
}

#[test]
fn training_is_reproducible() {
    let tmp = PathBuf::from(env!("CARGO_TARGET_TMPDIR"));
    let mgf = tmp.join("synthetic_corpus_repro.mgf");
    std::fs::write(&mgf, synthetic_mgf(200)).unwrap();
    let psms = load(&mgf);
    let cfg = small_config();

    let a = msgf_scorer::write_param(&counts::train(&psms, &cfg).0);
    let b = msgf_scorer::write_param(&counts::train(&psms, &cfg).0);
    assert_eq!(a, b, "counting is deterministic: same corpus => same bytes");
}

/// Trains on the real MassIVE-KB corpus when it has been fetched, and checks the shape of the
/// result. Skips on a clean checkout — `validation/data/` is not committed.
#[test]
fn trains_on_massive_kb_corpus_if_present() {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join("validation/data/training");
    let Ok(entries) = std::fs::read_dir(&dir) else {
        eprintln!(
            "skip: {} absent (run fetch_reference_data.sh --training)",
            dir.display()
        );
        return;
    };
    let mut files: Vec<PathBuf> = entries
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|x| x == "mgf"))
        .collect();
    files.sort();
    if files.is_empty() {
        eprintln!("skip: no corpus MGF in {}", dir.display());
        return;
    }

    let mut psms = Vec::new();
    let mut stats = CorpusStats::default();
    let filter = CorpusFilter::default();
    corpus::read_annotated_mgf(&files[0], &filter, &mut psms, &mut stats).unwrap();
    assert!(psms.len() > 1000, "expected a substantial corpus");

    let (model, _report, _scheme, _n) = counts::train(&psms, &TrainConfig::high_res_hcd_tryptic());
    let bytes = msgf_scorer::write_param(&model);
    assert_eq!(
        msgf_scorer::read_param(&bytes).unwrap(),
        model,
        "trained model round-trips"
    );
    // The dominant HCD ion types must come out of a real tryptic HCD corpus.
    let names: Vec<&str> = model
        .frag_off
        .iter()
        .flatten()
        .map(|f| f.name.as_str())
        .collect();
    for expect in ["S_1_19", "P_1_1"] {
        assert!(
            names.contains(&expect),
            "{expect} should be scored somewhere"
        );
    }
}
