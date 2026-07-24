//! Validates the `.param` encoder (`msgf_scorer::write_param`) against the reader.
//!
//! Two independent guarantees:
//!   1. `real_models_round_trip` — for every UC `.param` present, `read(write(read(f))) == read(f)`,
//!      and it reports whether the re-encoding is byte-for-byte identical. The `.param` data is
//!      gitignored; the test skips gracefully when it is absent (fresh clone / CI).
//!   2. `author_a_model_from_scratch` — builds a small [`ScoringModel`] with NO fetched data, writes
//!      it, reads it back, and scores a node. This is the clean-room authoring path a future
//!      `msgf-train` will emit into, so it must pass on any checkout, data or not.

use msgf_chem::Tolerance;
use msgf_scorer::{
    read_param, write_param, ErrorDist, FragOff, Partition, PrecursorOff, RankDist, ScoringModel,
};
use std::path::PathBuf;

fn repo(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

#[test]
fn real_models_round_trip() {
    let mdir = repo("validation/data/models");
    let Ok(dir) = std::fs::read_dir(&mdir) else {
        eprintln!("skip: {} absent (fetch_reference_data.sh)", mdir.display());
        return;
    };

    let mut checked = 0;
    for entry in dir.flatten() {
        let name = entry.file_name().to_string_lossy().into_owned();
        if !name.ends_with(".param") {
            continue;
        }
        let bytes = std::fs::read(entry.path()).unwrap();
        let m = read_param(&bytes).unwrap_or_else(|e| panic!("read {name}: {e}"));

        let out = write_param(&m);
        let m2 = read_param(&out).unwrap_or_else(|e| panic!("re-read {name}: {e}"));

        // The load-bearing guarantee: the encoder produces a stream the reader decodes identically.
        assert_eq!(m2, m, "struct round-trip differs for {name}");

        // Byte-identity is expected too, except for files with exact-zero ion-existence entries
        // (the reader floors those to 0.001, so the re-encode legitimately differs there).
        if out == bytes {
            eprintln!("ok {name}: byte-for-byte identical ({} bytes)", bytes.len());
        } else {
            let diffs = byte_diffs(&bytes, &out);
            eprintln!(
                "ok {name}: struct-exact, {diffs} byte diffs (ion-existence 0.0->0.001 floor)"
            );
        }
        checked += 1;
    }

    if checked == 0 {
        eprintln!(
            "skip: no .param files in {} (fetch_reference_data.sh)",
            mdir.display()
        );
    }
}

fn byte_diffs(a: &[u8], b: &[u8]) -> usize {
    let n = a.len().min(b.len());
    let mut d = a.len().abs_diff(b.len());
    for i in 0..n {
        if a[i] != b[i] {
            d += 1;
        }
    }
    d
}

#[test]
fn author_a_model_from_scratch() {
    // A minimal, hand-authored high-res-style model: one charge-2 partition scoring b (P_1_1) and
    // y (S_1_19) ions against a 4-column rank table. No `.param` bytes, no fetched data — this is
    // exactly what a trainer would assemble from confident PSMs, in miniature.
    //
    // `FragOff::name` must equal what the reader re-derives ("{P|S}_{charge}_{round(offset)}"),
    // and each rank row must have `max_rank + 1` columns, or the round-trip assertion below fails.
    let ncols = 4; // max_rank (3) + 1
    let model = ScoringModel {
        version: 1,
        activation: "HCD".into(),
        instrument: "HighRes".into(),
        enzyme: Some("Tryp".into()),
        protocol: None, // -> "Automatic"
        mme: Tolerance::ppm(10.0),
        apply_deconvolution: true,
        deconvolution_error_tolerance: 0.01,
        charge_histogram: vec![(2, 100)],
        num_segments: 1,
        partitions: vec![Partition {
            charge: 2,
            parent_mass: 0.0,
            seg: 0,
        }],
        precursor_off: vec![PrecursorOff {
            charge: 2,
            reduced_charge: 1,
            offset: 0.0,
            tol_ppm: true,
            tol_val: 10.0,
            frequency: 1.0,
        }],
        frag_off: vec![vec![
            FragOff {
                is_prefix: true,
                charge: 1,
                offset: 1.007_825, // proton -> round = 1 -> "P_1_1"
                frequency: 0.9,
                name: "P_1_1".into(),
            },
            FragOff {
                is_prefix: false,
                charge: 1,
                offset: 19.018_39, // y-ion offset -> round = 19 -> "S_1_19"
                frequency: 0.8,
                name: "S_1_19".into(),
            },
        ]],
        max_rank: 3,
        rank_dist: vec![RankDist {
            partition_index: 0,
            // ion rows in frag_off order, then noise last (the order the reader expects)
            ions: vec![
                ("P_1_1".into(), vec![0.50, 0.30, 0.15, 0.05]),
                ("S_1_19".into(), vec![0.40, 0.30, 0.20, 0.10]),
                ("noise".into(), vec![0.10, 0.20, 0.30, 0.40]),
            ],
        }],
        error_scaling_factor: 0,
        error_dist: Vec::<ErrorDist>::new(),
    };
    // sanity: rank rows are the right width
    for (_n, row) in &model.rank_dist[0].ions {
        assert_eq!(row.len(), ncols, "rank row must be max_rank+1 wide");
    }

    // write -> read must reproduce the model exactly
    let bytes = write_param(&model);
    let back = read_param(&bytes).expect("our own model must decode");
    assert_eq!(back, model, "hand-authored model failed to round-trip");

    // and it must actually score: node_score = ln(ionFreq[rank] / (noiseFreq[rank] * min(charge,segs)))
    // for b at rank 1: ln(0.50 / (0.10 * min(1,1))) = ln(5)
    let b_ion = &back.frag_off[0][0];
    let s = back.node_score(0, b_ion, 1);
    let expected = (0.50_f64 / (0.10 * 1.0)).ln() as f32;
    assert!(
        (s - expected).abs() < 1e-6,
        "node_score {s} != expected {expected}"
    );

    // missing-ion score uses the last (max_rank) column: ln(0.05 / (0.40 * 1))
    let miss = back.missing_ion_score(0, b_ion);
    let expected_miss = (0.05_f64 / (0.40 * 1.0)).ln() as f32;
    assert!(
        (miss - expected_miss).abs() < 1e-6,
        "missing_ion_score {miss} != expected {expected_miss}"
    );

    eprintln!(
        "ok: authored a {}-partition model from scratch, wrote {} bytes, decodes + scores",
        model.partitions.len(),
        bytes.len()
    );
}
