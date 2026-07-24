//! `PLAN2.md` TD-2 Gate 2: reproduce `edu.ucsd.msjava.fdr.TargetDecoyAnalysis` entry-for-entry.
//!
//! The oracle is `validation/golden/fdr/fdrmap_cases.golden.json`, dumped by
//! `validation/reference/java/DumpFdrMap.java` (regenerate with
//! `validation/reference/make_fdr_golden.sh`, which needs only the jar and a JVM — no spectra,
//! models or database). Being MS-GF+-derived it is **not committed**, so this test skips until you
//! generate it.
//!
//! Gate 1 (`golden_fdr.rs`) pins the F13 search columns but only ever sees two distinct q-values
//! (`PLAN2.md` §4). These 14 synthetic cases are built to separate the rules it cannot see: how a
//! run of equal decoy scores is charged, what happens when no target beats a decoy, and whether a
//! score sitting exactly on a threshold takes that threshold's q-value or the next one's — which
//! is why every threshold's immediate float neighbours are probed.

use msgf_fdr::QValueMap;
use serde_json::Value;
use std::path::PathBuf;

fn repo(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

/// Java writes every float as `Float.toString`, which round-trips exactly through Rust's parser.
fn parse_f32(v: &Value) -> f32 {
    v.as_str()
        .expect("floats are encoded as strings")
        .parse()
        .expect("Float.toString is parseable")
}

fn show(v: f32) -> String {
    format!("{v:e} (bits {:#010x})", v.to_bits())
}

#[test]
fn fdr_map_and_lookups_match_msgfplus() {
    let path = repo("validation/golden/fdr/fdrmap_cases.golden.json");
    if !path.exists() {
        eprintln!(
            "skip: {} absent (validation/reference/make_fdr_golden.sh)",
            path.display()
        );
        return;
    }
    let g: Value = serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
    let cases = g["cases"].as_array().expect("cases array");
    assert!(!cases.is_empty(), "golden has no cases");

    let mut checked_entries = 0usize;
    let mut checked_lookups = 0usize;
    let mut undefined_lookups = 0usize;
    let mut failures: Vec<String> = Vec::new();

    for case in cases {
        let name = case["name"].as_str().unwrap();
        let targets: Vec<f32> = case["targets"]
            .as_array()
            .unwrap()
            .iter()
            .map(parse_f32)
            .collect();
        let decoys: Vec<f32> = case["decoys"]
            .as_array()
            .unwrap()
            .iter()
            .map(parse_f32)
            .collect();
        let map = QValueMap::build(&targets, &decoys, 1.0);

        // 1. The threshold -> q-value map, entry for entry, including both sentinels.
        let want = case["map"].as_array().unwrap();
        let got = map.pairs();
        if got.len() != want.len() {
            failures.push(format!(
                "{name}: map has {} entries, MS-GF+ has {} ({:?} vs {:?})",
                got.len(),
                want.len(),
                got.iter().map(|p| p.0).collect::<Vec<_>>(),
                want.iter()
                    .map(|e| parse_f32(&e["key"]))
                    .collect::<Vec<_>>(),
            ));
            continue;
        }
        for (i, entry) in want.iter().enumerate() {
            let (want_key, want_q) = (parse_f32(&entry["key"]), parse_f32(&entry["q"]));
            let (got_key, got_q) = got[i];
            checked_entries += 1;
            if got_key.to_bits() != want_key.to_bits() {
                failures.push(format!(
                    "{name}: threshold {i} is {}, MS-GF+ has {}",
                    show(got_key),
                    show(want_key)
                ));
            } else if got_q.to_bits() != want_q.to_bits() {
                failures.push(format!(
                    "{name}: q at threshold {} is {}, MS-GF+ has {}",
                    show(got_key),
                    show(got_q),
                    show(want_q)
                ));
            }
        }

        // 2. Every probed lookup. `q` is null where Java's getPSMQValue has no answer at all
        // (it dereferences a null map entry); nothing to compare there.
        for lookup in case["lookups"].as_array().unwrap() {
            let score = parse_f32(&lookup["score"]);
            if lookup["q"].is_null() {
                undefined_lookups += 1;
                continue;
            }
            let want_q = parse_f32(&lookup["q"]);
            let got_q = map.q_value(score);
            checked_lookups += 1;
            if got_q.to_bits() != want_q.to_bits() {
                failures.push(format!(
                    "{name}: q_value({}) = {}, MS-GF+ has {}",
                    show(score),
                    show(got_q),
                    show(want_q)
                ));
            }
        }
    }

    assert!(
        failures.is_empty(),
        "{} divergences from MS-GF+ (first 10): {:#?}",
        failures.len(),
        failures.iter().take(10).collect::<Vec<_>>()
    );
    eprintln!(
        "FDR map: {} cases, {checked_entries} thresholds and {checked_lookups} lookups exact vs MS-GF+ \
         ({undefined_lookups} undefined in Java)",
        cases.len()
    );
}
