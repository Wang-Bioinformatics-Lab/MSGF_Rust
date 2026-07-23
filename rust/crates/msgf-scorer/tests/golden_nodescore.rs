//! Validates the scoring primitives (`node_score` / `missing_ion_score`) against
//! `validation/golden/models/node_scores.golden.json`, which holds MS-GF+'s own
//! `getNodeScore` / `getMissingIonScore` outputs for every partition, ion, and a range of
//! ranks (plus the MISSING bin), across all four high-res models. `.param` data is gitignored;
//! a missing model is skipped.

use msgf_scorer::FragOff;
use serde_json::Value;
use std::collections::HashMap;
use std::path::PathBuf;

fn repo(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

#[test]
fn node_scores_match() {
    let g: Value = serde_json::from_str(
        &std::fs::read_to_string(repo("validation/golden/models/node_scores.golden.json")).unwrap(),
    )
    .unwrap();

    let mut models_checked = 0;
    for (file, data) in g["models"].as_object().unwrap() {
        let param = repo("validation/data/models").join(file);
        if !param.exists() {
            eprintln!("skip {file}: data absent");
            continue;
        }
        let m = msgf_scorer::read_param_file(&param).unwrap();
        // per-partition ion lookup by name
        let by_pi: Vec<HashMap<&str, &FragOff>> = m
            .frag_off
            .iter()
            .map(|b| b.iter().map(|f| (f.name.as_str(), f)).collect())
            .collect();

        let rows = data["rows"].as_array().unwrap();
        let mut sum = 0.0f64;
        let mut worst = 0.0f64;
        for row in rows {
            let pi = row[0].as_u64().unwrap() as usize;
            let ion_name = row[1].as_str().unwrap();
            let ion_charge = row[2].as_i64().unwrap();
            let rank = row[3].as_str().unwrap();
            let expect = row[4].as_f64().unwrap();

            let fo = by_pi[pi]
                .get(ion_name)
                .unwrap_or_else(|| panic!("{file}: partition {pi} has no ion {ion_name}"));
            assert_eq!(fo.charge as i64, ion_charge, "{file} {ion_name} charge");

            let got = if rank == "MISSING" {
                m.missing_ion_score(pi, fo)
            } else {
                m.node_score(pi, fo, rank.parse().unwrap())
            } as f64;

            let d = (got - expect).abs();
            worst = worst.max(d);
            assert!(
                d <= 1e-4,
                "{file} p{pi} {ion_name} rank {rank}: {got} vs {expect} (Δ={d:.2e})"
            );
            sum += got;
        }

        assert_eq!(
            rows.len() as u64,
            data["count"].as_u64().unwrap(),
            "{file} row count"
        );
        assert!(
            (sum - data["score_sum"].as_f64().unwrap()).abs() <= 0.1,
            "{file} score_sum"
        );
        eprintln!(
            "ok {file}: {} node scores match (worst Δ={worst:.2e})",
            rows.len()
        );
        models_checked += 1;
    }
    assert!(
        models_checked > 0,
        "no .param data present to validate node scores"
    );
}
