//! Validates the RawScore *summation machinery* against MS-GF+'s FastScorer:
//! - peptide → cumulative nominal prefix masses (msgf-chem) vs the oracle's `nominal_prefix_masses`
//! - `raw_score_nodes(prefix, suffix, nominal)` vs the oracle's `node_only_score` (FastScorer)
//!
//! Node scores come from `f13_scored_spectrum.golden.json` (the arrays already validated bit-exact);
//! nominal masses + `node_only_score` come from `f13_rawscore.golden.json` (MS-GF+ FastScorer).
//! This isolates the summation + peptide-nominal port from the (separate) end-to-end search
//! question. Skipped if either golden is absent.

use msgf_chem::peptide;
use msgf_scorer::scored_spectrum::raw_score_nodes;
use serde_json::Value;
use std::collections::HashMap;
use std::path::PathBuf;

fn repo(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

#[test]
fn rawscore_summation_matches_fastscorer() {
    let ss_path = repo("validation/golden/rawscore/f13_scored_spectrum.golden.json");
    let rs_path = repo("validation/golden/rawscore/f13_rawscore.golden.json");
    if !ss_path.exists() || !rs_path.exists() {
        eprintln!("skip: rawscore goldens absent");
        return;
    }
    let ss: Value = serde_json::from_str(&std::fs::read_to_string(&ss_path).unwrap()).unwrap();
    let rs: Value = serde_json::from_str(&std::fs::read_to_string(&rs_path).unwrap()).unwrap();

    // prefix/suffix arrays by scan
    let mut arrays: HashMap<i64, (Vec<f32>, Vec<f32>)> = HashMap::new();
    for s in ss["spectra"].as_array().unwrap() {
        let scan = s["scan"].as_i64().unwrap();
        let pre = s["prefix_score"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_f64().unwrap() as f32)
            .collect();
        let suf = s["suffix_score"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_f64().unwrap() as f32)
            .collect();
        arrays.insert(scan, (pre, suf));
    }

    let (mut nominal_ok, mut node_ok, mut total) = (0, 0, 0);
    let mut nominal_bad = Vec::new();
    let mut node_bad = Vec::new();
    for r in rs["spectra"].as_array().unwrap() {
        let scan = r["scan"].as_i64().unwrap();
        let peptide = r["peptide"].as_str().unwrap();
        let exp_nominal: Vec<i32> = r["nominal_prefix_masses"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_i64().unwrap() as i32)
            .collect();
        let exp_node = r["node_only_score"].as_i64().unwrap() as i32;
        let Some((prefix, suffix)) = arrays.get(&scan) else {
            continue;
        };
        total += 1;

        // peptide -> nominal prefix masses (Rust)
        let residues = peptide::parse(peptide).expect("parse peptide");
        let got_nominal = peptide::nominal_prefix_masses(&residues);
        if got_nominal == exp_nominal {
            nominal_ok += 1;
        } else if nominal_bad.len() < 6 {
            nominal_bad.push(format!(
                "scan {scan} {peptide}: got {got_nominal:?} vs {exp_nominal:?}"
            ));
        }

        // node-only RawScore summation over MS-GF+'s own nominal masses
        let got_node = raw_score_nodes(prefix, suffix, &exp_nominal);
        if got_node == exp_node {
            node_ok += 1;
        } else if node_bad.len() < 8 {
            node_bad.push(format!(
                "scan {scan} {peptide}: got {got_node} vs {exp_node}"
            ));
        }
    }

    eprintln!(
        "peptide→nominal: {nominal_ok}/{total} exact; node-only summation: {node_ok}/{total} exact"
    );
    for b in &nominal_bad {
        eprintln!("  NOMINAL {b}");
    }
    for b in &node_bad {
        eprintln!("  NODESUM {b}");
    }
    assert_eq!(
        node_ok, total,
        "node-only summation must match FastScorer on all spectra"
    );
    assert_eq!(
        nominal_ok, total,
        "peptide→nominal must match on all spectra"
    );
}
