//! Validates the full RawScore path against MS-GF+'s FastScorer/DBScanScorer:
//! - peptide → cumulative nominal prefix masses vs the oracle's `nominal_prefix_masses`
//! - `raw_score_nodes` vs `node_only_score` (FastScorer, node only)
//! - `ScoredSpectrum::raw_score` (node + edge) vs `full_score` (DBScanScorer)
//!
//! Node/preprocessed peaks come from `f13_scored_spectrum.golden.json`; nominal masses and
//! node/full scores from `f13_rawscore.golden.json` (MS-GF+). Skipped if either golden is absent.

use msgf_chem::peptide;
use msgf_scorer::scored_spectrum::{raw_score_nodes, RankedPeak, ScoredSpectrum};
use serde_json::Value;
use std::collections::HashMap;
use std::path::PathBuf;

fn repo(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

struct Spec {
    charge: i32,
    parent_mass: f32,
    prefix: Vec<f32>,
    suffix: Vec<f32>,
    peaks: Vec<RankedPeak>,
}

#[test]
fn rawscore_matches_msgf() {
    let ss_path = repo("validation/golden/rawscore/f13_scored_spectrum.golden.json");
    let rs_path = repo("validation/golden/rawscore/f13_rawscore.golden.json");
    let param = repo("validation/data/models/HCD_QExactive_Tryp.param");
    if !ss_path.exists() || !rs_path.exists() || !param.exists() {
        eprintln!("skip: rawscore goldens/model absent");
        return;
    }
    let ss: Value = serde_json::from_str(&std::fs::read_to_string(&ss_path).unwrap()).unwrap();
    let rs: Value = serde_json::from_str(&std::fs::read_to_string(&rs_path).unwrap()).unwrap();
    let model = msgf_scorer::read_param_file(&param).unwrap();

    let mut specs: HashMap<i64, Spec> = HashMap::new();
    for s in ss["spectra"].as_array().unwrap() {
        let f = |k: &str| {
            s[k].as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_f64().unwrap() as f32)
                .collect::<Vec<_>>()
        };
        let peaks = s["peaks"]
            .as_array()
            .unwrap()
            .iter()
            .map(|p| RankedPeak {
                mz: p[0].as_f64().unwrap() as f32,
                intensity: p[1].as_f64().unwrap() as f32,
                rank: p[2].as_i64().unwrap() as i32,
            })
            .collect();
        specs.insert(
            s["scan"].as_i64().unwrap(),
            Spec {
                charge: s["charge"].as_i64().unwrap() as i32,
                parent_mass: s["precursor_mass"].as_f64().unwrap() as f32,
                prefix: f("prefix_score"),
                suffix: f("suffix_score"),
                peaks,
            },
        );
    }

    let (mut nominal_ok, mut node_ok, mut full_ok, mut total) = (0, 0, 0, 0);
    let mut full_bad = Vec::new();
    for r in rs["spectra"].as_array().unwrap() {
        let scan = r["scan"].as_i64().unwrap();
        let peptide = r["peptide"].as_str().unwrap();
        let Some(spec) = specs.get(&scan) else {
            continue;
        };
        total += 1;
        let exp_nominal: Vec<i32> = r["nominal_prefix_masses"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_i64().unwrap() as i32)
            .collect();
        let exp_node = r["node_only_score"].as_i64().unwrap() as i32;
        let exp_full = r["full_score"].as_i64().unwrap() as i32;

        let residues = peptide::parse(peptide).expect("parse peptide");
        let nominal = peptide::nominal_prefix_masses(&residues);
        let accurate = peptide::accurate_prefix_masses(&residues);
        let num_mods = peptide::num_mods(&residues) as i32;

        nominal_ok += (nominal == exp_nominal) as i32;
        node_ok += (raw_score_nodes(&spec.prefix, &spec.suffix, &exp_nominal) == exp_node) as i32;

        let sscored = ScoredSpectrum::from_ranked_peaks(
            &model,
            spec.charge,
            spec.parent_mass,
            spec.peaks.clone(),
        );
        let got_full = sscored.raw_score(&nominal, &accurate, num_mods);
        if got_full == exp_full {
            full_ok += 1;
        } else if full_bad.len() < 10 {
            full_bad.push(format!(
                "scan {scan} {peptide}: full got {got_full} vs {exp_full} (node {exp_node})"
            ));
        }
    }

    eprintln!(
        "nominal {nominal_ok}/{total}, node {node_ok}/{total}, FULL(node+edge) {full_ok}/{total}"
    );
    for b in &full_bad {
        eprintln!("  FULL {b}");
    }
    assert_eq!(nominal_ok, total, "peptide→nominal");
    assert_eq!(node_ok, total, "node-only summation");
    assert_eq!(
        full_ok, total,
        "full RawScore (node + edge) must match DBScanScorer"
    );
}
