//! Validates the Rust scoring stack against MS-GF+ on the **HighRes** model (the one the F13
//! search actually used), across three layers: preprocessing (raw F13 peaks → preprocessed
//! ranked peaks), node scoring (prefix/suffix arrays), and edge scoring (full RawScore). Any
//! discrepancy here is what costs the last decimal of the generating-function p-value.
//! Skipped if the HighRes goldens/model/data are absent.

use msgf_chem::peptide;
use msgf_scorer::preprocess::preprocess;
use msgf_scorer::scored_spectrum::{RankedPeak, ScoredSpectrum};
use serde_json::Value;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

fn repo(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

#[test]
fn highres_scoring_matches_msgf() {
    let ss_path = repo("validation/golden/rawscore/f13_scored_spectrum_highres.golden.json");
    let rs_path = repo("validation/golden/rawscore/f13_rawscore_highres.golden.json");
    let mgf = repo("validation/data/spectra/F13.mgf");
    let param = repo("validation/data/models/HCD_HighRes_Tryp.param");
    if !ss_path.exists() || !rs_path.exists() || !mgf.exists() || !param.exists() {
        eprintln!("skip: HighRes goldens/model/data absent");
        return;
    }
    let ss: Value = serde_json::from_str(&std::fs::read_to_string(&ss_path).unwrap()).unwrap();
    let rs: Value = serde_json::from_str(&std::fs::read_to_string(&rs_path).unwrap()).unwrap();
    let model = msgf_scorer::read_param_file(&param).unwrap();

    // raw peaks by scan
    let mut raw: HashMap<i64, Vec<(f32, f32)>> = HashMap::new();
    for s in msgf_io::MgfReader::new(BufReader::new(File::open(&mgf).unwrap())) {
        let s = s.unwrap();
        if let Some(scan) = s.scan.as_deref().and_then(|x| x.parse::<i64>().ok()) {
            raw.insert(
                scan,
                s.peaks
                    .iter()
                    .map(|p| (p.mz as f32, p.intensity as f32))
                    .collect(),
            );
        }
    }
    // HighRes rawscore golden by scan
    let mut rsmap: HashMap<i64, &Value> = HashMap::new();
    for r in rs["spectra"].as_array().unwrap() {
        rsmap.insert(r["scan"].as_i64().unwrap(), r);
    }

    let (mut prep_ok, mut node_ok, mut full_ok, mut total) = (0, 0, 0, 0);
    let mut worst_node = 0.0f64;
    for s in ss["spectra"].as_array().unwrap() {
        let scan = s["scan"].as_i64().unwrap();
        let charge = s["charge"].as_i64().unwrap() as i32;
        let parent_mass = s["precursor_mass"].as_f64().unwrap() as f32;
        let pep_nominal = s["peptide_mass_nominal"].as_i64().unwrap() as i32;
        total += 1;

        let gpeaks: Vec<RankedPeak> = s["peaks"]
            .as_array()
            .unwrap()
            .iter()
            .map(|p| RankedPeak {
                mz: p[0].as_f64().unwrap() as f32,
                intensity: p[1].as_f64().unwrap() as f32,
                rank: p[2].as_i64().unwrap() as i32,
            })
            .collect();

        // (1) preprocessing: my preprocess(raw) vs MS-GF+'s HighRes peaks — n_peaks + ranks
        let my_peaks = preprocess(&model, charge, parent_mass, &raw[&scan]);
        let prep_match = my_peaks.len() == gpeaks.len()
            && my_peaks
                .iter()
                .zip(&gpeaks)
                .all(|(a, b)| a.rank == b.rank && (a.mz - b.mz).abs() < 1e-2);
        prep_ok += prep_match as i32;

        // (2) node scoring on MS-GF+'s HighRes peaks: prefix/suffix arrays
        let scored = ScoredSpectrum::from_ranked_peaks(&model, charge, parent_mass, gpeaks.clone());
        let (prefix, suffix) = scored.prefix_suffix_scores(pep_nominal);
        let mut node_match = true;
        for (name, got, exp) in [
            ("prefix", &prefix, &s["prefix_score"]),
            ("suffix", &suffix, &s["suffix_score"]),
        ] {
            let exp = exp.as_array().unwrap();
            for nm in 1..got.len().min(exp.len()) {
                let d = (got[nm] as f64 - exp[nm].as_f64().unwrap()).abs();
                worst_node = worst_node.max(d);
                if d > 2e-3 {
                    node_match = false;
                }
            }
            let _ = name;
        }
        node_ok += node_match as i32;

        // (3) edge/raw scoring vs full_score (node+edge)
        if let Some(r) = rsmap.get(&scan) {
            let residues = peptide::parse(r["peptide"].as_str().unwrap()).unwrap();
            let nominal = peptide::nominal_prefix_masses(&residues);
            let accurate = peptide::accurate_prefix_masses(&residues);
            let num_mods = peptide::num_mods(&residues) as i32;
            let got = scored.raw_score(&nominal, &accurate, num_mods);
            full_ok += (got == r["full_score"].as_i64().unwrap() as i32) as i32;
        }
    }
    eprintln!("HighRes: preprocess {prep_ok}/{total}, node scores {node_ok}/{total} (worst Δ {worst_node:.2e}), full RawScore {full_ok}/{total}");
    assert_eq!(prep_ok, total, "HighRes preprocessing");
    assert_eq!(node_ok, total, "HighRes node scoring");
    assert_eq!(full_ok, total, "HighRes full RawScore");
}
