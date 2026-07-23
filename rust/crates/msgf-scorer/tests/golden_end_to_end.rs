//! End-to-end seam test: raw MGF peaks → `preprocess` → `ScoredSpectrum` → `prefixScore`/
//! `suffixScore`, checked against MS-GF+'s arrays. The two stage tests validate each half against
//! the same golden; this confirms they compose — i.e. the scored spectrum works on `preprocess`'s
//! actual output, not just MS-GF+'s. Skipped if the gitignored data/golden are absent.

use msgf_io::MgfReader;
use msgf_scorer::preprocess::preprocess;
use msgf_scorer::scored_spectrum::ScoredSpectrum;
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
fn raw_mgf_to_node_scores_matches_msgf() {
    let gpath = repo("validation/golden/rawscore/f13_scored_spectrum.golden.json");
    let mgf = repo("validation/data/spectra/F13.mgf");
    let param = repo("validation/data/models/HCD_QExactive_Tryp.param");
    if !gpath.exists() || !mgf.exists() || !param.exists() {
        eprintln!("skip: data/golden absent");
        return;
    }
    let g: Value = serde_json::from_str(&std::fs::read_to_string(&gpath).unwrap()).unwrap();
    let model = msgf_scorer::read_param_file(&param).unwrap();

    // index raw peaks by scan number
    let mut raw: HashMap<i64, Vec<(f32, f32)>> = HashMap::new();
    for sp in MgfReader::new(BufReader::new(File::open(&mgf).unwrap())) {
        let sp = sp.unwrap();
        if let Some(scan) = sp.scan.as_deref().and_then(|s| s.parse::<i64>().ok()) {
            raw.insert(
                scan,
                sp.peaks
                    .iter()
                    .map(|p| (p.mz as f32, p.intensity as f32))
                    .collect(),
            );
        }
    }

    let mut worst = 0.0f64;
    let mut worst_ctx = String::new();
    let mut n = 0;
    for sp in g["spectra"].as_array().unwrap() {
        let scan = sp["scan"].as_i64().unwrap();
        let charge = sp["charge"].as_i64().unwrap() as i32;
        let parent_mass = sp["precursor_mass"].as_f64().unwrap() as f32;
        let pep_nominal = sp["peptide_mass_nominal"].as_i64().unwrap() as i32;
        let Some(raw_peaks) = raw.get(&scan) else {
            panic!("scan {scan} not found in F13.mgf");
        };

        let peaks = preprocess(&model, charge, parent_mass, raw_peaks);
        let ss = ScoredSpectrum::from_ranked_peaks(&model, charge, parent_mass, peaks);
        let (prefix, suffix) = ss.prefix_suffix_scores(pep_nominal);

        for (name, got, exp) in [
            ("prefix", &prefix, &sp["prefix_score"]),
            ("suffix", &suffix, &sp["suffix_score"]),
        ] {
            let exp = exp.as_array().unwrap();
            for nm in 1..got.len().min(exp.len()) {
                let d = (got[nm] as f64 - exp[nm].as_f64().unwrap()).abs();
                if d > worst {
                    worst = d;
                    worst_ctx = format!("scan {scan} {name}[{nm}]");
                }
            }
        }
        n += 1;
    }
    eprintln!("end-to-end: {n} spectra raw→preprocess→score; worst Δ = {worst:.4e} ({worst_ctx})");
    assert!(
        worst < 2e-3,
        "end-to-end mismatch — worst Δ {worst:.4e} at {worst_ctx}"
    );
}
