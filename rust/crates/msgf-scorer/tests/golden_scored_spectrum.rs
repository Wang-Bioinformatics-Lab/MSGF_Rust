//! Validates the scored-spectrum node scoring against MS-GF+'s own `prefixScore`/`suffixScore`
//! arrays (`validation/golden/rawscore/f13_scored_spectrum.golden.json`, produced by the Java
//! ScoredSpectrumDumper). We feed MS-GF+'s *preprocessed* peaks into the Rust scorer, so this
//! isolates the node-scoring integration (ion m/z + peak lookup + segment/partition selection +
//! node_score) from spectrum preprocessing, which is ported separately.
//!
//! Skipped when the golden or the gitignored `.param`/data are absent.

use msgf_scorer::scored_spectrum::{RankedPeak, ScoredSpectrum};
use serde_json::Value;
use std::path::PathBuf;

fn repo(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

#[test]
fn scored_spectrum_matches_msgf() {
    let gpath = repo("validation/golden/rawscore/f13_scored_spectrum.golden.json");
    let param = repo("validation/data/models/HCD_QExactive_Tryp.param");
    if !gpath.exists() || !param.exists() {
        eprintln!("skip: golden or .param absent ({gpath:?})");
        return;
    }
    let g: Value = serde_json::from_str(&std::fs::read_to_string(&gpath).unwrap()).unwrap();
    let model = msgf_scorer::read_param_file(&param).unwrap();

    let mut worst = 0.0f64;
    let mut worst_ctx = String::new();
    let mut n_spectra = 0;
    let mut n_values = 0u64;

    for sp in g["spectra"].as_array().unwrap() {
        let charge = sp["charge"].as_i64().unwrap() as i32;
        let parent_mass = sp["precursor_mass"].as_f64().unwrap() as f32;
        let pep_nominal = sp["peptide_mass_nominal"].as_i64().unwrap() as i32;
        let scan = sp["scan"].as_i64().unwrap_or(-1);

        let peaks: Vec<RankedPeak> = sp["peaks"]
            .as_array()
            .unwrap()
            .iter()
            .map(|p| RankedPeak {
                mz: p[0].as_f64().unwrap() as f32,
                intensity: p[1].as_f64().unwrap() as f32,
                rank: p[2].as_i64().unwrap() as i32,
            })
            .collect();

        let ss = ScoredSpectrum::from_ranked_peaks(&model, charge, parent_mass, peaks);
        let (prefix, suffix) = ss.prefix_suffix_scores(pep_nominal);

        for (name, got, exp) in [
            ("prefix", &prefix, &sp["prefix_score"]),
            ("suffix", &suffix, &sp["suffix_score"]),
        ] {
            let exp = exp.as_array().unwrap();
            let n = got.len().min(exp.len());
            for nm in 1..n {
                let e = exp[nm].as_f64().unwrap();
                let d = (got[nm] as f64 - e).abs();
                if d > worst {
                    worst = d;
                    worst_ctx =
                        format!("scan {scan} c{charge} {name}[{nm}]: got {} vs {e}", got[nm]);
                }
                n_values += 1;
            }
        }
        n_spectra += 1;
    }

    eprintln!(
        "checked {n_spectra} spectra, {n_values} node scores; worst Δ = {worst:.4e} ({worst_ctx})"
    );
    assert!(
        worst < 2e-3,
        "scored-spectrum mismatch — worst Δ {worst:.4e} at {worst_ctx}"
    );
}
