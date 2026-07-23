//! Validates `msgf_scorer::preprocess` against MS-GF+'s own preprocessed spectra
//! (`validation/golden/rawscore/f13_scored_spectrum.golden.json`, produced by the Java
//! ScoredSpectrumDumper). For each golden spectrum we read its raw peaks from `F13.mgf`
//! (matched by scan), run `preprocess`, and require the result to match the golden peak list:
//! same count, per-peak m/z within 1e-3 and intensity within 1e-2, and **rank exactly**.
//!
//! Skipped when the golden, the `.param`, or the gitignored `F13.mgf` are absent.

use msgf_scorer::preprocess::preprocess;
use serde_json::Value;
use std::collections::HashMap;
use std::path::PathBuf;

fn repo(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

#[test]
fn preprocess_matches_msgf() {
    let gpath = repo("validation/golden/rawscore/f13_scored_spectrum.golden.json");
    let param = repo("validation/data/models/HCD_QExactive_Tryp.param");
    let mgf = repo("validation/data/spectra/F13.mgf");
    if !gpath.exists() || !param.exists() || !mgf.exists() {
        eprintln!("skip: golden, .param, or F13.mgf absent ({gpath:?}, {mgf:?})");
        return;
    }

    let g: Value = serde_json::from_str(&std::fs::read_to_string(&gpath).unwrap()).unwrap();
    let model = msgf_scorer::read_param_file(&param).unwrap();

    // Index raw MGF spectra by scan number (SCANS=).
    let specs = msgf_io::read_mgf_file(&mgf).unwrap();
    let by_scan: HashMap<i64, Vec<(f32, f32)>> = specs
        .iter()
        .filter_map(|s| {
            let scan = s.scan.as_ref()?.parse::<i64>().ok()?;
            let peaks = s
                .peaks
                .iter()
                .map(|p| (p.mz as f32, p.intensity as f32))
                .collect();
            Some((scan, peaks))
        })
        .collect();

    let mut n_spectra = 0;
    let mut n_peaks = 0u64;
    let mut mismatches: Vec<String> = Vec::new();
    let mut worst_mz = 0.0f64;
    let mut worst_int = 0.0f64;

    for sp in g["spectra"].as_array().unwrap() {
        let scan = sp["scan"].as_i64().unwrap();
        let charge = sp["charge"].as_i64().unwrap() as i32;
        let precursor_mass = sp["precursor_mass"].as_f64().unwrap() as f32;

        let Some(raw) = by_scan.get(&scan) else {
            mismatches.push(format!("scan {scan}: not found in F13.mgf"));
            continue;
        };

        let got = preprocess(&model, charge, precursor_mass, raw);
        let exp = sp["peaks"].as_array().unwrap();

        if got.len() != exp.len() {
            mismatches.push(format!(
                "scan {scan} c{charge}: peak count {} vs golden {}",
                got.len(),
                exp.len()
            ));
            continue;
        }

        for (i, (gp, ep)) in got.iter().zip(exp.iter()).enumerate() {
            let emz = ep[0].as_f64().unwrap();
            let eint = ep[1].as_f64().unwrap();
            let erank = ep[2].as_i64().unwrap() as i32;

            let dmz = (gp.mz as f64 - emz).abs();
            let dint = (gp.intensity as f64 - eint).abs();
            worst_mz = worst_mz.max(dmz);
            worst_int = worst_int.max(dint);

            if (dmz > 1e-3 || dint > 1e-2 || gp.rank != erank) && mismatches.len() < 40 {
                mismatches.push(format!(
                    "scan {scan} c{charge} peak[{i}]: got (mz={:.5}, int={:.4}, rank={}) \
                         vs golden (mz={emz:.5}, int={eint:.4}, rank={erank})  \
                         [Δmz={dmz:.2e} Δint={dint:.2e} Δrank={}]",
                    gp.mz,
                    gp.intensity,
                    gp.rank,
                    gp.rank - erank
                ));
            }
            n_peaks += 1;
        }
        n_spectra += 1;
    }

    eprintln!(
        "checked {n_spectra} spectra, {n_peaks} peaks; worst Δmz = {worst_mz:.3e}, worst Δint = {worst_int:.3e}"
    );
    assert!(
        mismatches.is_empty(),
        "preprocess mismatch ({} problem(s)):\n{}",
        mismatches.len(),
        mismatches.join("\n")
    );
}
