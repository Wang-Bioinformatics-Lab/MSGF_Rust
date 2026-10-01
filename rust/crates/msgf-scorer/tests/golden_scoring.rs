//! Spectrum preparation and RawScore against MS-GF+'s own outputs on F13, for both the
//! QExactive and the HighRes model.
//!
//! - **Preparation**: raw F13 peaks → [`msgf_scorer::preprocess`] must reproduce MS-GF+'s prepared
//!   peak list (`f13_scored_spectrum*.golden.json`): same count, m/z within 1e-3, ranks exactly.
//! - **RawScore**: raw F13 spectrum → [`msgf_scorer::prepare`] → the node + edge match score of the
//!   golden peptide must equal MS-GF+'s `full_score` (`f13_rawscore*.golden.json`) exactly.
//!
//! These goldens come from running MS-GF+ and are not committed (`LICENSING.md` §2); build them
//! with `validation/reference/build_all_golden.sh --with-java`. Skipped when absent. Only black-box
//! outputs are compared: nothing here depends on how MS-GF+ computes them.

use msgf_chem::peptide;
use msgf_scorer::{
    match_and_terminal_parts, prepare, preprocess, Candidate, Cleavage, Ms2, PreprocessParams,
    Residue, ScoreModel,
};
use serde_json::Value;
use std::collections::HashMap;
use std::path::PathBuf;

fn repo(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

/// Raw F13 spectra by scan: (precursor m/z, peaks).
fn f13() -> Option<HashMap<i64, (f64, Vec<(f64, f64)>)>> {
    let mgf = repo("validation/data/spectra/F13.mgf");
    let specs = msgf_io::read_mgf_file(&mgf).ok()?;
    Some(
        specs
            .iter()
            .filter_map(|s| {
                let scan = s.scan.as_ref()?.parse::<i64>().ok()?;
                let peaks = s.peaks.iter().map(|p| (p.mz, p.intensity)).collect();
                Some((scan, (s.precursor_mz?, peaks)))
            })
            .collect(),
    )
}

fn load(rel: &str) -> Option<Value> {
    let p = repo(rel);
    p.exists()
        .then(|| serde_json::from_str(&std::fs::read_to_string(&p).unwrap()).unwrap())
}

fn check(model_file: &str, scored_golden: &str, raw_golden: &str) {
    let param = repo(model_file);
    let (Some(ss), Some(rs), Some(raw)) = (load(scored_golden), load(raw_golden), f13()) else {
        eprintln!("skip: {scored_golden} / {raw_golden} / F13.mgf absent");
        return;
    };
    if !param.exists() {
        eprintln!("skip: {model_file} absent");
        return;
    }
    let file = msgf_scorer::read_param_file(&param).unwrap();
    let params = PreprocessParams::from_param(&file);
    let model = ScoreModel::from_file(&param).unwrap();

    // (1) preparation
    let (mut prep_ok, mut prep_total) = (0, 0);
    for s in ss["spectra"].as_array().unwrap() {
        let scan = s["scan"].as_i64().unwrap();
        let charge = s["charge"].as_i64().unwrap() as i32;
        let parent_mass = s["precursor_mass"].as_f64().unwrap() as f32;
        let Some((_, peaks)) = raw.get(&scan) else {
            continue;
        };
        let peaks: Vec<(f32, f32)> = peaks.iter().map(|&(m, i)| (m as f32, i as f32)).collect();
        let got = preprocess(&params, charge, parent_mass, &peaks);
        let exp = s["peaks"].as_array().unwrap();
        prep_total += 1;
        let same = got.len() == exp.len()
            && got.iter().zip(exp).all(|(g, e)| {
                g.rank == e[2].as_i64().unwrap() as i32
                    && (g.mz as f64 - e[0].as_f64().unwrap()).abs() <= 1e-3
            });
        prep_ok += same as i32;
    }

    // (2) RawScore (node + edge, no terminus terms)
    let (mut full_ok, mut full_total) = (0, 0);
    let mut bad = Vec::new();
    let charges: HashMap<i64, i32> = ss["spectra"]
        .as_array()
        .unwrap()
        .iter()
        .map(|s| {
            (
                s["scan"].as_i64().unwrap(),
                s["charge"].as_i64().unwrap() as i32,
            )
        })
        .collect();
    for r in rs["spectra"].as_array().unwrap() {
        let scan = r["scan"].as_i64().unwrap();
        let (Some(&charge), Some((mz, peaks))) = (charges.get(&scan), raw.get(&scan)) else {
            continue;
        };
        let pep = r["peptide"].as_str().unwrap();
        let residues: Vec<Residue> = peptide::parse(pep)
            .unwrap()
            .iter()
            .map(|x| Residue {
                letter: x.aa,
                delta: x.mod_delta,
            })
            .collect();
        let ms2 = Ms2 {
            peaks,
            precursor_mz: *mz,
            charge,
        };
        let Some(prep) = prepare(&model, &ms2) else {
            continue;
        };
        let cand = Candidate {
            residues: &residues,
            n_term_credit: true,
        };
        let (m, _) = match_and_terminal_parts(&prep, &cand, &Cleavage::disabled()).unwrap();
        let want = r["full_score"].as_i64().unwrap() as i32;
        full_total += 1;
        if m == want {
            full_ok += 1;
        } else if bad.len() < 10 {
            bad.push(format!("scan {scan} {pep}: {m} vs {want}"));
        }
    }
    eprintln!(
        "{model_file}: preparation {prep_ok}/{prep_total}, RawScore {full_ok}/{full_total} {bad:?}"
    );
    assert_eq!(prep_ok, prep_total, "prepared peak lists");
    assert_eq!(full_ok, full_total, "RawScore (node + edge)");
}

#[test]
fn qexactive_scoring_matches_msgf() {
    check(
        "validation/data/models/HCD_QExactive_Tryp.param",
        "validation/golden/rawscore/f13_scored_spectrum.golden.json",
        "validation/golden/rawscore/f13_rawscore.golden.json",
    );
}

#[test]
fn highres_scoring_matches_msgf() {
    check(
        "validation/data/models/HCD_HighRes_Tryp.param",
        "validation/golden/rawscore/f13_scored_spectrum_highres.golden.json",
        "validation/golden/rawscore/f13_rawscore_highres.golden.json",
    );
}
