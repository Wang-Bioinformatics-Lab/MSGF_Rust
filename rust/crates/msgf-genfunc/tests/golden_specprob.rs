//! Validates the generating function (DeNovoScore + SpecEValue) against MS-GF+'s own
//! `f13_specprob.golden.json` — built on the **HighRes** model (the one the F13 search used), with
//! DB-composition amino-acid probabilities, oxidised M and the golden's isotope-error sink range.
//! Raw F13 spectrum → `prepare` → `score_distribution` → DeNovoScore and the tail at MS-GF+'s
//! RawScore. Only black-box outputs are compared. Skipped if goldens/model/data are absent (the
//! golden comes from running MS-GF+ and is not committed; see `LICENSING.md` §2).

use msgf_genfunc::{score_distribution, Cleavage, NullModel};
use msgf_scorer::{prepare, Ms2, ScoreModel};
use serde_json::Value;
use std::collections::HashMap;
use std::path::PathBuf;

fn repo(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

/// iPRG-2013 human FASTA composition: residue → probability.
fn iprg_probs() -> HashMap<u8, f64> {
    [
        (b'G', 0.065416),
        (b'A', 0.069428),
        (b'S', 0.083673),
        (b'P', 0.063069),
        (b'V', 0.059874),
        (b'T', 0.053651),
        (b'C', 0.022017),
        (b'L', 0.098978),
        (b'I', 0.043441),
        (b'N', 0.035999),
        (b'D', 0.048123),
        (b'Q', 0.048243),
        (b'K', 0.057485),
        (b'E', 0.071678),
        (b'M', 0.021972),
        (b'H', 0.025877),
        (b'F', 0.035796),
        (b'R', 0.056718),
        (b'Y', 0.026244),
        (b'W', 0.012320),
    ]
    .into_iter()
    .collect()
}

#[test]
fn generating_function_matches_golden() {
    let ss_path = repo("validation/golden/rawscore/f13_scored_spectrum.golden.json");
    let sp_path = repo("validation/golden/rawscore/f13_specprob.golden.json");
    let mgf = repo("validation/data/spectra/F13.mgf");
    let param = repo("validation/data/models/HCD_HighRes_Tryp.param"); // the model the F13 search used
    if !ss_path.exists() || !sp_path.exists() || !mgf.exists() || !param.exists() {
        eprintln!("skip: goldens/model/data absent");
        return;
    }
    let ss: Value = serde_json::from_str(&std::fs::read_to_string(&ss_path).unwrap()).unwrap();
    let sp: Value = serde_json::from_str(&std::fs::read_to_string(&sp_path).unwrap()).unwrap();
    let model = ScoreModel::from_file(&param).unwrap();

    // precursor charge per scan, from the scored-spectrum golden
    let charge_of: HashMap<i64, i32> = ss["spectra"]
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
    // raw spectra by scan: (precursor m/z, peaks)
    let raw: HashMap<i64, (f64, Vec<(f64, f64)>)> = msgf_io::read_mgf_file(&mgf)
        .unwrap()
        .iter()
        .filter_map(|s| {
            let scan = s.scan.as_deref()?.parse::<i64>().ok()?;
            Some((
                scan,
                (
                    s.precursor_mz?,
                    s.peaks.iter().map(|p| (p.mz, p.intensity)).collect(),
                ),
            ))
        })
        .collect();

    let probs = iprg_probs();
    // variable oxidation on M (iprg-2013_Mods.txt): an extra residue at M's DB probability
    let base = NullModel::from_probs(
        |l| probs[&l],
        &[(b'M', 15.994915)],
        Cleavage::trypsin(),
        (0, 1),
    );

    let (mut denovo_ok, mut spec_ok, mut total) = (0, 0, 0);
    let mut worst = String::new();
    for e in sp["spectra"].as_array().unwrap() {
        let scan = e["scan"].as_i64().unwrap();
        let (Some(&charge), Some((mz, peaks))) = (charge_of.get(&scan), raw.get(&scan)) else {
            continue;
        };
        let range = e["mass_index_range"].as_array().unwrap();
        let (lo_sink, hi_sink) = (
            range[0].as_i64().unwrap() as i32,
            range[1].as_i64().unwrap() as i32,
        );
        let g_raw = e["raw_score"].as_i64().unwrap() as i32;
        let g_denovo = e["denovo_score"].as_i64().unwrap() as i32;
        let g_spec = e["spec_prob"].as_f64().unwrap();
        total += 1;

        let ms2 = Ms2 {
            peaks,
            precursor_mz: *mz,
            charge,
        };
        let Some(prep) = prepare(&model, &ms2) else {
            continue;
        };
        // The golden's sinks, as an isotope range around this spectrum's N0.
        let mut null = base.clone();
        null.isotope = (prep.n0() - hi_sink, prep.n0() - lo_sink);
        let Some(tail) = score_distribution(&prep, &null, None) else {
            continue;
        };

        let my_denovo = tail.best_possible();
        let my_spec = tail.tail_mass(g_raw);
        let dlog = if g_spec > 0.0 && my_spec > 0.0 {
            (my_spec / g_spec).log10()
        } else {
            f64::NAN
        };
        denovo_ok += (my_denovo == g_denovo) as i32;
        spec_ok += (dlog.abs() <= 1e-4) as i32;
        if (my_denovo != g_denovo || dlog.abs() > 1e-4) && worst.len() < 900 {
            worst.push_str(&format!(
                "\n  scan {scan}: denovo {my_denovo} vs {g_denovo} | spec {my_spec:.4e} vs {g_spec:.4e} (Δlog {dlog:.3})"
            ));
        }
    }
    eprintln!("DeNovoScore exact: {denovo_ok}/{total}; SpecEValue: {spec_ok}/{total}{worst}");
    assert!(total >= 25, "expected ~30 matched spectra");
    assert_eq!(denovo_ok, total, "DeNovoScore must match MS-GF+ exactly");
    assert_eq!(
        spec_ok, total,
        "SpecEValue must match MS-GF+ to |Δlog10| <= 1e-4"
    );
}
