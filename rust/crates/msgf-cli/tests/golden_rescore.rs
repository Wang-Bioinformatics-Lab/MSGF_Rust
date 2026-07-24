//! End-to-end test of `msgf rescore` against MS-GF+'s own `f13_specprob.golden.json`.
//!
//! Builds a PSM list (scan, peptide, charge) from the golden, runs the compiled `msgf` binary with
//! the HighRes model + iPRG DB-composition amino-acid probabilities + oxidation-on-M (the exact
//! configuration the F13 search used), and asserts the CLI reproduces MS-GF+ RawScore and
//! DeNovoScore exactly and SpecEValue within tolerance. Skipped if goldens/model/data are absent.

use serde_json::Value;
use std::collections::HashMap;
use std::path::PathBuf;
use std::process::Command;

fn repo(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

/// iPRG-2013 human FASTA composition (DBScanner.setAminoAcidProbabilities).
const IPRG: &[(&str, f64)] = &[
    ("G", 0.065416),
    ("A", 0.069428),
    ("S", 0.083673),
    ("P", 0.063069),
    ("V", 0.059874),
    ("T", 0.053651),
    ("C", 0.022017),
    ("L", 0.098978),
    ("I", 0.043441),
    ("N", 0.035999),
    ("D", 0.048123),
    ("Q", 0.048243),
    ("K", 0.057485),
    ("E", 0.071678),
    ("M", 0.021972),
    ("H", 0.025877),
    ("F", 0.035796),
    ("R", 0.056718),
    ("Y", 0.026244),
    ("W", 0.012320),
];

#[test]
fn rescore_cli_matches_golden() {
    let sp_path = repo("validation/golden/rawscore/f13_specprob.golden.json");
    let mgf = repo("validation/data/spectra/F13.mgf");
    let param = repo("validation/data/models/HCD_HighRes_Tryp.param");
    if !sp_path.exists() || !mgf.exists() || !param.exists() {
        eprintln!("skip: goldens/model/data absent");
        return;
    }
    let sp: Value = serde_json::from_str(&std::fs::read_to_string(&sp_path).unwrap()).unwrap();
    let entries = sp["spectra"].as_array().unwrap();

    // Write the PSM list and the aa-probabilities file into the test temp dir.
    let tmp = PathBuf::from(env!("CARGO_TARGET_TMPDIR"));
    let psms = tmp.join("psms.tsv");
    let probs = tmp.join("iprg.tsv");
    let mut psm_txt = String::from("scan\tpeptide\tcharge\n");
    let mut expected: HashMap<String, (i64, i64, f64)> = HashMap::new();
    for e in entries {
        let scan = e["scan"].as_i64().unwrap().to_string();
        psm_txt.push_str(&format!(
            "{}\t{}\t{}\n",
            scan,
            e["peptide"].as_str().unwrap(),
            e["charge"].as_i64().unwrap()
        ));
        expected.insert(
            scan,
            (
                e["raw_score"].as_i64().unwrap(),
                e["denovo_score"].as_i64().unwrap(),
                e["spec_prob"].as_f64().unwrap(),
            ),
        );
    }
    std::fs::write(&psms, psm_txt).unwrap();
    let probs_txt: String = IPRG.iter().map(|(r, p)| format!("{r}\t{p}\n")).collect();
    std::fs::write(&probs, probs_txt).unwrap();

    // Run the compiled binary.
    let out = tmp.join("rescored.tsv");
    let status = Command::new(env!("CARGO_BIN_EXE_msgf"))
        .args(["rescore", "-s"])
        .arg(&mgf)
        .arg("-p")
        .arg(&param)
        .arg("-i")
        .arg(&psms)
        .arg("--aa-probs")
        .arg(&probs)
        .arg("--ox-m")
        .arg("-o")
        .arg(&out)
        .status()
        .expect("run msgf rescore");
    assert!(status.success(), "msgf rescore exited with failure");

    // Compare the output rows to the golden.
    let text = std::fs::read_to_string(&out).unwrap();
    let mut lines = text.lines();
    let header = lines.next().unwrap();
    let col: HashMap<&str, usize> = header
        .split('\t')
        .enumerate()
        .map(|(i, c)| (c, i))
        .collect();
    let (mut raw_ok, mut den_ok, mut spec_ok, mut n) = (0, 0, 0, 0);
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        let scan = f[col["scan"]];
        let (g_raw, g_den, g_spec) = expected[scan];
        let raw: i64 = f[col["raw_score"]].parse().unwrap();
        let den: i64 = f[col["denovo_score"]].parse().unwrap();
        let spec: f64 = f[col["spec_evalue"]].parse().unwrap();
        n += 1;
        raw_ok += (raw == g_raw) as i32;
        den_ok += (den == g_den) as i32;
        let dlog = if spec > 0.0 && g_spec > 0.0 {
            (spec / g_spec).log10().abs()
        } else {
            f64::INFINITY
        };
        spec_ok += (dlog <= 0.05) as i32;
    }
    eprintln!("rescore CLI: raw {raw_ok}/{n}, denovo {den_ok}/{n}, spec {spec_ok}/{n}");
    assert_eq!(n, entries.len() as i32, "every PSM should be scored");
    assert_eq!(raw_ok, n, "RawScore must match MS-GF+ exactly");
    assert_eq!(den_ok, n, "DeNovoScore must match MS-GF+ exactly");
    assert_eq!(spec_ok, n, "SpecEValue must match MS-GF+ within tolerance");
}
