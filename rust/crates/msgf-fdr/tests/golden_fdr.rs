//! `PLAN2.md` TD-2 Gate 1: reproduce MS-GF+'s `QValue` and `PepQValue` columns exactly.
//!
//! The oracle is `validation/golden/iprg2013_F13.golden.json`. Being MS-GF+-derived it is **not
//! committed** — regenerate it with `validation/reference/generate_golden.sh` (jar + spectra +
//! FASTA); this test skips until you do.
//!
//! The golden holds MS-GF+'s `-unroll 1` output: 4133 rows, one per *protein occurrence*. FDR
//! counts *matches*, so the rows are first rolled back up into 1610 unique PSMs keyed by
//! `(spec_file, scan, charge, peptide)`, with each match carrying every protein it hit — that list
//! is what decides decoy status (`PLAN2.md` §1.3: decoy iff **every** occurrence is a decoy).

use msgf_fdr::{is_decoy_match, peptide_key, PsmRecord, TargetDecoyAnalysis};
use serde_json::Value;
use std::collections::BTreeMap;
use std::path::PathBuf;

fn repo(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

const DECOY_PREFIX: &str = "XXX_";

#[test]
fn f13_q_values_match_msgfplus() {
    let path = repo("validation/golden/iprg2013_F13.golden.json");
    if !path.exists() {
        eprintln!("skip: {} absent", path.display());
        return;
    }
    let g: Value = serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
    let rows = g["psms"].as_array().expect("psms array");
    assert_eq!(rows.len(), 4133, "golden shape changed");

    // Roll the unrolled rows back up into one record per match.
    struct Match {
        score: f32,
        peptide: String,
        proteins: Vec<String>,
        want_q: f32,
        want_pep_q: f32,
    }
    let mut by_key: BTreeMap<(String, String, i64, String), Match> = BTreeMap::new();
    for r in rows {
        let key = (
            r["spec_file"].as_str().unwrap().to_string(),
            r["scan"].as_str().unwrap().to_string(),
            r["charge"].as_i64().unwrap(),
            r["peptide"].as_str().unwrap().to_string(),
        );
        let protein = r["protein"].as_str().unwrap().to_string();
        by_key
            .entry(key)
            .and_modify(|m| m.proteins.push(protein.clone()))
            .or_insert_with(|| Match {
                score: r["spec_evalue"].as_f64().unwrap() as f32,
                peptide: r["peptide"].as_str().unwrap().to_string(),
                proteins: vec![protein],
                want_q: r["qvalue"].as_f64().unwrap() as f32,
                want_pep_q: r["pep_qvalue"].as_f64().unwrap() as f32,
            });
    }
    let matches: Vec<Match> = by_key.into_values().collect();
    assert_eq!(matches.len(), 1610, "expected 1610 unique PSMs");

    let records: Vec<PsmRecord> = matches
        .iter()
        .map(|m| PsmRecord {
            score: m.score,
            peptide: peptide_key(&m.peptide),
            is_decoy: is_decoy_match(m.proteins.iter().map(String::as_str), DECOY_PREFIX),
        })
        .collect();
    let n_decoy = records.iter().filter(|r| r.is_decoy).count();
    assert_eq!(
        n_decoy, 765,
        "decoy labelling changed (PLAN2 §4 measured 765)"
    );

    let tda = TargetDecoyAnalysis::new(&records, 1.0);

    let (mut q_ok, mut pep_ok) = (0usize, 0usize);
    let (mut q_bad, mut pep_bad) = (Vec::new(), Vec::new());
    for (m, rec) in matches.iter().zip(&records) {
        let got_q = tda.psm_q_value(rec.score);
        if got_q == m.want_q {
            q_ok += 1;
        } else if q_bad.len() < 5 {
            q_bad.push(format!("{}: got {got_q} want {}", m.peptide, m.want_q));
        }
        let got_pep = tda.pep_q_value(&rec.peptide, rec.score);
        if got_pep == m.want_pep_q {
            pep_ok += 1;
        } else if pep_bad.len() < 5 {
            pep_bad.push(format!(
                "{}: got {got_pep} want {}",
                m.peptide, m.want_pep_q
            ));
        }
    }

    assert_eq!(
        q_ok,
        matches.len(),
        "QValue mismatches (first few): {q_bad:?}"
    );
    assert_eq!(
        pep_ok,
        matches.len(),
        "PepQValue mismatches (first few): {pep_bad:?}"
    );
    eprintln!(
        "F13 target-decoy: {}/{} QValue and PepQValue exact vs MS-GF+",
        q_ok,
        matches.len()
    );
}

/// The F13 corpus is degenerate as an FDR oracle (`PLAN2.md` §4): MS-GF+ itself reports q = 1 for
/// 4132 of 4133 rows. This test pins that fact so the gate above is never mistaken for proof that
/// the estimator is right across its range — that needs PLAN2's TD-2 Gate 2 Java probe.
#[test]
fn f13_is_a_degenerate_fdr_oracle() {
    let path = repo("validation/golden/iprg2013_F13.golden.json");
    if !path.exists() {
        eprintln!("skip: {} absent", path.display());
        return;
    }
    let g: Value = serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
    let rows = g["psms"].as_array().unwrap();
    let distinct: std::collections::BTreeSet<String> = rows
        .iter()
        .map(|r| format!("{:?}", r["qvalue"].as_f64().unwrap()))
        .collect();
    assert_eq!(
        distinct.len(),
        2,
        "F13 was expected to yield only q = 0 and q = 1; got {distinct:?}"
    );
    let n_one = rows
        .iter()
        .filter(|r| r["qvalue"].as_f64().unwrap() == 1.0)
        .count();
    assert_eq!(n_one, 4132, "F13 q-value distribution changed");
}
