//! Validates msgf-chem against the frozen chemistry golden in `validation/golden/chemistry/`.
//! That golden is derived from authoritative atomic masses and guarded against published
//! peptide calibrants, so matching it proves the Rust mass model is correct.

use serde_json::Value;
use std::path::PathBuf;

fn golden_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../validation/golden/chemistry")
        .canonicalize()
        .expect("validation/golden/chemistry must exist (committed)")
}

fn load(name: &str) -> Value {
    let p = golden_dir().join(name);
    serde_json::from_str(&std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("read {p:?}: {e}")))
        .unwrap()
}

fn approx(a: f64, b: f64, tol: f64, ctx: &str) {
    assert!(
        (a - b).abs() <= tol,
        "{ctx}: {a} vs {b} (Δ={:.3e} > {tol:.0e})",
        (a - b).abs()
    );
}

#[test]
fn constants_match() {
    let g = load("constants.golden.json");
    let atoms = &g["atoms"];
    approx(msgf_chem::mass::H, atoms["H"].as_f64().unwrap(), 1e-9, "H");
    approx(msgf_chem::mass::C, atoms["C"].as_f64().unwrap(), 1e-9, "C");
    approx(msgf_chem::mass::N, atoms["N"].as_f64().unwrap(), 1e-9, "N");
    approx(msgf_chem::mass::O, atoms["O"].as_f64().unwrap(), 1e-9, "O");
    approx(msgf_chem::mass::S, atoms["S"].as_f64().unwrap(), 1e-9, "S");
    approx(msgf_chem::mass::P, atoms["P"].as_f64().unwrap(), 1e-9, "P");
    approx(
        msgf_chem::mass::PROTON,
        g["proton"].as_f64().unwrap(),
        1e-9,
        "proton",
    );
    approx(
        msgf_chem::mass::WATER,
        g["H2O"].as_f64().unwrap(),
        1e-9,
        "H2O",
    );
}

#[test]
fn residue_masses_match() {
    let g = load("residue_masses.golden.json");
    let residues = g["residues"].as_object().unwrap();
    for (aa, rec) in residues {
        let expect = rec["mass"].as_f64().unwrap();
        let got = msgf_chem::residue_mass(aa.as_bytes()[0]).unwrap();
        approx(got, expect, 1e-6, &format!("residue {aa}"));
    }
    assert_eq!(residues.len(), 20, "expected 20 standard residues");
}

#[test]
fn peptide_masses_match() {
    let g = load("peptide_masses.golden.json");
    for pep in g["peptides"].as_array().unwrap() {
        let seq = pep["sequence"].as_str().unwrap();
        let neutral = msgf_chem::peptide_neutral_mass(seq).unwrap();
        approx(
            neutral,
            pep["neutral_mass"].as_f64().unwrap(),
            1e-4,
            &format!("{seq} neutral"),
        );
        for (z, key) in [(1u32, "1+"), (2, "2+"), (3, "3+")] {
            let got = msgf_chem::mz(neutral, z);
            approx(
                got,
                pep["mz"][key].as_f64().unwrap(),
                1e-4,
                &format!("{seq} {key}"),
            );
        }
    }
}

#[test]
fn fragment_ions_match() {
    let g = load("fragment_ions.golden.json");
    for fp in g["peptides"].as_array().unwrap() {
        let seq = fp["peptide"].as_str().unwrap();
        let b = msgf_chem::b_ions(seq);
        for (k, ion) in fp["b_ions"].as_array().unwrap().iter().enumerate() {
            approx(
                b[k].mz1,
                ion["z1"].as_f64().unwrap(),
                1e-4,
                &format!("{seq} b{} z1", b[k].index),
            );
            approx(
                b[k].mz2,
                ion["z2"].as_f64().unwrap(),
                1e-4,
                &format!("{seq} b{} z2", b[k].index),
            );
        }
        let y = msgf_chem::y_ions(seq);
        for (k, ion) in fp["y_ions"].as_array().unwrap().iter().enumerate() {
            approx(
                y[k].mz1,
                ion["z1"].as_f64().unwrap(),
                1e-4,
                &format!("{seq} y{} z1", y[k].index),
            );
            approx(
                y[k].mz2,
                ion["z2"].as_f64().unwrap(),
                1e-4,
                &format!("{seq} y{} z2", y[k].index),
            );
        }
    }
}
