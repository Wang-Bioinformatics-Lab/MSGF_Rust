//! The CLI must score with **no model on the command line** — that is what makes the tool usable
//! (and shippable) without any UC-licensed `.param`.
//!
//! Self-contained on purpose: the spectrum is synthesised here, so this test runs on a clean
//! checkout with `validation/data/` absent.

use msgf_chem::{mass, peptide};
use std::path::PathBuf;
use std::process::Command;

const PEPTIDE: &str = "SAMPLERPEPTIDEK";

/// A b/y spectrum for `PEPTIDE` at charge 2, intense enough to score.
fn synthetic_mgf() -> String {
    let residues = peptide::parse(PEPTIDE).unwrap();
    let acc = peptide::accurate_prefix_masses(&residues);
    let pep_mass = acc[acc.len() - 1] + mass::WATER;

    let mut peaks: Vec<(f64, f64)> = Vec::new();
    for (k, &prefix) in acc[..acc.len() - 1].iter().enumerate() {
        let suffix = acc[acc.len() - 1] - prefix;
        peaks.push((
            suffix + mass::WATER + mass::PROTON,
            2000.0 - k as f64 * 10.0,
        ));
        peaks.push((prefix + mass::PROTON, 900.0 - k as f64 * 10.0));
    }
    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let mut s = String::from("BEGIN IONS\n");
    s.push_str(&format!(
        "TITLE=synthetic\nPEPMASS={}\nCHARGE=2+\nSCANS=1\n",
        (pep_mass + 2.0 * mass::PROTON) / 2.0
    ));
    for (mz, it) in peaks {
        s.push_str(&format!("{mz:.5} {it:.1}\n"));
    }
    s.push_str("END IONS\n");
    s
}

#[test]
fn rescore_uses_the_bundled_model_when_no_param_is_given() {
    let tmp = PathBuf::from(env!("CARGO_TARGET_TMPDIR"));
    let mgf = tmp.join("bundled_default.mgf");
    let psms = tmp.join("bundled_default_psms.tsv");
    std::fs::write(&mgf, synthetic_mgf()).unwrap();
    std::fs::write(&psms, format!("scan\tpeptide\tcharge\n1\t{PEPTIDE}\t2\n")).unwrap();

    let out = Command::new(env!("CARGO_BIN_EXE_msgf"))
        .args(["rescore", "-s"])
        .arg(&mgf)
        .arg("-i")
        .arg(&psms)
        .output()
        .expect("run msgf rescore");
    assert!(
        out.status.success(),
        "rescore failed without --param:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );

    // The run must say which model produced the numbers.
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains(msgf_scorer::bundled::NAME) && stderr.contains("MassIVE-KB"),
        "the bundled model should be announced on stderr, got: {stderr}"
    );

    let stdout = String::from_utf8_lossy(&out.stdout);
    let mut lines = stdout.lines();
    let header: Vec<&str> = lines.next().expect("header").split('\t').collect();
    let row: Vec<&str> = lines.next().expect("one scored PSM").split('\t').collect();
    let col = |c: &str| row[header.iter().position(|h| *h == c).unwrap()];

    let raw: i64 = col("raw_score").parse().unwrap();
    let denovo: i64 = col("denovo_score").parse().unwrap();
    let evalue: f64 = col("spec_evalue").parse().unwrap();
    assert!(
        raw > 0,
        "a complete b/y ladder should score positive, got {raw}"
    );
    assert!(
        denovo >= raw,
        "DeNovoScore is the best path, so it cannot be below RawScore"
    );
    assert!(
        evalue > 0.0 && evalue < 1e-3,
        "SpecEValue should be small and finite for a full ladder, got {evalue}"
    );
}
