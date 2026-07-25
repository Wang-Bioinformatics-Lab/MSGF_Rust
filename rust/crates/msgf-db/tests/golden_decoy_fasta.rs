//! `plans/PLAN2.md` TD-1 gate: our target-decoy FASTA must be **byte-identical** to the one MS-GF+'s
//! `-tda 1` produced.
//!
//! The reference `.revCat.fasta` files live in the gitignored `validation/data/` (UC-licensed, see
//! `CLAUDE.md`), so this test skips gracefully when they are absent — a fresh clone still passes.
//! Run `validation/fetch_reference_data.sh` (add `--full` for the human database) to enable it.
//!
//! The two references disagree on line endings because they were generated on different platforms:
//! MS-GF+ re-emits every line with the JVM's `line.separator`, so the Windows-generated
//! `Tryp_Pig_Bov` file is CRLF throughout and the Linux-generated `iprg2013_human` one is LF (its
//! CRLF input is normalised down). [`LineSep`] selects which, and both are checked here.

use msgf_db::decoy::{write_decoy_database, DecoyOptions, LineSep, Output};
use std::path::PathBuf;

fn repo(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

fn check(stem: &str, line_sep: LineSep) {
    let src = repo(&format!("validation/data/fasta/{stem}.fasta"));
    let want = repo(&format!("validation/data/fasta/{stem}.revCat.fasta"));
    if !src.exists() || !want.exists() {
        eprintln!("skip: {stem} reference FASTA absent (run validation/fetch_reference_data.sh)");
        return;
    }
    let out = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(format!("{stem}.revCat.fasta"));
    let opts = DecoyOptions {
        line_sep,
        output: Output::Concatenated,
        ..Default::default()
    };
    let n = write_decoy_database(&src, &out, &opts).unwrap();

    let got = std::fs::read(&out).unwrap();
    let expected = std::fs::read(&want).unwrap();
    assert_eq!(
        got.len(),
        expected.len(),
        "{stem}: length differs ({} vs {} bytes) after writing {n} decoys",
        got.len(),
        expected.len()
    );
    if got != expected {
        let at = got
            .iter()
            .zip(&expected)
            .position(|(a, b)| a != b)
            .unwrap_or(0);
        let lo = at.saturating_sub(60);
        panic!(
            "{stem}: first byte difference at offset {at}\n  got: {:?}\n want: {:?}",
            String::from_utf8_lossy(&got[lo..(at + 60).min(got.len())]),
            String::from_utf8_lossy(&expected[lo..(at + 60).min(expected.len())]),
        );
    }
    eprintln!(
        "{stem}: {n} decoys, {} bytes byte-identical to MS-GF+",
        got.len()
    );
}

#[test]
fn tryp_pig_bov_is_byte_identical() {
    // Windows-generated reference: CRLF throughout.
    check("Tryp_Pig_Bov", LineSep::Crlf);
}

#[test]
fn iprg2013_human_is_byte_identical() {
    // Linux-generated reference: LF, with the input's CRLF normalised away. Needs `--full` data.
    check("iprg2013_human", LineSep::Lf);
}

/// The load-time sanity gates MS-GF+ applies must accept a real concatenated database.
#[test]
fn reference_database_passes_validation() {
    let path = repo("validation/data/fasta/Tryp_Pig_Bov.revCat.fasta");
    if !path.exists() {
        eprintln!("skip: reference FASTA absent");
        return;
    }
    let db = msgf_db::fasta::ProteinDb::read(&path, msgf_db::fasta::DEFAULT_DECOY_PREFIX).unwrap();
    let accessions: Vec<String> = db.proteins.iter().map(|p| p.name.clone()).collect();
    msgf_db::decoy::validate_concatenated(&accessions, "XXX").expect("reference DB must validate");
    assert_eq!(
        db.n_decoys() * 2,
        db.proteins.len(),
        "half the DB should be decoys"
    );
}
