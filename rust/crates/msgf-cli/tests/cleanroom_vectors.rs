//! Byte-identity of the `msgf` binary against the clean-room test vectors.
//!
//! The scorer and generating function were replaced on 2026-09-30 by a clean-room implementation
//! (`docs/cleanroom/`). The test-vector package written for that process holds the recorded output
//! of the previous MSGF_Rust release (`0fb0738`) on 3,000 HeLa spectra and 15 synthetic edge cases:
//! seven `msgf rescore` configurations and one `msgf search`. This test reruns the current binary
//! on the same inputs and requires **byte-identical** stdout, and identical skip lines on stderr.
//!
//! The package (~20 MB) is not vendored. Point `MSGF_CLEANROOM_VECTORS` at its `test_vectors/`
//! directory to run the rescore sets; also set `MSGF_CLEANROOM_FASTA` to UniProt `UP000005640`
//! (human, SHA-256 `2329a517…84d5a0`) to run the search set. Without them the test skips.

use std::path::{Path, PathBuf};
use std::process::Command;

fn vectors() -> Option<PathBuf> {
    let dir = PathBuf::from(std::env::var_os("MSGF_CLEANROOM_VECTORS")?);
    dir.join("rescore").is_dir().then_some(dir)
}

fn msgf(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_msgf"))
        .args(args)
        .output()
        .expect("run msgf")
}

fn path(p: &Path) -> &str {
    p.to_str().expect("utf-8 path")
}

#[test]
fn rescore_sets_are_byte_identical() {
    let Some(v) = vectors() else {
        eprintln!("skip: MSGF_CLEANROOM_VECTORS not set");
        return;
    };
    let probs = v.join("inputs/human_aa_probs.tsv");
    let sets: [(&str, &str, &str, &[&str]); 7] = [
        ("hela", "hela_r01_sub.mgf", "hela_uniform_ti0_1", &[]),
        (
            "hela",
            "hela_r01_sub.mgf",
            "hela_uniform_ti0_0",
            &["--ti", "0,0"],
        ),
        (
            "hela",
            "hela_r01_sub.mgf",
            "hela_uniform_ti-1_2",
            &["--ti", "-1,2"],
        ),
        (
            "hela",
            "hela_r01_sub.mgf",
            "hela_composition_oxm_ti0_1",
            &["COMPOSITION", "--ox-m"],
        ),
        ("synthetic", "synthetic.mgf", "synthetic_uniform_ti0_1", &[]),
        (
            "synthetic",
            "synthetic.mgf",
            "synthetic_uniform_ti-1_2",
            &["--ti", "-1,2"],
        ),
        (
            "synthetic",
            "synthetic.mgf",
            "synthetic_composition_oxm_ti0_1",
            &["COMPOSITION", "--ox-m"],
        ),
    ];
    let tmp = PathBuf::from(env!("CARGO_TARGET_TMPDIR"));
    for (psm_set, mgf, name, extra) in sets {
        let spectra = v.join("spectra").join(mgf);
        let psms = v.join(format!("inputs/psms_{psm_set}.tsv"));
        let out = tmp.join(format!("{name}.tsv"));
        let mut args = vec![
            "rescore",
            "-s",
            path(&spectra),
            "-i",
            path(&psms),
            "-o",
            path(&out),
        ];
        for &a in extra {
            if a == "COMPOSITION" {
                args.extend(["--aa-probs", path(&probs)]);
            } else {
                args.push(a);
            }
        }
        let run = msgf(&args);
        assert!(run.status.success(), "{name}: rescore failed");
        let got = std::fs::read(&out).unwrap();
        let want = std::fs::read(v.join(format!("rescore/{name}.tsv"))).unwrap();
        assert!(
            got == want,
            "{name}: output differs from the recorded vector"
        );
        let want_err = std::fs::read_to_string(v.join(format!("rescore/{name}.stderr"))).unwrap();
        assert_eq!(
            String::from_utf8_lossy(&run.stderr),
            want_err,
            "{name}: stderr differs"
        );
        eprintln!("ok: {name} byte-identical ({} bytes)", got.len());
    }
}

#[test]
fn search_set_is_byte_identical() {
    let (Some(v), Some(fasta)) = (vectors(), std::env::var_os("MSGF_CLEANROOM_FASTA")) else {
        eprintln!("skip: MSGF_CLEANROOM_VECTORS / MSGF_CLEANROOM_FASTA not set");
        return;
    };
    let fasta = PathBuf::from(fasta);
    let spectra = v.join("spectra/hela_r01_sub.mgf");
    let out = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("search_n5.tsv");
    let run = msgf(&[
        "search",
        "-s",
        path(&spectra),
        "-d",
        path(&fasta),
        "--tda",
        "--fixed-mod",
        "C+57.021464",
        "--var-mod",
        "M+15.994915",
        "-t",
        "10ppm",
        "--ti",
        "0,1",
        "-c",
        "2",
        "-n",
        "5",
        "--threads",
        "16",
        "-o",
        path(&out),
    ]);
    assert!(
        run.status.success(),
        "search failed: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let got = std::fs::read(&out).unwrap();
    let want = std::fs::read(v.join("search_hela_r01/search_n5.tsv")).unwrap();
    assert!(
        got == want,
        "search output differs from the recorded vector"
    );
    eprintln!("ok: search_n5 byte-identical ({} bytes)", got.len());
}
