//! Validates the MGF reader against `validation/golden/spectra/*.json`: per-spectrum charge,
//! peak count and canonical peak-list hash, plus whole-file spectrum count, total peaks, and a
//! rolling hash over every spectrum. mzML goldens are skipped until the mzML reader lands.
//! Data files live under the gitignored `validation/data/`; missing files are skipped, not failed.

use msgf_io::MgfReader;
use serde_json::Value;
use sha1::{Digest, Sha1};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

fn repo_subdir(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

fn sha1_hex(bytes: &[u8]) -> String {
    let mut h = Sha1::new();
    h.update(bytes);
    h.finalize().iter().map(|b| format!("{b:02x}")).collect()
}

#[test]
fn mgf_golden_matches() {
    let gdir = repo_subdir("validation/golden/spectra");
    let mut checked_files = 0;

    for entry in std::fs::read_dir(&gdir)
        .expect("golden/spectra must exist")
        .flatten()
    {
        let g: Value =
            serde_json::from_str(&std::fs::read_to_string(entry.path()).unwrap()).unwrap();
        if g["format"] != "mgf" {
            continue; // mzML handled once the mzML reader exists
        }
        let file = g["file"].as_str().unwrap();
        let data = repo_subdir("validation/data/spectra").join(file);
        if !data.exists() {
            eprintln!("skip {file}: data absent (fetch_reference_data.sh)");
            continue;
        }

        // stored per-spectrum expectations, keyed by index
        let mut stored: HashMap<u64, &Value> = HashMap::new();
        for s in g["spectra"].as_array().unwrap() {
            stored.insert(s["index"].as_u64().unwrap(), s);
        }

        let reader = MgfReader::new(BufReader::new(File::open(&data).unwrap()));
        let mut count: u64 = 0;
        let mut total_peaks: u64 = 0;
        let mut roll = Sha1::new();

        for spec in reader {
            let spec = spec.unwrap();
            let peak_hash = sha1_hex(spec.canonical_peak_string().as_bytes());
            count += 1;
            total_peaks += spec.n_peaks() as u64;
            roll.update(peak_hash.as_bytes());

            if let Some(exp) = stored.get(&(spec.index as u64)) {
                assert_eq!(
                    spec.n_peaks() as u64,
                    exp["n_peaks"].as_u64().unwrap(),
                    "{file}#{} n_peaks",
                    spec.index
                );
                let exp_charge = exp["charge"].as_i64().map(|v| v as i32);
                assert_eq!(spec.charge, exp_charge, "{file}#{} charge", spec.index);
                assert_eq!(
                    peak_hash,
                    exp["peaks_sha1"].as_str().unwrap(),
                    "{file}#{} peaks_sha1",
                    spec.index
                );
            }
        }

        assert_eq!(count, g["n_spectra"].as_u64().unwrap(), "{file} n_spectra");
        assert_eq!(
            total_peaks,
            g["total_peaks"].as_u64().unwrap(),
            "{file} total_peaks"
        );
        let roll_hex: String = roll.finalize().iter().map(|b| format!("{b:02x}")).collect();
        assert_eq!(
            roll_hex,
            g["rolling_peak_sha1"].as_str().unwrap(),
            "{file} rolling hash"
        );
        eprintln!("ok {file}: {count} spectra, {total_peaks} peaks");
        checked_files += 1;
    }

    assert!(
        checked_files > 0,
        "no MGF golden data present to validate against"
    );
}
