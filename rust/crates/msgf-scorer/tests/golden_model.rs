//! Validates the `.param` reader against `validation/golden/models/*.model.golden.json`, each
//! derived from MS-GF+'s own authoritative `writeParametersPlainText` dump.
//!
//! A successful `read_param_file` already proves alignment (the reader checks the trailing
//! `0x7FFFFFFF` sentinel). On top of that we compare identity, scalars, the charge histogram,
//! partitions, precursor/fragment offsets, and rank/error distributions (grand sums + counts +
//! exact first-partition samples). The `.param` data is gitignored; missing files are skipped.

use msgf_scorer::ScoringModel;
use serde_json::Value;
use std::collections::HashMap;
use std::path::PathBuf;

fn repo(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}
fn approx(a: f64, b: f64, tol: f64, ctx: &str) {
    assert!(
        (a - b).abs() <= tol,
        "{ctx}: {a} vs {b} (Δ={:.3e} > {tol:.0e})",
        (a - b).abs()
    );
}

#[test]
fn all_high_res_models_match() {
    let gdir = repo("validation/golden/models");
    let mut validated = 0;

    for entry in std::fs::read_dir(&gdir)
        .expect("golden/models must exist")
        .flatten()
    {
        let name = entry.file_name().to_string_lossy().into_owned();
        if !name.ends_with(".model.golden.json") {
            continue;
        }
        let g: Value =
            serde_json::from_str(&std::fs::read_to_string(entry.path()).unwrap()).unwrap();
        let param = repo("validation/data/models").join(g["file"].as_str().unwrap());
        if !param.exists() {
            eprintln!("skip {}: data absent (fetch_reference_data.sh)", g["file"]);
            continue;
        }
        // Ok => the trailing 0x7FFFFFFF sentinel matched => the whole stream stayed aligned.
        let m =
            msgf_scorer::read_param_file(&param).unwrap_or_else(|e| panic!("parse {name}: {e}"));
        validate(&m, &g, &name);
        validated += 1;
    }

    assert!(validated > 0, "no .param data present to validate against");
}

fn validate(m: &ScoringModel, g: &Value, tag: &str) {
    // identity + scalars
    assert_eq!(
        m.activation,
        g["activation"].as_str().unwrap(),
        "{tag} activation"
    );
    assert_eq!(
        m.instrument,
        g["instrument"].as_str().unwrap(),
        "{tag} instrument"
    );
    assert_eq!(
        m.enzyme.as_deref().unwrap_or(""),
        g["enzyme"].as_str().unwrap(),
        "{tag} enzyme"
    );
    assert_eq!(
        m.protocol_name(),
        g["protocol"].as_str().unwrap(),
        "{tag} protocol"
    );
    approx(
        m.mme.value,
        g["mme"]["value"].as_f64().unwrap(),
        1e-6,
        "mme",
    );
    assert_eq!(
        m.mme.unit == msgf_chem::Unit::Ppm,
        g["mme"]["ppm"].as_bool().unwrap(),
        "{tag} mme ppm"
    );
    assert_eq!(
        m.apply_deconvolution,
        g["apply_deconvolution"].as_bool().unwrap(),
        "{tag} deconv flag"
    );
    approx(
        m.deconvolution_error_tolerance as f64,
        g["deconvolution_error_tolerance"].as_f64().unwrap(),
        1e-6,
        "deconv tol",
    );
    assert_eq!(
        m.num_segments as i64,
        g["num_segments"].as_i64().unwrap(),
        "{tag} num_segments"
    );

    // charge histogram (order-independent)
    let hist: HashMap<i64, i64> = m
        .charge_histogram
        .iter()
        .map(|&(c, n)| (c as i64, n as i64))
        .collect();
    for row in g["charge_histogram"].as_array().unwrap() {
        let c = row[0].as_i64().unwrap();
        assert_eq!(
            hist.get(&c),
            Some(&row[1].as_i64().unwrap()),
            "{tag} charge_histogram[{c}]"
        );
    }
    assert_eq!(
        hist.len(),
        g["charge_histogram"].as_array().unwrap().len(),
        "{tag} hist size"
    );

    // partitions (same sorted order as MS-GF+; charge/seg exact, parent_mass within f32 slop)
    let gparts = g["partitions"].as_array().unwrap();
    assert_eq!(m.partitions.len(), gparts.len(), "{tag} n partitions");
    for (p, gp) in m.partitions.iter().zip(gparts) {
        assert_eq!(
            p.charge as i64,
            gp[0].as_i64().unwrap(),
            "{tag} partition charge"
        );
        assert_eq!(p.seg as i64, gp[1].as_i64().unwrap(), "{tag} partition seg");
        approx(
            p.parent_mass as f64,
            gp[2].as_f64().unwrap(),
            0.05,
            "partition parent_mass",
        );
    }

    // precursor offsets
    let gpoff = g["precursor_off"].as_array().unwrap();
    assert_eq!(m.precursor_off.len(), gpoff.len(), "{tag} n precursor_off");
    for (o, go) in m.precursor_off.iter().zip(gpoff) {
        assert_eq!(
            o.charge as i64,
            go[0].as_i64().unwrap(),
            "{tag} precursor charge"
        );
        assert_eq!(
            o.reduced_charge as i64,
            go[1].as_i64().unwrap(),
            "{tag} reduced charge"
        );
        approx(
            o.offset as f64,
            go[2].as_f64().unwrap(),
            1e-4,
            "precursor offset",
        );
        approx(
            o.frequency as f64,
            go[3].as_f64().unwrap(),
            1e-6,
            "precursor freq",
        );
    }

    // fragment offsets: grand aggregates + first-partition sample
    let frag_entries: usize = m.frag_off.iter().map(|b| b.len()).sum();
    let frag_freq_sum: f64 = m
        .frag_off
        .iter()
        .flatten()
        .map(|f| f.frequency as f64)
        .sum();
    let frag_off_sum: f64 = m.frag_off.iter().flatten().map(|f| f.offset as f64).sum();
    assert_eq!(
        frag_entries as i64,
        g["frag_off_total_entries"].as_i64().unwrap(),
        "{tag} frag entries"
    );
    approx(
        frag_freq_sum,
        g["frag_off_freq_sum"].as_f64().unwrap(),
        1e-2,
        "frag freq sum",
    );
    approx(
        frag_off_sum,
        g["frag_off_offset_sum"].as_f64().unwrap(),
        1e-2,
        "frag offset sum",
    );
    let gsample = g["frag_off_sample"].as_object().unwrap();
    let sample0: HashMap<&str, &msgf_scorer::FragOff> =
        m.frag_off[0].iter().map(|f| (f.name.as_str(), f)).collect();
    assert_eq!(sample0.len(), gsample.len(), "{tag} frag sample size");
    for (name, vals) in gsample {
        let f = sample0
            .get(name.as_str())
            .unwrap_or_else(|| panic!("{tag} missing frag ion {name}"));
        approx(
            f.frequency as f64,
            vals[0].as_f64().unwrap(),
            1e-6,
            &format!("frag {name} freq"),
        );
        approx(
            f.offset as f64,
            vals[1].as_f64().unwrap(),
            1e-4,
            &format!("frag {name} offset"),
        );
    }

    // rank distributions: grand aggregates + first-partition sample (exact arrays)
    assert_eq!(
        m.max_rank as i64,
        g["max_rank"].as_i64().unwrap(),
        "{tag} max_rank"
    );
    let rank_floats: usize = m
        .rank_dist
        .iter()
        .flat_map(|r| r.ions.iter())
        .map(|(_, v)| v.len())
        .sum();
    let rank_sum: f64 = m
        .rank_dist
        .iter()
        .flat_map(|r| r.ions.iter())
        .flat_map(|(_, v)| v.iter())
        .map(|&x| x as f64)
        .sum();
    assert_eq!(
        rank_floats as i64,
        g["rank_dist_total_floats"].as_i64().unwrap(),
        "{tag} rank floats"
    );
    approx(
        rank_sum,
        g["rank_dist_freq_sum"].as_f64().unwrap(),
        0.1,
        "rank freq sum",
    );
    let grank = g["rank_dist_sample"].as_object().unwrap();
    let rank0: HashMap<&str, &Vec<f32>> = m.rank_dist[0]
        .ions
        .iter()
        .map(|(n, v)| (n.as_str(), v))
        .collect();
    assert_eq!(rank0.len(), grank.len(), "{tag} rank sample ions");
    for (name, arr) in grank {
        let got = rank0
            .get(name.as_str())
            .unwrap_or_else(|| panic!("{tag} missing rank ion {name}"));
        let exp = arr.as_array().unwrap();
        assert_eq!(got.len(), exp.len(), "rank {name} len");
        for (k, e) in exp.iter().enumerate() {
            approx(
                got[k] as f64,
                e.as_f64().unwrap(),
                1e-6,
                &format!("rank {name}[{k}]"),
            );
        }
    }

    // error distributions
    assert_eq!(
        m.error_scaling_factor as i64,
        g["error_scaling_factor"].as_i64().unwrap(),
        "{tag} esf"
    );
    let sig_sum: f64 = m
        .error_dist
        .iter()
        .flat_map(|e| e.signal.iter())
        .map(|&x| x as f64)
        .sum();
    let noise_sum: f64 = m
        .error_dist
        .iter()
        .flat_map(|e| e.noise.iter())
        .map(|&x| x as f64)
        .sum();
    let ionex_sum: f64 = m
        .error_dist
        .iter()
        .flat_map(|e| e.ion_existence.iter())
        .map(|&x| x as f64)
        .sum();
    approx(
        sig_sum,
        g["error_signal_sum"].as_f64().unwrap(),
        1e-2,
        "error signal sum",
    );
    approx(
        noise_sum,
        g["error_noise_sum"].as_f64().unwrap(),
        1e-2,
        "error noise sum",
    );
    approx(
        ionex_sum,
        g["ion_existence_sum"].as_f64().unwrap(),
        1e-2,
        "ion existence sum",
    );
    let es = &g["error_sample"];
    let e0 = &m.error_dist[0];
    for (k, e) in es["signal"].as_array().unwrap().iter().enumerate() {
        approx(
            e0.signal[k] as f64,
            e.as_f64().unwrap(),
            1e-6,
            &format!("err signal[{k}]"),
        );
    }
    for (k, e) in es["ion_existence"].as_array().unwrap().iter().enumerate() {
        approx(
            e0.ion_existence[k] as f64,
            e.as_f64().unwrap(),
            1e-6,
            &format!("ion_existence[{k}]"),
        );
    }

    eprintln!(
        "ok {}: {} partitions, {rank_floats} rank floats, aligned to terminator",
        g["file"],
        m.partitions.len()
    );
}
