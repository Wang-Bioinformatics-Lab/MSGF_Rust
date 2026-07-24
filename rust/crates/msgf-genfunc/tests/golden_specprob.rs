//! Validates the generating function (DeNovoScore + SpecEValue p-value) against MS-GF+'s own
//! `f13_specprob.golden.json` — built on the **HighRes** model (the one the F13 search used), with
//! DB-composition amino-acid probabilities and the isotope-error sink range. Reproduces the search:
//! re-preprocess F13 raw peaks (HighRes model) → scored spectrum → reverse graph → DP → DeNovoScore
//! and spectral probability. Skipped if goldens/model/data are absent.

use msgf_genfunc::graph::{build_reverse_graph, standard_aa_nominal, Aa};
use msgf_genfunc::{compute, merge_group, Cleavage};
use msgf_io::MgfReader;
use msgf_scorer::preprocess::preprocess;
use msgf_scorer::scored_spectrum::ScoredSpectrum;
use serde_json::Value;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

fn repo(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

/// iPRG-2013 human FASTA composition (`DBScanner.setAminoAcidProbabilities`): residue → probability.
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
    let model = msgf_scorer::read_param_file(&param).unwrap();

    // precursor (charge, neutral mass) per scan — model-independent, from the scored-spectrum golden
    let mut prec: HashMap<i64, (i32, f32)> = HashMap::new();
    for s in ss["spectra"].as_array().unwrap() {
        prec.insert(
            s["scan"].as_i64().unwrap(),
            (
                s["charge"].as_i64().unwrap() as i32,
                s["precursor_mass"].as_f64().unwrap() as f32,
            ),
        );
    }
    // raw peaks by scan
    let mut raw: HashMap<i64, Vec<(f32, f32)>> = HashMap::new();
    for s in MgfReader::new(BufReader::new(File::open(&mgf).unwrap())) {
        let s = s.unwrap();
        if let Some(scan) = s.scan.as_deref().and_then(|x| x.parse::<i64>().ok()) {
            raw.insert(
                scan,
                s.peaks
                    .iter()
                    .map(|p| (p.mz as f32, p.intensity as f32))
                    .collect(),
            );
        }
    }

    let probs = iprg_probs();
    let mut aa: Vec<Aa> = standard_aa_nominal()
        .into_iter()
        .map(|(r, n)| Aa {
            residue: r,
            nominal: n,
            accurate_mass: msgf_chem::residue_mass(r).unwrap() as f32,
            prob: probs[&r],
        })
        .collect();
    // variable oxidation on M (iprg-2013_Mods.txt): an extra amino acid at M's DB probability
    let m_ox = msgf_chem::residue_mass(b'M').unwrap() + 15.994915;
    aa.push(Aa {
        residue: b'M',
        nominal: msgf_chem::scaling::nominal_bin(m_ox as f32),
        accurate_mass: m_ox as f32,
        prob: probs[&b'M'],
    });
    let prob_cleavage = probs[&b'K'] + probs[&b'R'];

    let (mut denovo_ok, mut spec_ok, mut total) = (0, 0, 0);
    let mut worst = String::new();
    for e in sp["spectra"].as_array().unwrap() {
        let scan = e["scan"].as_i64().unwrap();
        let (Some(&(charge, parent_mass)), Some(raw_peaks)) = (prec.get(&scan), raw.get(&scan))
        else {
            continue;
        };
        let complement = e["peptide_mass_nominal"].as_i64().unwrap() as i32;
        let range = e["mass_index_range"].as_array().unwrap();
        let sinks: Vec<i32> =
            (range[0].as_i64().unwrap() as i32..=range[1].as_i64().unwrap() as i32).collect();
        let g_raw = e["raw_score"].as_i64().unwrap() as i32;
        let g_denovo = e["denovo_score"].as_i64().unwrap() as i32;
        let g_spec = e["spec_prob"].as_f64().unwrap();
        total += 1;

        let _ = complement; // isotope range is handled per-mass below
        let peaks = preprocess(&model, charge, parent_mass, raw_peaks);
        let scored = ScoredSpectrum::from_ranked_peaks(&model, charge, parent_mass, peaks);
        let cleave = Cleavage {
            credit: 2,
            penalty: -11,
            prob_cleavage_sites: prob_cleavage,
        };
        // GeneratingFunctionGroup: one graph per candidate peptide mass (its own complement), merged.
        // Tables AND edges are candidate-independent, so build both once for the largest candidate
        // and only recompute node scores per candidate — the shared path used in production.
        let max_p = sinks.iter().copied().max().unwrap_or(0);
        let tables = scored.tables(max_p);
        let (mut graph, _) = build_reverse_graph(&scored, &tables, max_p, &[max_p], &aa, 2, -11);
        let mut gfs = Vec::new();
        for &p in &sinks {
            graph.recompute_node_scores(&tables, p, &[p]);
            if let Some(gf) = compute(&graph, &[p as usize], Some(cleave)) {
                gfs.push(gf);
            }
        }
        let Some(gf) = merge_group(&gfs) else {
            continue;
        };

        // direct ScoreDist diff for the sampled spectra
        if let Some(sd) = e.get("score_dist_sample") {
            let gmin = sd["min_score"].as_i64().unwrap() as i32;
            let gprobs: Vec<f64> = sd["probs"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_f64().unwrap())
                .collect();
            let (mut worst_d, mut worst_s) = (0.0f64, 0);
            for (i, &gp) in gprobs.iter().enumerate() {
                let score = gmin + i as i32;
                let mp = if score >= gf.dist.min_score
                    && (score - gf.dist.min_score) < gf.dist.probs.len() as i32
                {
                    gf.dist.probs[(score - gf.dist.min_score) as usize]
                } else {
                    0.0
                };
                if (mp - gp).abs() > worst_d {
                    worst_d = (mp - gp).abs();
                    worst_s = score;
                }
            }
            eprintln!("  [dist scan {scan}] my[min {} len {}] vs golden[min {} len {}]; worst Δprob {worst_d:.3e} at score {worst_s}",
                gf.dist.min_score, gf.dist.probs.len(), gmin, gprobs.len());
        }

        let my_denovo = gf.max_score();
        let my_spec = gf.spectral_probability(g_raw);
        let dlog = if g_spec > 0.0 && my_spec > 0.0 {
            (my_spec / g_spec).log10()
        } else {
            f64::NAN
        };
        denovo_ok += (my_denovo == g_denovo) as i32;
        // SpecEValue matches to f64 accumulation noise (distributions agree to ~2e-8)
        spec_ok += (dlog.abs() <= 1e-4) as i32;
        if (my_denovo != g_denovo || dlog.abs() > 1e-4) && worst.len() < 900 {
            worst.push_str(&format!(
                "\n  scan {scan}: denovo {my_denovo} vs {g_denovo} | spec {my_spec:.4e} vs {g_spec:.4e} (Δlog {dlog:.3})"
            ));
        }
    }
    eprintln!("DeNovoScore exact: {denovo_ok}/{total}; SpecEValue exact: {spec_ok}/{total}{worst}");
    assert!(total >= 25, "expected ~30 matched spectra");
    assert_eq!(denovo_ok, total, "DeNovoScore must match MS-GF+ exactly");
    assert_eq!(
        spec_ok, total,
        "SpecEValue (p-value) must match MS-GF+ exactly"
    );
}
