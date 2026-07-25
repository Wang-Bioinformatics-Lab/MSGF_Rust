//! Pins the compact graph representation (`edge_aa` + `aa_prob` instead of a per-edge `f64`
//! `edge_prob`) and the reusable-buffer builder **bit-exact** against the previous builder.
//!
//! `reference_build_reverse_graph` below is a verbatim copy of the pre-change
//! `graph::build_reverse_graph` — five freshly allocated `Vec`s, a two-pass edge count, an
//! unconditional `edge_score_with` per edge, and a per-edge `f64` probability. The test asserts:
//!
//!  1. every CSR array is element-identical, and the probability the DP would read for edge `e`
//!     has the identical `f64::to_bits` pattern in both;
//!  2. the DP's `ScoreDist` over both graphs is identical **by `to_bits`** — not by epsilon — for
//!     every score cell, over the real F13 isotope-error mass group;
//!  3. a reused (dirty) `Graph` produces exactly the same arrays as a fresh one, so buffer reuse
//!     cannot leak stale state.
//!
//! Skipped when `validation/data/` is absent (the data-absence contract).

use msgf_chem::{mass, scaling};
use msgf_genfunc::graph::{
    build_reverse_graph, build_reverse_graph_into, standard_aa_nominal, Aa, PeptideCleavage,
};
use msgf_genfunc::{compute, Cleavage, Graph};
use msgf_io::MgfReader;
use msgf_scorer::preprocess::preprocess;
use msgf_scorer::scored_spectrum::{ScoredSpectrum, SpectrumTables};
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

fn repo(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

/// The pre-change builder's output: `(node_score, edge_start, edge_prev, edge_score, edge_prob)`.
type RefGraph = (Vec<i32>, Vec<u32>, Vec<u32>, Vec<i32>, Vec<f64>);

/// The pre-change builder, copied verbatim. Returns the CSR arrays with the per-edge `f64`
/// probability the DP used to read directly. **Do not "clean up"** — its value is being the old
/// code, so any divergence in the new builder shows up as a test failure rather than a silent
/// fidelity change.
fn reference_build_reverse_graph(
    scored: &ScoredSpectrum,
    tables: &SpectrumTables,
    complement_mass: i32,
    sinks: &[i32],
    aa: &[Aa],
    cleavage: PeptideCleavage<'_>,
) -> RefGraph {
    let graph_max = sinks
        .iter()
        .copied()
        .max()
        .unwrap_or(complement_mass)
        .max(complement_mass);
    let n = graph_max as usize + 1;
    let is_sink = |m: i32| sinks.contains(&m);

    let mut node_score = vec![0i32; n];
    for m in 1..graph_max {
        if is_sink(m) || m >= complement_mass {
            continue;
        }
        node_score[m as usize] = msgf_chem::round_half_up(
            tables.prefix[(complement_mass - m) as usize] + tables.suffix[m as usize],
        );
    }

    let mut edge_start = vec![0u32; n + 1];
    for m in 1..=graph_max {
        let mut c = 0u32;
        for a in aa {
            if m - a.nominal >= 0 {
                c += 1;
            }
        }
        edge_start[m as usize] = c;
    }
    let mut acc = 0u32;
    for slot in edge_start.iter_mut() {
        let c = *slot;
        *slot = acc;
        acc += c;
    }
    let total = acc as usize;

    let node_mass = &tables.node_mass;
    let mut edge_prev = vec![0u32; total];
    let mut edge_score = vec![0i32; total];
    let mut edge_prob = vec![0f64; total];
    for m in 1..=graph_max {
        let mut pos = edge_start[m as usize] as usize;
        for a in aa {
            let prev = m - a.nominal;
            if prev < 0 {
                continue;
            }
            let mut es = scored.edge_score_with(
                node_mass[m as usize],
                node_mass[prev as usize],
                a.accurate_mass,
            );
            if prev == 0 {
                es += cleavage.score(a.residue);
            }
            edge_prev[pos] = prev as u32;
            edge_score[pos] = es;
            edge_prob[pos] = a.prob;
            pos += 1;
        }
    }
    (node_score, edge_start, edge_prev, edge_score, edge_prob)
}

/// The reference arrays as a `Graph`, so the DP can be run over them. `from_adj` interns the
/// probabilities by exact bit pattern, so the values the DP reads are the reference's own `f64`s.
fn reference_graph(
    node_score: &[i32],
    edge_start: &[u32],
    edge_prev: &[u32],
    edge_score: &[i32],
    edge_prob: &[f64],
) -> Graph {
    let adj: Vec<msgf_genfunc::AdjNode> = (0..node_score.len())
        .map(|i| {
            let (e0, e1) = (edge_start[i] as usize, edge_start[i + 1] as usize);
            (
                node_score[i],
                (e0..e1)
                    .map(|e| (edge_prev[e] as usize, edge_score[e], edge_prob[e]))
                    .collect(),
            )
        })
        .collect();
    Graph::from_adj(&adj)
}

#[test]
fn compact_graph_is_bit_exact_against_the_previous_builder() {
    let mgf = repo("validation/data/spectra/F13.mgf");
    let param = repo("validation/data/models/HCD_HighRes_Tryp.param");
    if !mgf.exists() || !param.exists() {
        eprintln!("skip: validation/data absent");
        return;
    }
    let model = msgf_scorer::read_param_file(&param).expect("HCD_HighRes_Tryp.param");
    let aa: Vec<Aa> = standard_aa_nominal()
        .into_iter()
        .map(|(residue, nominal)| Aa {
            residue,
            nominal,
            accurate_mass: msgf_chem::residue_mass(residue).unwrap() as f32,
            // Deliberately non-uniform: uniform probabilities would hide any mix-up of the
            // per-edge amino-acid index (I/L share nominal 113, K/Q share 128).
            prob: 0.03 + (residue as f64 % 7.0) * 0.011,
        })
        .collect();
    let cleave = Cleavage {
        credit: 2,
        penalty: -11,
        prob_cleavage_sites: 0.10,
    };

    // A dirty, reused graph: proves buffer reuse cannot leak stale state into the next spectrum.
    let mut reused = Graph::default();
    let mut checked = 0usize;

    let file = File::open(&mgf).expect("F13.mgf");
    for spec in MgfReader::new(BufReader::new(file)) {
        let spec = spec.expect("mgf record");
        let (Some(charge), Some(mz)) = (spec.charge, spec.precursor_mz) else {
            continue;
        };
        let parent_mass = mz as f32 * charge as f32 - charge as f32 * mass::PROTON as f32;
        let pep_nominal = scaling::nominal_bin(parent_mass - mass::WATER as f32);
        if !(200..=6000).contains(&pep_nominal) {
            continue;
        }
        let peaks: Vec<(f32, f32)> = spec
            .peaks
            .iter()
            .map(|p| (p.mz as f32, p.intensity as f32))
            .collect();
        let ranked = preprocess(&model, charge, parent_mass, &peaks);
        let scored = ScoredSpectrum::from_ranked_peaks(&model, charge, parent_mass, ranked);
        let tables = scored.tables(pep_nominal);

        let (ref_ns, ref_start, ref_prev, ref_score, ref_prob) = reference_build_reverse_graph(
            &scored,
            &tables,
            pep_nominal,
            &[pep_nominal],
            &aa,
            PeptideCleavage::TRYPSIN,
        );
        let (fresh, _) = build_reverse_graph(
            &scored,
            &tables,
            pep_nominal,
            &[pep_nominal],
            &aa,
            PeptideCleavage::TRYPSIN,
        );
        build_reverse_graph_into(
            &mut reused,
            &scored,
            &tables,
            pep_nominal,
            &[pep_nominal],
            &aa,
            PeptideCleavage::TRYPSIN,
        );

        for g in [&fresh, &reused] {
            assert_eq!(g.node_score, ref_ns, "node_score");
            assert_eq!(g.edge_start, ref_start, "edge_start");
            assert_eq!(g.edge_prev, ref_prev, "edge_prev");
            assert_eq!(g.edge_score, ref_score, "edge_score");
            assert_eq!(g.edge_aa.len(), ref_prob.len(), "edge count");
            for (e, &want) in ref_prob.iter().enumerate() {
                let got = g.aa_prob[g.edge_aa[e] as usize];
                assert_eq!(
                    got.to_bits(),
                    want.to_bits(),
                    "edge {e}: probability bit pattern differs ({got} vs {want})"
                );
            }
        }

        // End-to-end: the DP over the reference arrays vs. over the compact graph, cell by cell,
        // by bit pattern — over the real `-ti 0,1` isotope mass group.
        let mut ref_g = reference_graph(&ref_ns, &ref_start, &ref_prev, &ref_score, &ref_prob);
        let mut new_g = fresh;
        for p in (pep_nominal - 1..=pep_nominal).filter(|&p| p > 0) {
            ref_g.recompute_node_scores(&tables, p, &[p]);
            new_g.recompute_node_scores(&tables, p, &[p]);
            assert_eq!(ref_g.node_score, new_g.node_score, "recomputed node scores");
            let a = compute(&ref_g, &[p as usize], Some(cleave));
            let b = compute(&new_g, &[p as usize], Some(cleave));
            match (a, b) {
                (None, None) => {}
                (Some(a), Some(b)) => {
                    assert_eq!(a.dist.min_score, b.dist.min_score, "ScoreDist min_score");
                    assert_eq!(a.dist.probs.len(), b.dist.probs.len(), "ScoreDist width");
                    for (i, (x, y)) in a.dist.probs.iter().zip(&b.dist.probs).enumerate() {
                        assert_eq!(
                            x.to_bits(),
                            y.to_bits(),
                            "score {} ({x:e} vs {y:e}) is not bit-identical",
                            a.dist.min_score + i as i32
                        );
                    }
                }
                _ => panic!("one path produced a distribution and the other did not"),
            }
        }

        checked += 1;
        if checked == 40 {
            break;
        }
    }
    assert!(checked > 0, "no F13 spectrum was checked");
    eprintln!("graph bit-exactness: {checked}/{checked} spectra identical (to_bits)");
}

/// The reuse path must not depend on the incoming buffer sizes: a large graph followed by a small
/// one (and vice versa) has to leave exactly the arrays a fresh build would.
#[test]
fn reuse_survives_shrinking_and_growing() {
    let mgf = repo("validation/data/spectra/F13.mgf");
    let param = repo("validation/data/models/HCD_HighRes_Tryp.param");
    if !mgf.exists() || !param.exists() {
        eprintln!("skip: validation/data absent");
        return;
    }
    let model = msgf_scorer::read_param_file(&param).expect("HCD_HighRes_Tryp.param");
    let aa: Vec<Aa> = standard_aa_nominal()
        .into_iter()
        .map(|(residue, nominal)| Aa {
            residue,
            nominal,
            accurate_mass: msgf_chem::residue_mass(residue).unwrap() as f32,
            prob: 0.05,
        })
        .collect();

    let file = File::open(&mgf).expect("F13.mgf");
    let spec = MgfReader::new(BufReader::new(file))
        .filter_map(|s| s.ok())
        .find(|s| s.charge.is_some() && s.precursor_mz.is_some())
        .expect("one usable F13 spectrum");
    let (charge, mz) = (spec.charge.unwrap(), spec.precursor_mz.unwrap());
    let parent_mass = mz as f32 * charge as f32 - charge as f32 * mass::PROTON as f32;
    let peaks: Vec<(f32, f32)> = spec
        .peaks
        .iter()
        .map(|p| (p.mz as f32, p.intensity as f32))
        .collect();
    let ranked = preprocess(&model, charge, parent_mass, &peaks);
    let scored = ScoredSpectrum::from_ranked_peaks(&model, charge, parent_mass, ranked);
    let tables = scored.tables(3000);

    let mut reused = Graph::default();
    // Big, then small, then big again — every ordering the reuse path can see.
    for &m in &[2500i32, 400, 1200, 2500, 300] {
        build_reverse_graph_into(
            &mut reused,
            &scored,
            &tables,
            m,
            &[m],
            &aa,
            PeptideCleavage::TRYPSIN,
        );
        let (fresh, _) =
            build_reverse_graph(&scored, &tables, m, &[m], &aa, PeptideCleavage::TRYPSIN);
        assert_eq!(reused.node_score, fresh.node_score, "node_score at {m}");
        assert_eq!(reused.edge_start, fresh.edge_start, "edge_start at {m}");
        assert_eq!(reused.edge_prev, fresh.edge_prev, "edge_prev at {m}");
        assert_eq!(reused.edge_score, fresh.edge_score, "edge_score at {m}");
        assert_eq!(reused.edge_aa, fresh.edge_aa, "edge_aa at {m}");
        assert_eq!(reused.aa_prob.len(), fresh.aa_prob.len(), "aa_prob at {m}");
        for (x, y) in reused.aa_prob.iter().zip(&fresh.aa_prob) {
            assert_eq!(x.to_bits(), y.to_bits(), "aa_prob bits at {m}");
        }
    }
}
