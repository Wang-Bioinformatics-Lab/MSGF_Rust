//! De novo graph builder: turns a `ScoredSpectrum` into the (node, edge) lists the DP consumes.
//!
//! Node score at nominal mass `m` = `round(node_score(m, prefix) + node_score(pep − m, suffix))`;
//! edge score from `ScoredSpectrum`'s validated edge scoring; each edge weighted by the amino-acid
//! background probability (uniform 1/20 in MS-GF+). The full builder — reverse graph for C-terminal
//! enzymes, sink tolerance, and source/neighboring cleavage credits — is layered on next and
//! validated against MS-GF+'s own DeNovoScore / SpecEValue oracle.

use crate::Graph;
use msgf_chem::{residue_mass, scaling};
use msgf_scorer::scored_spectrum::{ScoredSpectrum, SpectrumTables};

/// Uniform amino-acid background probability in MS-GF+ (`AminoAcid.probability = 0.05`).
pub const AA_PROB: f64 = 0.05;

/// The 20 standard amino acids as `(residue, nominal residue mass)` for building graph edges.
pub fn standard_aa_nominal() -> Vec<(u8, i32)> {
    b"GASPVTCLINDQKEMHFRYW"
        .iter()
        .map(|&aa| {
            (
                aa,
                scaling::nominal_bin(residue_mass(aa).expect("standard residue") as f32),
            )
        })
        .collect()
}

/// An amino acid for graph edges: residue (unmodified, for cleavage checks), its nominal and
/// accurate residue masses (including any modification), and its (DB-composition) probability.
#[derive(Debug, Clone, Copy)]
pub struct Aa {
    pub residue: u8,
    pub nominal: i32,
    pub accurate_mass: f32,
    pub prob: f64,
}

/// The 20 standard amino acids with uniform probability (`AA_PROB`), for de novo without a DB.
pub fn standard_aa() -> Vec<Aa> {
    standard_aa_nominal()
        .into_iter()
        .map(|(residue, nominal)| Aa {
            residue,
            nominal,
            accurate_mass: residue_mass(residue).expect("standard residue") as f32,
            prob: AA_PROB,
        })
        .collect()
}

/// Build the **reverse** de novo graph (C-terminal enzyme, e.g. trypsin): nodes indexed by suffix
/// nominal mass, node 0 = source, the masses in `sinks` are sinks (the precursor + isotope-error
/// mass range). `complement_mass` is the peptide mass used for the node-score complement.
///
/// Node score at intermediate `m` = `round(node_score(complement−m, prefix) + node_score(m, suffix))`;
/// source and sinks are 0. Each edge `prev = m − aa` carries the validated `ScoredSpectrum::edge_score`
/// weighted by the amino acid's probability; the source edge (C-terminal residue) also gets the
/// peptide cleavage credit/penalty for K/R.
///
/// `tables` holds the per-node quantities (main-ion mass, prefix/suffix node scores) precomputed
/// once per spectrum via [`ScoredSpectrum::tables`] — none depend on `complement_mass`, so the
/// isotope-error candidate graphs reuse one instance instead of redoing the peak lookups. It must
/// cover `0..=max(sinks, complement_mass)`.
pub fn build_reverse_graph(
    scored: &ScoredSpectrum,
    tables: &SpectrumTables,
    complement_mass: i32,
    sinks: &[i32],
    aa: &[Aa],
    peptide_credit: i32,
    peptide_penalty: i32,
) -> (Graph, Vec<usize>) {
    let graph_max = sinks
        .iter()
        .copied()
        .max()
        .unwrap_or(complement_mass)
        .max(complement_mass);
    let n = graph_max as usize + 1;
    // `sinks` is tiny (the isotope range), so a linear check beats allocating a HashSet.
    let is_sink = |m: i32| sinks.contains(&m);

    let mut node_score = vec![0i32; n];
    for m in 1..graph_max {
        if is_sink(m) || m >= complement_mass {
            continue; // sinks and out-of-range nodes score 0
        }
        node_score[m as usize] = msgf_chem::round_half_up(
            tables.prefix[(complement_mass - m) as usize] + tables.suffix[m as usize],
        );
    }

    // CSR edges, built in two passes with no per-node allocation. Pass 1: count incoming edges per
    // node (one per amino acid whose predecessor mass is >= 0) and prefix-sum into `edge_start`.
    let mut edge_start = vec![0u32; n + 1];
    for m in 1..=graph_max {
        let mut c = 0u32;
        for a in aa {
            if m - a.nominal >= 0 {
                c += 1;
            }
        }
        edge_start[m as usize] = c; // stash the per-node count, converted to offsets below
    }
    let mut acc = 0u32;
    for slot in edge_start.iter_mut() {
        let c = *slot;
        *slot = acc;
        acc += c;
    }
    let total = acc as usize;

    // Pass 2: fill edges in amino-acid order (the DP's summation order — preserved for bit-exact
    // reproduction). Edge scores are **candidate-independent** — they depend only on the shared
    // node masses and the amino acid, never on `complement_mass` or which node is the sink — so the
    // edge arrays are built once and reused across the isotope-error candidates (see
    // [`Graph::recompute_node_scores`]). The sink's edges (errorScore 0, `setBackwardEdgesFromSink`)
    // are zeroed in the DP instead of here, since the sink differs per candidate.
    //
    // Node masses come from the shared `tables` (one peak lookup per node per spectrum, not one per
    // incident edge per graph). `m, prev` are always in-range, so the `edge_score` bounds guard
    // never fires; `edge_score_with` on the resolved masses is equivalent.
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
            let mut es =
                scored.edge_score_with(node_mass[m as usize], node_mass[prev as usize], a.accurate_mass);
            if prev == 0 {
                es += if a.residue == b'K' || a.residue == b'R' {
                    peptide_credit
                } else {
                    peptide_penalty
                };
            }
            edge_prev[pos] = prev as u32;
            edge_score[pos] = es;
            edge_prob[pos] = a.prob;
            pos += 1;
        }
    }

    let graph = Graph {
        node_score,
        edge_start,
        edge_prev,
        edge_score,
        edge_prob,
    };
    (graph, sinks.iter().map(|&s| s as usize).collect())
}

impl Graph {
    /// Recompute the node scores for a different candidate peptide mass (`complement_mass` / `sinks`)
    /// **without rebuilding the edges** — the edges are candidate-independent (see
    /// [`build_reverse_graph`]). This is how one built graph serves every isotope-error candidate:
    /// build once for the largest candidate, then call this per candidate before the DP. Node scores
    /// for masses above the candidate are zeroed; the DP only visits nodes up to the candidate's sink.
    pub fn recompute_node_scores(
        &mut self,
        tables: &SpectrumTables,
        complement_mass: i32,
        sinks: &[i32],
    ) {
        let graph_max = sinks
            .iter()
            .copied()
            .max()
            .unwrap_or(complement_mass)
            .max(complement_mass);
        let is_sink = |m: i32| sinks.contains(&m);
        for s in self.node_score.iter_mut() {
            *s = 0;
        }
        for m in 1..graph_max {
            if is_sink(m) || m >= complement_mass {
                continue;
            }
            self.node_score[m as usize] = msgf_chem::round_half_up(
                tables.prefix[(complement_mass - m) as usize] + tables.suffix[m as usize],
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn aa_nominal_masses() {
        let m: HashMap<u8, i32> = standard_aa_nominal().into_iter().collect();
        assert_eq!(m.len(), 20);
        assert_eq!(m[&b'G'], 57); // glycine
        assert_eq!(m[&b'A'], 71); // alanine
        assert_eq!(m[&b'W'], 186); // tryptophan
        assert_eq!(m[&b'K'], 128);
        assert_eq!(m[&b'R'], 156);
    }
}
