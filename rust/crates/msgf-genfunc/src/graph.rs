//! De novo graph builder: turns a `ScoredSpectrum` into the (node, edge) lists the DP consumes.
//!
//! Node score at nominal mass `m` = `round(node_score(m, prefix) + node_score(pep − m, suffix))`;
//! edge score from `ScoredSpectrum`'s validated edge scoring; each edge weighted by the amino-acid
//! background probability (uniform 1/20 in MS-GF+). The full builder — reverse graph for C-terminal
//! enzymes, sink tolerance, and source/neighboring cleavage credits — is layered on next and
//! validated against MS-GF+'s own DeNovoScore / SpecEValue oracle.

use crate::{Edge, Node};
use msgf_chem::{residue_mass, scaling};
use msgf_scorer::scored_spectrum::ScoredSpectrum;

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
pub fn build_reverse_graph(
    scored: &ScoredSpectrum,
    complement_mass: i32,
    sinks: &[i32],
    aa: &[Aa],
    peptide_credit: i32,
    peptide_penalty: i32,
) -> (Vec<Node>, Vec<usize>) {
    let graph_max = sinks
        .iter()
        .copied()
        .max()
        .unwrap_or(complement_mass)
        .max(complement_mass);
    let mut nodes: Vec<Node> = vec![Node::default(); graph_max as usize + 1];
    let sink_set: std::collections::HashSet<i32> = sinks.iter().copied().collect();
    let max_n = graph_max + 1;

    for m in 1..graph_max {
        if sink_set.contains(&m) || m >= complement_mass {
            continue; // sinks and out-of-range nodes score 0
        }
        nodes[m as usize].node_score = msgf_chem::round_half_up(
            scored.node_score(complement_mass - m, true) + scored.node_score(m, false),
        );
    }
    for m in 1..=graph_max {
        for a in aa {
            let prev = m - a.nominal;
            if prev < 0 {
                continue;
            }
            let mut es = scored.edge_score(m, prev, a.accurate_mass, max_n);
            if prev == 0 {
                es += if a.residue == b'K' || a.residue == b'R' {
                    peptide_credit
                } else {
                    peptide_penalty
                };
            }
            nodes[m as usize].edges.push(Edge {
                prev: prev as usize,
                edge_score: es,
                aa_prob: a.prob,
            });
        }
    }
    (nodes, sinks.iter().map(|&s| s as usize).collect())
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
