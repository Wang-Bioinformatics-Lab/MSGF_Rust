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

/// Build the **reverse** de novo graph (C-terminal enzyme, e.g. trypsin): nodes indexed by suffix
/// nominal mass `0..=pep_nominal`, node 0 = source, `pep_nominal` = sink.
///
/// Node score at `m` = `round(node_score(pep−m, prefix) + node_score(m, suffix))`; source and sink
/// are 0 (`computeNodeScores`). Each edge `prev = m − aa` carries the validated
/// `ScoredSpectrum::edge_score` weighted by [`AA_PROB`]; the first edge (from the source — the
/// C-terminal residue) also gets the peptide cleavage credit/penalty for K/R.
pub fn build_reverse_graph(
    scored: &ScoredSpectrum,
    pep_nominal: i32,
    peptide_credit: i32,
    peptide_penalty: i32,
) -> (Vec<Node>, Vec<usize>) {
    let n = pep_nominal.max(0) as usize + 1;
    let mut nodes: Vec<Node> = vec![Node::default(); n];
    let aa = standard_aa_nominal();
    let max_n = pep_nominal + 1;

    for m in 1..pep_nominal {
        nodes[m as usize].node_score =
            (scored.node_score(pep_nominal - m, true) + scored.node_score(m, false)).round() as i32;
    }
    for m in 1..=pep_nominal {
        for &(residue, aa_nom) in &aa {
            let prev = m - aa_nom;
            if prev < 0 {
                continue;
            }
            let aa_real = residue_mass(residue).expect("standard residue") as f32;
            let mut es = scored.edge_score(m, prev, aa_real, max_n);
            if prev == 0 {
                es += if residue == b'K' || residue == b'R' {
                    peptide_credit
                } else {
                    peptide_penalty
                };
            }
            nodes[m as usize].edges.push(Edge {
                prev: prev as usize,
                edge_score: es,
                aa_prob: AA_PROB,
            });
        }
    }
    (nodes, vec![pep_nominal as usize])
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
