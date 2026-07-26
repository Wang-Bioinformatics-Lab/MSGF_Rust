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

/// Cleavage scoring for the peptide's own C-terminus, applied to the graph's source edge.
///
/// `cleave_at` are the enzyme's cleavage residues — [`PeptideCleavage::TRYPSIN`] (`KR`, +2/-11) is
/// the configuration validated bit-exact against MS-GF+. A residue in the set earns `credit`, any
/// other earns `penalty`; use [`PeptideCleavage::NONE`] to switch cleavage scoring off entirely,
/// which is what an unspecific enzyme needs.
#[derive(Debug, Clone, Copy)]
pub struct PeptideCleavage<'a> {
    pub cleave_at: &'a [u8],
    pub credit: i32,
    pub penalty: i32,
}

impl PeptideCleavage<'_> {
    /// Trypsin: cleaves after K/R, +2 credit and -11 penalty. The validated default.
    pub const TRYPSIN: PeptideCleavage<'static> = PeptideCleavage {
        cleave_at: b"KR",
        credit: 2,
        penalty: -11,
    };

    /// No cleavage scoring — every source edge scores 0.
    pub const NONE: PeptideCleavage<'static> = PeptideCleavage {
        cleave_at: &[],
        credit: 0,
        penalty: 0,
    };

    /// The credit or penalty this residue earns at the peptide's C-terminus.
    #[inline]
    pub fn score(&self, residue: u8) -> i32 {
        if self.cleave_at.contains(&residue) {
            self.credit
        } else {
            self.penalty
        }
    }
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
///
/// `cleavage` scores the source edge — the peptide's C-terminal residue. See [`PeptideCleavage`].
pub fn build_reverse_graph(
    scored: &ScoredSpectrum,
    tables: &SpectrumTables,
    complement_mass: i32,
    sinks: &[i32],
    aa: &[Aa],
    cleavage: PeptideCleavage<'_>,
) -> (Graph, Vec<usize>) {
    let mut graph = Graph::default();
    let sink_idx = build_reverse_graph_into(
        &mut graph,
        scored,
        tables,
        complement_mass,
        sinks,
        aa,
        cleavage,
    );
    (graph, sink_idx)
}

/// Fill the CSR offsets `edge_start[0..=graph_max + 1]` and return the total edge count.
///
/// The incoming-edge count of node `m` is a pure function of the node mass — the number of amino
/// acids with `nominal <= m` — so it saturates at `|aa|` for every node at or above the heaviest
/// residue (only ~186 of the ~1,500 nodes are below it). Count those few by scanning the alphabet
/// and fill the rest from the saturated count, rather than rescanning the alphabet for every node
/// and prefix-summing in a second pass.
///
/// Node 0 is the source and always has zero incoming edges. That is why `saturated_from` is clamped
/// to at least 1: an alphabet whose nominal masses are all `<= 0` drives it to 0, and without the
/// clamp the saturating loop would credit the source with `|aa|` edges and shift every subsequent
/// offset by that much (`fused_offsets_match_the_naive_count` pins this).
fn fill_edge_offsets(edge_start: &mut [u32], graph_max: i32, aa: &[Aa]) -> u32 {
    let full = aa.len() as u32;
    let max_nominal = aa.iter().map(|a| a.nominal).max().unwrap_or(0).max(0);
    let saturated_from = max_nominal.min(graph_max + 1).max(1);
    let mut acc = 0u32;
    edge_start[0] = 0;
    for m in 1..saturated_from {
        edge_start[m as usize] = acc;
        let mut c = 0u32;
        for a in aa {
            if m - a.nominal >= 0 {
                c += 1;
            }
        }
        acc += c;
    }
    for m in saturated_from..=graph_max {
        edge_start[m as usize] = acc;
        acc += full;
    }
    edge_start[graph_max as usize + 1] = acc;
    acc
}

/// [`build_reverse_graph`] into a caller-owned [`Graph`], so a whole run reuses one set of buffers.
///
/// Identical output to `build_reverse_graph` — every element of every array is overwritten, so the
/// incoming contents of `graph` never affect the result. This is the graph-side analogue of
/// [`crate::DpScratch`]: the five CSR buffers are ~0.5 MB per spectrum on the F13 grid, which is
/// large enough that glibc services them with `mmap`/`munmap` and the writes fault in fresh zero
/// pages every time. Hold one `Graph` across the whole run and the per-spectrum allocation is nil.
pub fn build_reverse_graph_into(
    graph: &mut Graph,
    scored: &ScoredSpectrum,
    tables: &SpectrumTables,
    complement_mass: i32,
    sinks: &[i32],
    aa: &[Aa],
    cleavage: PeptideCleavage<'_>,
) -> Vec<usize> {
    assert!(
        aa.len() <= u16::MAX as usize,
        "graph edge alphabet is limited to {} amino acids (got {})",
        u16::MAX,
        aa.len()
    );
    let graph_max = sinks
        .iter()
        .copied()
        .max()
        .unwrap_or(complement_mass)
        .max(complement_mass)
        .max(0);
    let n = graph_max as usize + 1;

    // Node scores: exactly the per-candidate computation, so share the one implementation.
    graph.node_score.resize(n, 0);
    graph.recompute_node_scores(tables, complement_mass, sinks);

    // The probability table the per-edge `edge_aa` index addresses.
    graph.aa_prob.clear();
    graph.aa_prob.extend(aa.iter().map(|a| a.prob));

    // CSR offsets. The incoming-edge count of node `m` is a pure function of the node mass — the
    // number of amino acids with `nominal <= m` — so it saturates at `|aa|` for every node at or
    // above the heaviest residue (~186 of the ~1500 nodes are below it). Count those few by
    // scanning the alphabet and fill the rest from the saturated count, instead of rescanning the
    // alphabet for every node and then prefix-summing in a second pass.
    graph.edge_start.resize(n + 1, 0);
    let total = fill_edge_offsets(&mut graph.edge_start, graph_max, aa) as usize;

    // Fill edges in amino-acid order (the DP's summation order — preserved for bit-exact
    // reproduction). Edge scores are **candidate-independent** — they depend only on the shared
    // node masses and the amino acid, never on `complement_mass` or which node is the sink — so the
    // edge arrays are built once and reused across the isotope-error candidates (see
    // [`Graph::recompute_node_scores`]). The sink's edges (errorScore 0, `setBackwardEdgesFromSink`)
    // are zeroed in the DP instead of here, since the sink differs per candidate.
    //
    // Node masses come from the shared `tables` (one peak lookup per node per spectrum, not one per
    // incident edge per graph). `m, prev` are always in-range, so the `edge_score` bounds guard
    // never fires; `edge_score_with` on the resolved masses is equivalent.
    graph.edge_prev.resize(total, 0);
    graph.edge_score.resize(total, 0);
    graph.edge_aa.resize(total, 0);
    if total > 0 {
        // `getEdgeScoreInt` only consults the mass error when *both* endpoints matched a peak
        // (`index == 3`); otherwise it is `round(ionExistenceScore(index))`, a per-spectrum
        // constant. Most nodes on this grid match no peak, so resolve those three constants once
        // and keep the arithmetic for the edges that actually need it. Same integers either way.
        let base = [
            scored.edge_score_with(-1.0, -1.0, 0.0), // index 0: neither endpoint matched
            scored.edge_score_with(0.0, -1.0, 0.0),  // index 1: cur matched, prev did not
            scored.edge_score_with(-1.0, 0.0, 0.0),  // index 2: prev matched, cur did not
        ];
        let node_mass = &tables.node_mass;
        let Graph {
            edge_start,
            edge_prev,
            edge_score,
            edge_aa,
            ..
        } = graph;
        for m in 1..=graph_max {
            let mut pos = edge_start[m as usize] as usize;
            let cur_mass = node_mass[m as usize];
            let cur_matched = cur_mass >= 0.0;
            for (ai, a) in aa.iter().enumerate() {
                let prev = m - a.nominal;
                if prev < 0 {
                    continue;
                }
                let prev_mass = node_mass[prev as usize];
                let mut es = if cur_matched && prev_mass >= 0.0 {
                    scored.edge_score_with(cur_mass, prev_mass, a.accurate_mass)
                } else if cur_matched {
                    base[1]
                } else if prev_mass >= 0.0 {
                    base[2]
                } else {
                    base[0]
                };
                if prev == 0 {
                    es += cleavage.score(a.residue);
                }
                edge_prev[pos] = prev as u32;
                edge_score[pos] = es;
                edge_aa[pos] = ai as u16;
                pos += 1;
            }
        }
    }

    sinks.iter().map(|&s| s as usize).collect()
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
        // The loop below assigns every node in `1..graph_max`, so only the source, the skipped
        // nodes, and the tail above `graph_max` need clearing — a blanket memset of the whole array
        // (which is then immediately overwritten) is wasted work, not correctness.
        let ns = &mut self.node_score;
        if let Some(s) = ns.first_mut() {
            *s = 0;
        }
        for m in 1..graph_max {
            ns[m as usize] = if is_sink(m) || m >= complement_mass {
                0 // sinks and out-of-range nodes score 0
            } else {
                msgf_chem::round_half_up(
                    tables.prefix[(complement_mass - m) as usize] + tables.suffix[m as usize],
                )
            };
        }
        let tail = (graph_max.max(0) as usize).min(ns.len());
        for s in ns[tail..].iter_mut() {
            *s = 0;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    fn aa_with(nominals: &[i32]) -> Vec<Aa> {
        nominals
            .iter()
            .map(|&nominal| Aa {
                residue: b'A',
                nominal,
                accurate_mass: nominal as f32,
                prob: 0.05,
            })
            .collect()
    }

    /// The saturating offset fill must agree with a naive per-node count on every alphabet shape,
    /// including degenerate ones. A zero/negative-nominal alphabet drives `saturated_from` to 0 and
    /// used to credit the source (node 0) with a full edge list, shifting every offset and leaving
    /// one stale edge slot at the front.
    #[test]
    fn fused_offsets_match_the_naive_count() {
        let cases: Vec<(&str, Vec<i32>)> = vec![
            (
                "standard",
                standard_aa_nominal().into_iter().map(|(_, n)| n).collect(),
            ),
            ("all-zero nominal", vec![0]),
            ("zero and positive", vec![0, 5, 9]),
            ("negative nominal", vec![-3, -1]),
            ("single heavy", vec![186]),
            ("heavier than the graph", vec![500]),
            ("duplicate nominals", vec![113, 113, 128, 128]),
        ];
        for (name, nominals) in cases {
            let aa = aa_with(&nominals);
            for graph_max in [1i32, 2, 57, 100, 200] {
                let n = graph_max as usize + 1;
                // Naive reference: node 0 has no incoming edges; node m has one per `nominal <= m`.
                let mut want = vec![0u32; n + 1];
                let mut acc = 0u32;
                for m in 0..=graph_max {
                    want[m as usize] = acc;
                    if m >= 1 {
                        acc += aa.iter().filter(|a| m - a.nominal >= 0).count() as u32;
                    }
                }
                want[n] = acc;

                let mut got = vec![u32::MAX; n + 1];
                let total = fill_edge_offsets(&mut got, graph_max, &aa);
                assert_eq!(got, want, "offsets differ: {name}, graph_max={graph_max}");
                assert_eq!(total, acc, "total differs: {name}, graph_max={graph_max}");
                assert_eq!(got[1], 0, "source must have no incoming edges: {name}");
            }
        }
    }

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
