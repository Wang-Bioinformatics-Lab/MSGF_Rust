//! Test-only reference for the fixed N-terminal label (`NullModel::nterm_delta`): the plain,
//! unblocked generating function with the label as first written in DIA_Proteomics_Rust (branch
//! `dda-tmt`, `rust/src/dda/specprob/null.rs` at commit `26079085`; MIT OR Apache-2.0, same
//! author, the clean-room code this crate came from, extended there). The fast [`super::build`]
//! must reproduce it bit for bit.

#![allow(clippy::all)]

use super::{nominal, NullModel, NullTail, PreparedSpectrum, SinkDetail};

struct RefEdge {
    nom: usize,
    mass: f32,
    prob: f64,
    term: i32,
}

pub(super) fn ref_build(
    prep: &PreparedSpectrum,
    null: &NullModel,
    min_query: Option<i32>,
    detail: bool,
) -> Option<NullTail> {
    let edges: Vec<RefEdge> = null
        .alphabet
        .iter()
        .filter_map(|a| {
            let n = nominal(a.mass);
            (n >= 1).then(|| RefEdge {
                nom: n as usize,
                mass: a.mass as f32,
                prob: a.prob,
                term: null.cleavage.term(a.letter),
            })
        })
        .collect();
    let cl = &null.cleavage;
    let (kmin, kmax) = if cl.enabled {
        (cl.credit.min(cl.penalty), cl.credit.max(cl.penalty))
    } else {
        (0, 0)
    };
    let pi = null.mixture_pi();
    // Edges into the sink (the N-terminal residue) with the fixed N-terminal label (none = the
    // regular alphabet lands on the sink, as without a label).
    let nedges: Option<Vec<RefEdge>> = (null.nterm_delta != 0.0).then(|| {
        null.alphabet
            .iter()
            .filter_map(|a| {
                (nominal(a.mass) >= 1).then(|| {
                    let m = a.mass + null.nterm_delta;
                    RefEdge {
                        nom: nominal(m).max(1) as usize,
                        mass: m as f32,
                        prob: a.prob,
                        term: null.cleavage.term(a.letter),
                    }
                })
            })
            .collect()
    });

    let n0 = prep.n0();
    let (ilo, ihi) = null.isotope;
    let mut per_sink: Vec<(i32, i32, Vec<f64>, i32, Vec<i32>)> = Vec::new(); // sink, G lowest, G, best, v
    for p in (n0 - ihi)..=(n0 - ilo) {
        if p <= 0 {
            continue;
        }
        let cut = min_query.map(|q| q - kmax);
        let Some(sd) = ref_sink_dp(prep, &edges, nedges.as_deref(), p, cut) else {
            continue;
        };
        let (lo, g) = if cl.enabled {
            let glo = sd.low + kmin;
            let ghi = sd.best + kmax;
            let mut g = vec![0f64; (ghi - glo + 1) as usize];
            let dget = |s: i32| -> f64 {
                let i = s - sd.low;
                if i < 0 || i as usize >= sd.cells.len() {
                    0.0
                } else {
                    sd.cells[i as usize]
                }
            };
            for (i, c) in g.iter_mut().enumerate() {
                let s = glo + i as i32;
                *c = pi * dget(s - cl.credit) + (1.0 - pi) * dget(s - cl.penalty);
            }
            (glo, g)
        } else {
            (sd.low, sd.cells)
        };
        per_sink.push((
            p,
            lo,
            g,
            sd.best,
            if detail { sd.vertex } else { Vec::new() },
        ));
    }
    if per_sink.is_empty() {
        return None;
    }
    let lowest = per_sink.iter().map(|s| s.1).min().unwrap();
    let top = per_sink
        .iter()
        .map(|s| s.1 + s.2.len() as i32 - 1)
        .max()
        .unwrap();
    let mut mass = vec![0f64; (top - lowest + 1) as usize];
    for s in &per_sink {
        let off = (s.1 - lowest) as usize;
        for (i, v) in s.2.iter().enumerate() {
            mass[off + i] += v;
        }
    }
    let best = per_sink.iter().map(|s| s.3 + kmax).max().unwrap();
    let sinks = if detail {
        per_sink
            .into_iter()
            .map(
                |(sink, lowest, mass, best_path, vertex_weights)| SinkDetail {
                    sink,
                    vertex_weights,
                    lowest,
                    mass,
                    best_path,
                },
            )
            .collect()
    } else {
        Vec::new()
    };
    Some(NullTail {
        lowest,
        mass,
        best,
        min_query,
        sinks,
    })
}

struct RefSinkDist {
    low: i32,
    cells: Vec<f64>,
    best: i32,
    vertex: Vec<i32>,
}

/// Unmixed distribution D_p of one sink. `cut`: drop cells whose best completion is below it.
/// `nedges`: the edges into the sink when a fixed N-terminal label shifts them (same order and
/// length as `edges`); `None` = the regular alphabet.
fn ref_sink_dp(
    prep: &PreparedSpectrum,
    edges: &[RefEdge],
    nedges: Option<&[RefEdge]>,
    p: i32,
    cut: Option<i32>,
) -> Option<RefSinkDist> {
    let pu = p as usize;
    let na = edges.len();
    // The alphabet of the edges INTO vertex m.
    let into = |m: usize| -> &[RefEdge] {
        match nedges {
            Some(n) if m == pu => n,
            _ => edges,
        }
    };
    // Vertex weights (index = suffix nominal mass) and anchor masses.
    let mut v = vec![0i32; pu + 1];
    for m in 1..pu {
        v[m] = prep.vertex(p, p - m as i32);
    }
    let alpha: Vec<f32> = (0..=p).map(|k| prep.alpha(k)).collect();
    // Edge scores e[m * na + a] (into m, label a); valid when m >= nom_a.
    let mut e = vec![0i32; (pu + 1) * na];
    for m in 1..pu {
        for (a, ed) in edges.iter().enumerate() {
            if m >= ed.nom {
                let mp = m - ed.nom;
                let mut s = prep.edge_between(alpha[m], alpha[mp], ed.mass);
                if mp == 0 {
                    s += ed.term;
                }
                e[m * na + a] = s;
            }
        }
    }
    // Forward: reachability and structural support [lo, hi].
    let mut reach = vec![false; pu + 1];
    let mut lo = vec![i32::MAX; pu + 1];
    let mut hi = vec![i32::MIN; pu + 1];
    reach[0] = true;
    lo[0] = 0;
    hi[0] = 0;
    for m in 1..=pu {
        for (a, ed) in into(m).iter().enumerate() {
            if m >= ed.nom && reach[m - ed.nom] {
                let mp = m - ed.nom;
                let sh = v[m] + e[m * na + a];
                reach[m] = true;
                lo[m] = lo[m].min(lo[mp] + sh);
                hi[m] = hi[m].max(hi[mp] + sh);
            }
        }
    }
    if !reach[pu] {
        return None;
    }
    // Backward: best completion to the sink.
    let mut rem = vec![i32::MIN; pu + 1];
    rem[pu] = 0;
    for m in (0..pu).rev() {
        if !reach[m] {
            continue;
        }
        for (a, ed) in edges.iter().enumerate() {
            let to = m + ed.nom;
            // With a label, a regular edge never lands on the sink (the shifted ones below do).
            if to <= pu && !(nedges.is_some() && to == pu) && rem[to] != i32::MIN {
                let cand = v[to] + e[to * na + a] + rem[to];
                if cand > rem[m] {
                    rem[m] = cand;
                }
            }
        }
        if let Some(n) = nedges {
            for (a, ed) in n.iter().enumerate() {
                if m + ed.nom == pu {
                    let cand = v[pu] + e[pu * na + a] + rem[pu];
                    if cand > rem[m] {
                        rem[m] = cand;
                    }
                }
            }
        }
    }
    let best = hi[pu];
    let cut = cut.map(|c| c.min(best));
    // Stored range per vertex.
    let mut low = vec![0i32; pu + 1];
    let mut start = vec![0usize; pu + 2];
    let mut total = 0usize;
    for m in 0..=pu {
        start[m] = total;
        if reach[m] && rem[m] != i32::MIN {
            let l = match cut {
                Some(c) => lo[m].max(c - rem[m]),
                None => lo[m],
            };
            low[m] = l;
            if hi[m] >= l {
                total += (hi[m] - l + 1) as usize;
            }
        }
    }
    start[pu + 1] = total;
    let mut arena = vec![0f64; total];
    arena[0] = 1.0; // D_0 = {0: 1}; low[0] = 0 always (the cut never exceeds the best path).
    for m in 1..=pu {
        let (ms, me) = (start[m], start[m + 1]);
        if ms == me {
            continue;
        }
        let lm = low[m];
        for (a, ed) in into(m).iter().enumerate() {
            if m < ed.nom {
                continue;
            }
            let mp = m - ed.nom;
            let (ps, pe) = (start[mp], start[mp + 1]);
            if ps == pe {
                continue;
            }
            let sh = v[m] + e[m * na + a];
            // Source cell i (score low[mp] + i) lands at target index low[mp] + i + sh - lm.
            let base = low[mp] + sh - lm;
            let skip = if base < 0 { (-base) as usize } else { 0 };
            let w = ed.prob;
            let (head, tail) = arena.split_at_mut(ms);
            let src = &head[ps..pe];
            let dst = &mut tail[..me - ms];
            if skip >= src.len() {
                continue;
            }
            let off = (base + skip as i32) as usize;
            for (d, s) in dst[off..].iter_mut().zip(&src[skip..]) {
                *d += w * s;
            }
        }
    }
    let cells = arena[start[pu]..start[pu + 1]].to_vec();
    Some(RefSinkDist {
        low: low[pu],
        cells,
        best,
        vertex: v,
    })
}
