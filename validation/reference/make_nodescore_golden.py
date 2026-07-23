#!/usr/bin/env python3
"""
make_nodescore_golden.py — distill the authoritative MS-GF+ node-score dumps (from the Java
ScoreDumper, which calls NewRankScorer.getNodeScore / getMissingIonScore) into a golden for the
Rust scoring primitives (ScoringModel::node_score / missing_ion_score).

Every (partition, ion, sampled rank | MISSING) -> score row is kept, per model, plus a grand
score sum as a cheap global check. Usage:
  make_nodescore_golden.py <m1.nodescore.tsv> [<m2.tsv> ...] -o <out.golden.json>
"""
import argparse, json, os


def parse_tsv(path):
    rows = []
    total = 0.0
    with open(path) as fh:
        next(fh)  # header
        for line in fh:
            pi, charge, seg, pmass, ion, ion_charge, rank, score = line.rstrip("\n").split("\t")
            score = float(score)
            rows.append([int(pi), ion, int(ion_charge), rank, score])
            total += score
    # model file name: strip ".nodescore" -> ".param"
    base = os.path.basename(path).replace(".nodescore.tsv", ".param")
    return base, {"count": len(rows), "score_sum": total, "rows": rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tsv", nargs="+")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    models = {}
    for t in a.tsv:
        name, data = parse_tsv(t)
        models[name] = data
        print(f"  {name}: {data['count']} rows, score_sum={data['score_sum']:.4f}")
    out = {
        "note": "MS-GF+ getNodeScore/getMissingIonScore per (partition, ion, rank); rank 'MISSING' = absent peak bin",
        "ranks_sampled": [1, 2, 3, 5, 10, 50, 100, 149, 150, 151, "MISSING"],
        "row_format": ["partition_index", "ion", "ion_charge", "rank", "score"],
        "compare": {"score": "abs 1e-4", "score_sum": "abs 1e-1"},
        "models": models,
    }
    json.dump(out, open(a.out, "w"))
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
