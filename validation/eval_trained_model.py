#!/usr/bin/env python3
"""Head-to-head evaluation of two `.param` scoring models on the same spectra.

Both models rescore an identical PSM list containing, for every spectrum, the identified peptide
(**target**) and a mass-identical shuffled version of it (**decoy**). A model is good insofar as it
pushes targets away from decoys, so the metrics are decoy-referenced and need no external FDR
oracle:

* `target<decoy`   — fraction of spectra where the target's SpecEValue beats its own decoy's
* `median gap`     — median log10(decoy SpecEValue / target SpecEValue)
* `IDs @ 1% decoy` — targets scoring better than the 1st-percentile decoy SpecEValue (an ID count
                     at a decoy-calibrated threshold: the ID-parity gate, computed per model)

Usage:
    eval_trained_model.py library --mgf held_out.mgf --models A.param B.param [--n 2000]
    eval_trained_model.py f13     --tsv golden.tsv --mgf F13.mgf --models A.param B.param
"""
from __future__ import annotations

import argparse
import bisect
import math
import random
import re
import subprocess
import sys
import time
from pathlib import Path

MOD = re.compile(r"[+-]\d+(?:\.\d+)?")
BIN = Path(__file__).resolve().parent.parent / "rust/target/release/msgf"


def normalize_peptide(seq: str) -> str | None:
    """`SEQ=` form -> the CLI's peptide form (a leading N-term delta moves onto residue 1)."""
    seq = seq.strip()
    m = MOD.match(seq)
    nterm = ""
    if m:
        nterm, seq = m.group(0), seq[m.end() :]
    if not seq or not seq[0].isalpha():
        return None
    out = seq[0] + nterm + seq[1:]
    return out if all(c.isalpha() or c in "+-." or c.isdigit() for c in out) else None


def split_residues(pep: str) -> list[str]:
    """Split a peptide into residue+mod tokens."""
    toks, i = [], 0
    while i < len(pep):
        if not pep[i].isalpha():
            return []
        j = i + 1
        while j < len(pep) and not pep[j].isalpha():
            j += 1
        toks.append(pep[i:j])
        i = j
    return toks


def shuffle_decoy(pep: str, rng: random.Random) -> str | None:
    """Mass-identical decoy: permute all but the C-terminal residue (with its mods attached)."""
    toks = split_residues(pep)
    if len(toks) < 5:
        return None
    body, last = toks[:-1], toks[-1]
    for _ in range(20):
        shuffled = body[:]
        rng.shuffle(shuffled)
        if shuffled != body:
            return "".join(shuffled) + last
    return None


def read_library(mgf: Path, n: int, seed: int, n_decoys: int = 1):
    """Sample annotated spectra from a library MGF -> (spectrum blocks, psm rows)."""
    blocks, cur, scan, seq = [], [], None, None
    for line in open(mgf):
        if line.startswith("BEGIN IONS"):
            cur, scan, seq = [line], None, None
        elif line.startswith("END IONS"):
            cur.append(line)
            if scan and seq:
                blocks.append((scan, seq, cur))
        else:
            if line.startswith("SCANS="):
                scan = line[6:].strip()
            elif line.startswith("SEQ="):
                seq = line[4:].strip()
            if cur is not None:
                cur.append(line)
    rng = random.Random(seed)
    rng.shuffle(blocks)
    out_blocks, rows = [], []
    for scan, seq, block in blocks:
        pep = normalize_peptide(seq)
        if not pep:
            continue
        charge = next((l[7:].strip().rstrip("+") for l in block if l.startswith("CHARGE=")), None)
        decoys = {d for d in (shuffle_decoy(pep, rng) for _ in range(n_decoys)) if d and d != pep}
        if not charge or not decoys:
            continue
        out_blocks.append(block)
        rows.append((scan, pep, charge, "target"))
        for d in decoys:
            rows.append((scan, d, charge, "decoy"))
        if len(out_blocks) >= n:
            break
    return out_blocks, rows


def read_f13(tsv: Path, n: int, seed: int, n_decoys: int = 1):
    """PSM rows from the MS-GF+ result TSV (one best match per scan) + shuffled decoys."""
    rng = random.Random(seed)
    best: dict[str, tuple[float, str, str]] = {}
    with open(tsv) as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        ci = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            scan, pep, charge = f[ci["ScanNum"]], f[ci["Peptide"]], f[ci["Charge"]]
            ev = float(f[ci["SpecEValue"]])
            pep = pep.split(".", 1)[-1].rsplit(".", 1)[0] if pep.count(".") >= 2 else pep
            if scan not in best or ev < best[scan][0]:
                best[scan] = (ev, pep, charge)
    rows, msgf_ev = [], {}
    for scan, (ev, pep, charge) in sorted(best.items(), key=lambda kv: int(kv[0])):
        decoys = {d for d in (shuffle_decoy(pep, rng) for _ in range(n_decoys)) if d and d != pep}
        if not decoys:
            continue
        msgf_ev[scan] = ev
        rows.append((scan, pep, charge, "target"))
        for d in decoys:
            rows.append((scan, d, charge, "decoy"))
        if len(msgf_ev) >= n:
            break
    return rows, msgf_ev


def run_rescore(model: Path, mgf: Path, psms: Path, out: Path, extra: list[str]) -> float:
    t0 = time.time()
    cmd = [str(BIN), "rescore", "-s", str(mgf), "-p", str(model), "-i", str(psms), "-o", str(out)]
    r = subprocess.run(cmd + extra, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"rescore failed for {model}:\n{r.stderr}")
    return time.time() - t0


def analyse(out: Path, labels: dict[tuple[str, str], str], subset: set[str] | None = None):
    """-> (per-scan target/decoy SpecEValue+RawScore, aggregate metrics)."""
    tgt: dict[str, tuple[float, int]] = {}
    dec: dict[str, tuple[float, int]] = {}
    with open(out) as fh:
        ci = {c: i for i, c in enumerate(fh.readline().rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            scan, pep = f[ci["scan"]], f[ci["peptide"]]
            ev, raw = float(f[ci["spec_evalue"]]), int(f[ci["raw_score"]])
            kind = labels.get((scan, pep))
            if subset is not None and scan not in subset:
                continue
            if kind == "target":
                tgt[scan] = (ev, raw)
            elif kind == "decoy" and (scan not in dec or ev < dec[scan][0]):
                dec[scan] = (ev, raw)  # best (smallest) decoy for the spectrum
    both = sorted(set(tgt) & set(dec))
    wins = sum(tgt[s][0] < dec[s][0] for s in both)
    gaps = [
        math.log10(dec[s][0] / tgt[s][0])
        for s in both
        if tgt[s][0] > 0 and dec[s][0] > 0
    ]
    gaps.sort()
    dvals = sorted(dec[s][0] for s in both)
    thresh = dvals[max(0, int(0.01 * len(dvals)) - 1)] if dvals else 0.0
    ids = sum(tgt[s][0] < thresh for s in both)
    tev = sorted(math.log10(tgt[s][0]) for s in both if tgt[s][0] > 0)
    return (
        tgt,
        dec,
        {
            "scans": len(both),
            "median_target_log10_specE": tev[len(tev) // 2] if tev else float("nan"),
            "target<decoy": wins / len(both) if both else 0.0,
            "median_gap_log10": gaps[len(gaps) // 2] if gaps else float("nan"),
            "decoy_1pct_threshold": thresh,
            "ids_at_1pct_decoy": ids,
            "id_rate": ids / len(both) if both else 0.0,
            "median_target_raw": sorted(v[1] for v in tgt.values())[len(tgt) // 2] if tgt else 0,
            "median_decoy_raw": sorted(v[1] for v in dec.values())[len(dec) // 2] if dec else 0,
        },
    )


def spearman(x: list[float], y: list[float]) -> float:
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        for pos, i in enumerate(order):
            r[i] = pos
        return r

    rx, ry = ranks(x), ranks(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["library", "f13"])
    ap.add_argument("--mgf", required=True, type=Path)
    ap.add_argument("--tsv", type=Path)
    ap.add_argument("--models", nargs="+", required=True, type=Path)
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--decoys", type=int, default=1, help="decoys per spectrum")
    ap.add_argument("--seed", type=int, default=20260724)
    ap.add_argument("--work", type=Path, default=Path("/tmp/msgf_eval"))
    ap.add_argument("--aa-probs", type=Path)
    ap.add_argument("--ox-m", action="store_true")
    args = ap.parse_args()

    args.work.mkdir(parents=True, exist_ok=True)
    msgf_ev: dict[str, float] = {}
    if args.mode == "library":
        blocks, rows = read_library(args.mgf, args.n, args.seed, args.decoys)
        mgf = args.work / "eval.mgf"
        with open(mgf, "w") as fh:
            for b in blocks:
                fh.writelines(b)
    else:
        rows, msgf_ev = read_f13(args.tsv, args.n, args.seed, args.decoys)
        mgf = args.mgf

    psms = args.work / "eval_psms.tsv"
    with open(psms, "w") as fh:
        fh.write("scan\tpeptide\tcharge\n")
        for scan, pep, charge, _ in rows:
            fh.write(f"{scan}\t{pep}\t{charge}\n")
    labels = {(s, p): k for s, p, _, k in rows}
    n_spec = sum(1 for r in rows if r[3] == "target")
    print(f"{n_spec} spectra, {len(rows)} PSMs ({args.decoys} shuffled decoy(s) per target)")

    extra = []
    if args.aa_probs:
        extra += ["--aa-probs", str(args.aa_probs)]
    if args.ox_m:
        extra.append("--ox-m")

    results = {}
    for model in args.models:
        out = args.work / f"rescored_{model.stem}.tsv"
        secs = run_rescore(model, mgf, psms, out, extra)
        tgt, dec, m = analyse(out, labels)
        m["seconds"] = secs
        m["psms_per_sec"] = len(rows) / secs if secs else 0
        results[model.stem] = (tgt, dec, m)

    keys = [
        "scans",
        "median_target_log10_specE",
        "target<decoy",
        "median_gap_log10",
        "ids_at_1pct_decoy",
        "id_rate",
        "median_target_raw",
        "median_decoy_raw",
        "seconds",
        "psms_per_sec",
    ]
    names = list(results)
    w = max(24, max(len(n) for n in names) + 2)
    print("\n" + " " * 22 + "".join(f"{n:>{w}}" for n in names))
    for k in keys:
        vals = []
        for n in names:
            v = results[n][2][k]
            vals.append(f"{v:>{w}.4f}" if isinstance(v, float) else f"{v:>{w}}")
        print(f"{k:<22}" + "".join(vals))

    # Stratify F13 by MS-GF+'s own confidence: the peptides were *chosen* by the reference model,
    # so any deficit that lives only in the low-confidence stratum is selection bias, not model
    # quality.
    if msgf_ev:
        conf = sorted(msgf_ev.items(), key=lambda kv: kv[1])
        cut = max(1, len(conf) // 5)
        strata = {"top 20% by MS-GF+ SpecEValue": {s for s, _ in conf[:cut]},
                  "bottom 80%": {s for s, _ in conf[cut:]}}
        for label, subset in strata.items():
            print(f"\n  [{label}]  n={len(subset)}")
            for n in names:
                out = args.work / f"rescored_{n}.tsv"
                _, _, m = analyse(out, labels, subset)
                print(f"    {n:<28} target<decoy {m['target<decoy']:.3f}  "
                      f"median gap {m['median_gap_log10']:6.2f}  "
                      f"median target log10 E {m['median_target_log10_specE']:7.2f}  "
                      f"IDs@1%decoy {m['ids_at_1pct_decoy']:>5}")

    # --- calibration: SpecEValue is the probability that a random peptide of the same mass scores
    # at least this well, so for the shuffled decoys P(SpecEValue < x) should be about x. A model
    # whose curve sits above the diagonal is over-confident (its null is too pessimistic).
    print("\nSpecEValue calibration on decoys — empirical P(E < x), ideal = x")
    xs = [1e-4, 1e-6, 1e-8, 1e-10]
    print(f"{'x':>12}" + "".join(f"{n:>26}" for n in names))
    for x in xs:
        cells = []
        for n in names:
            out = args.work / f"rescored_{n}.tsv"
            vals = []
            with open(out) as fh:
                ci = {c: i for i, c in enumerate(fh.readline().rstrip("\n").split("\t"))}
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if labels.get((f[ci["scan"]], f[ci["peptide"]])) == "decoy":
                        vals.append(float(f[ci["spec_evalue"]]))
            frac = sum(v < x for v in vals) / max(len(vals), 1)
            cells.append(f"{frac:>18.2e} ({frac / x:>4.1f}x)")
        print(f"{x:>12.0e}" + "".join(f"{c:>26}" for c in cells))

    if len(names) == 2:
        a, b = names
        ta, tb = results[a][0], results[b][0]
        common = sorted(set(ta) & set(tb))
        ea = [math.log10(ta[s][0]) for s in common if ta[s][0] > 0 and tb[s][0] > 0]
        eb = [math.log10(tb[s][0]) for s in common if ta[s][0] > 0 and tb[s][0] > 0]
        ra = [ta[s][1] for s in common]
        rb = [tb[s][1] for s in common]
        print(
            f"\nagreement on targets (n={len(ea)}): "
            f"spearman log10 SpecEValue {spearman(ea, eb):.4f}, "
            f"spearman RawScore {spearman(ra, rb):.4f}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
