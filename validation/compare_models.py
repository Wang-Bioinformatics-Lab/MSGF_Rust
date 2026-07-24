#!/usr/bin/env python3
"""Compare two MS-GF+ `.param` scoring models table-for-table.

An independent decoder written from `docs/param-format.md` (deliberately *not* a binding to the
Rust reader — decoding the file twice, in two languages, is what makes the comparison evidence).

Usage:
    python3 validation/compare_models.py A.param B.param [--json out.json]

Prints identity/config differences, the partition schemes, ion-type selection overlap, and the
quantity that actually drives scoring: the node score `ln(ion[rank] / (noise[rank] * k))` that
`ScoringModel::score_from_table` computes.
"""
from __future__ import annotations

import argparse
import json
import math
import struct
import sys
from dataclasses import dataclass, field


class R:
    """Big-endian reader mirroring Java DataInputStream."""

    def __init__(self, b: bytes):
        self.b, self.p = b, 0

    def take(self, n: int) -> bytes:
        s = self.b[self.p : self.p + n]
        if len(s) != n:
            raise EOFError(f"short read at {self.p}")
        self.p += n
        return s

    def i32(self) -> int:
        return struct.unpack(">i", self.take(4))[0]

    def f32(self) -> float:
        return struct.unpack(">f", self.take(4))[0]

    def u8(self) -> int:
        return self.take(1)[0]

    def boolean(self) -> bool:
        return self.u8() != 0

    def jstring(self) -> str | None:
        n = self.u8()
        if n == 0:
            return None
        return self.take(2 * n).decode("utf-16-be")


@dataclass
class Model:
    path: str = ""
    version: int = 0
    activation: str = ""
    instrument: str = ""
    enzyme: str | None = None
    protocol: str | None = None
    mme_ppm: bool = False
    mme: float = 0.0
    deconv: bool = False
    deconv_tol: float = 0.0
    charge_hist: list = field(default_factory=list)
    num_segments: int = 0
    partitions: list = field(default_factory=list)  # (charge, mass, seg)
    precursor_off: list = field(default_factory=list)
    frag_off: list = field(default_factory=list)  # per partition: [(name, prefix, charge, off, freq)]
    max_rank: int = 0
    rank_dist: dict = field(default_factory=dict)  # partition index -> {name: [floats]}
    esf: int = 0
    error_dist: list = field(default_factory=list)  # per partition: (signal, noise, ionexist)
    size: int = 0


def read_param(path: str) -> Model:
    raw = open(path, "rb").read()
    r = R(raw)
    m = Model(path=path, size=len(raw))
    m.version = r.i32()
    m.activation = r.jstring()
    m.instrument = r.jstring()
    m.enzyme = r.jstring()
    m.protocol = r.jstring()
    m.mme_ppm, m.mme = r.boolean(), r.f32()
    m.deconv, m.deconv_tol = r.boolean(), r.f32()

    m.charge_hist = [(r.i32(), r.i32()) for _ in range(r.i32())]

    n = r.i32()
    m.num_segments = r.i32()
    parts = [(r.i32(), r.f32(), r.i32()) for _ in range(n)]
    # canonical order is (charge, seg, mass)
    m.partitions = sorted({(c, s, mass) for c, mass, s in parts})
    m.partitions = [(c, mass, s) for c, s, mass in m.partitions]

    m.precursor_off = [
        (r.i32(), r.i32(), r.f32(), r.boolean(), r.f32(), r.f32()) for _ in range(r.i32())
    ]

    for _ in m.partitions:
        blk = []
        for _ in range(r.i32()):
            prefix, charge, off, freq = r.boolean(), r.i32(), r.f32(), r.f32()
            name = f"{'P' if prefix else 'S'}_{charge}_{math.floor(off + 0.5):.0f}"
            blk.append((name, prefix, charge, off, freq))
        m.frag_off.append(blk)

    m.max_rank = r.i32()
    cols = m.max_rank + 1
    for pi, blk in enumerate(m.frag_off):
        if not blk:
            continue
        rows = {}
        for name, *_ in blk:
            rows[name] = [r.f32() for _ in range(cols)]
        rows["noise"] = [r.f32() for _ in range(cols)]
        m.rank_dist[pi] = rows

    m.esf = r.i32()
    if m.esf > 0:
        w = 2 * m.esf + 1
        for _ in m.partitions:
            sig = [r.f32() for _ in range(w)]
            noi = [r.f32() for _ in range(w)]
            ie = [r.f32() for _ in range(4)]
            m.error_dist.append((sig, noi, [x if x != 0.0 else 0.001 for x in ie]))

    term = r.i32()
    assert term == 0x7FFFFFFF, f"bad terminator {term:#x} at {r.p - 4}"
    assert r.p == len(raw), f"{len(raw) - r.p} trailing bytes"
    return m


def node_score(m: Model, pi: int, name: str, rank: int) -> float | None:
    """`ScoringModel::score_from_table` — the number the scorer actually adds up.

    `rank` is 1-based; pass `rank = max_rank + 1` for the 'ion absent' bin.
    """
    rows = m.rank_dist.get(pi)
    if not rows or name not in rows:
        return None
    idx = m.max_rank if rank > m.max_rank + 0.5 and rank == m.max_rank + 1 else (
        m.max_rank - 1 if rank > m.max_rank else rank - 1
    )
    ion_charge = next(c for n, _, c, _, _ in m.frag_off[pi] if n == name)
    denom = rows["noise"][idx] * min(ion_charge, m.num_segments)
    if denom <= 0 or rows[name][idx] <= 0:
        return None
    return math.log(rows[name][idx] / denom)


def partition_for(m: Model, charge: int, mass: float, seg: int) -> int | None:
    """The scorer's TreeSet floor lookup over (charge, seg, parent_mass)."""

    def floor(c, s, ms):
        best = None
        for i, (pc, pm, ps) in enumerate(m.partitions):
            if (pc, ps, pm) <= (c, s, ms):
                best = i
        return best

    i = floor(charge, seg, mass)
    if i is None:
        return floor(m.partitions[0][0], seg, mass) if m.partitions else None
    mc = m.partitions[i][0]
    return i if mc == charge else floor(mc, seg, mass)


def describe(m: Model) -> dict:
    ion_hist: dict[str, int] = {}
    for blk in m.frag_off:
        for name, *_ in blk:
            ion_hist[name] = ion_hist.get(name, 0) + 1
    charges = sorted({c for c, _, _ in m.partitions})
    return {
        "file": m.path,
        "bytes": m.size,
        "identity": f"{m.activation}/{m.instrument}/{m.enzyme}/{m.protocol or 'Automatic'}",
        "mme": f"{m.mme} {'ppm' if m.mme_ppm else 'Da'}",
        "deconvolution": (m.deconv, m.deconv_tol),
        "num_segments": m.num_segments,
        "max_rank": m.max_rank,
        "error_scaling_factor": m.esf,
        "partitions": len(m.partitions),
        "partitions_per_charge": {c: len({p[1] for p in m.partitions if p[0] == c}) for c in charges},
        "training_spectra": sum(n for _, n in m.charge_hist),
        "precursor_offsets": len(m.precursor_off),
        "ion_types": dict(sorted(ion_hist.items(), key=lambda kv: -kv[1])),
        "mean_ions_per_partition": round(
            sum(len(b) for b in m.frag_off) / max(len(m.frag_off), 1), 2
        ),
    }


def curve(m: Model, charge: int, mass: float, seg: int, name: str, ranks) -> list:
    pi = partition_for(m, charge, mass, seg)
    if pi is None:
        return []
    return [node_score(m, pi, name, r) for r in ranks]


def error_width(m: Model, pi: int) -> tuple[float, float]:
    """Fraction of the signal error distribution inside ±1 and ±5 bins of centre."""
    sig = m.error_dist[pi][0]
    c = len(sig) // 2
    return (
        sum(sig[c - 1 : c + 2]) / max(sum(sig), 1e-12),
        sum(sig[c - 5 : c + 6]) / max(sum(sig), 1e-12),
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--json")
    args = ap.parse_args()

    A, B = read_param(args.a), read_param(args.b)
    da, db = describe(A), describe(B)

    w = 34
    print(f"{'':<26}{'A: ' + A.path.split('/')[-1]:<{w}}{'B: ' + B.path.split('/')[-1]:<{w}}")
    print("-" * (26 + 2 * w))
    for k in da:
        if k == "file":
            continue
        va, vb = str(da[k]), str(db[k])
        if len(va) > w - 2 or len(vb) > w - 2:
            print(f"{k:<26}\n    A  {va}\n    B  {vb}")
        else:
            print(f"{k:<26}{va:<{w}}{vb:<{w}}")

    # --- ion-type selection overlap
    ia, ib = set(da["ion_types"]), set(db["ion_types"])
    print(f"\nion types: shared {sorted(ia & ib)}")
    print(f"           only A {sorted(ia - ib)}   only B {sorted(ib - ia)}")

    # --- the scoring curve, in a mid-range charge-2 partition of each model
    ranks = [1, 2, 3, 5, 10, 20, 50, 100, 150, 151]  # 151 = the 'absent' bin
    print("\nnode score ln(ion/noise) for charge 2, parent mass 1200, segment 0")
    print(f"{'rank':>6}  " + "".join(f"{r:>9}" for r in ranks[:-1]) + f"{'absent':>9}")
    rows = {}
    for label, m in (("A", A), ("B", B)):
        for name in ("S_1_19", "P_1_1"):
            c = curve(m, 2, 1200.0, 0, name, ranks)
            if not c:
                continue
            rows[f"{label}:{name}"] = c
            print(
                f"{label}:{name:<8}"
                + "".join(f"{v:9.2f}" if v is not None else f"{'-':>9}" for v in c)
            )

    # --- are the two models' tables even on the same scale?
    print("\nrank-row scale (charge 2, mass 1200, segment 0): a row sums to the average number of")
    print("scored sites per spectrum; 'hit' is the fraction of those sites that matched a peak")
    for label, m in (("A", A), ("B", B)):
        pi = partition_for(m, 2, 1200.0, 0)
        rows = m.rank_dist.get(pi, {})
        for name in ("S_1_19", "P_1_1", "noise"):
            if name not in rows:
                continue
            row = rows[name]
            total, present = sum(row), sum(row[:-1])
            print(
                f"  {label}:{name:<8} sum {total:6.3f}  present {present:6.3f}  "
                f"absent {row[-1]:6.3f}  hit {present / total:5.3f}"
            )

    # --- mass-error distribution sharpness (the high-res term)
    print("\nsignal mass-error concentration (charge 2, mass 1200, last segment)")
    for label, m in (("A", A), ("B", B)):
        pi = partition_for(m, 2, 1200.0, m.num_segments - 1)
        if pi is not None and m.error_dist:
            e1, e5 = error_width(m, pi)
            ie = m.error_dist[pi][2]
            print(
                f"  {label}: within +-1 bin {e1:6.3f}   within +-5 bins {e5:6.3f}   "
                f"ion-existence {[round(x, 4) for x in ie]}"
            )

    # --- precursor offsets
    print("\nprecursor offsets (charge 2, reduced charge 0)")
    for label, m in (("A", A), ("B", B)):
        got = [(round(o, 4), round(f, 3)) for c, rc, o, _, _, f in m.precursor_off if c == 2 and rc == 0]
        print(f"  {label}: {got}")

    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"a": da, "b": db, "curves": rows}, fh, indent=2)
        print(f"\nwrote {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
