#!/usr/bin/env python3
"""
run_regression.py — the MSGF_Rust regression suite.

Runs TODAY with no Java and no Rust: it re-derives every golden fixture from the raw data and
authoritative constants and asserts they match. This catches (a) data drift, (b) golden
corruption, and (c) constants regressions. Once the Rust CLI exists, the same golden files
become its oracle (the compare semantics are recorded inside each golden json).

Exit code 0 = all checks pass; nonzero = at least one failure (CI-friendly).
"""
import json, os, sys, math, struct, importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
GOLDEN = os.path.join(ROOT, "golden")
DATA = os.path.join(ROOT, "data")

PASS, FAIL = 0, 0
def check(name, ok, detail=""):
    global PASS, FAIL
    if ok: PASS += 1
    else:  FAIL += 1; print(f"  FAIL: {name}  {detail}")

def load(p): return json.load(open(p))
def close(a, b, tol): return abs(a - b) <= tol

# import the reference parsers so we re-read raw spectra the same way the golden was built
def _import(mod, path):
    spec = importlib.util.spec_from_file_location(mod, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
spectra_mod = _import("make_spectra_golden", os.path.join(ROOT, "reference", "make_spectra_golden.py"))


def test_chemistry():
    cdir = os.path.join(GOLDEN, "chemistry")
    if not os.path.isdir(cdir): print("  (chemistry golden absent — run reference/make_chemistry_golden.py)"); return
    C = load(os.path.join(cdir, "constants.golden.json"))
    atoms, proton, H2O = C["atoms"], C["proton"], C["H2O"]
    check("constants.H2O", close(H2O, 2*atoms["H"] + atoms["O"], 1e-9))

    R = load(os.path.join(cdir, "residue_masses.golden.json"))["residues"]
    for aa, rec in R.items():
        m = sum(atoms[el]*n for el, n in rec["formula"].items())
        check(f"residue[{aa}]", close(m, rec["mass"], 1e-6), f"{m} vs {rec['mass']}")

    P = load(os.path.join(cdir, "peptide_masses.golden.json"))["peptides"]
    for pep in P:
        seq = pep["sequence"]
        neutral = sum(R[a]["mass"] for a in seq) + H2O
        check(f"peptide[{seq}].neutral", close(neutral, pep["neutral_mass"], 1e-4))
        for z, key in [(1,"1+"),(2,"2+"),(3,"3+")]:
            check(f"peptide[{seq}].{key}", close((neutral+z*proton)/z, pep["mz"][key], 1e-4))
        if "published_mh" in pep:
            check(f"calibrant[{seq}] vs published", close(pep["mz"]["1+"], pep["published_mh"], 0.003),
                  f"{pep['mz']['1+']:.4f} vs {pep['published_mh']}")

    F = load(os.path.join(cdir, "fragment_ions.golden.json"))["peptides"]
    for fp in F:
        seq = fp["peptide"]
        pre = 0.0
        for b in fp["b_ions"]:
            pre += R[seq[b["i"]-1]]["mass"]
            check(f"b[{seq},{b['i']}].z1", close(pre+proton, b["z1"], 1e-4))
        suf, n = 0.0, len(seq)
        for y in fp["y_ions"]:
            suf += R[seq[n-y["j"]]]["mass"]
            check(f"y[{seq},{y['j']}].z1", close(suf+H2O+proton, y["z1"], 1e-4))


def test_spectra():
    sdir = os.path.join(GOLDEN, "spectra")
    if not os.path.isdir(sdir): return
    parsers = {"mgf": spectra_mod.parse_mgf, "mzML": spectra_mod.parse_mzml}
    for fn in sorted(os.listdir(sdir)):
        g = load(os.path.join(sdir, fn))
        path = os.path.join(DATA, "spectra", g["file"])
        if not os.path.exists(path): print(f"  (skip {g['file']} — data absent)"); continue
        stored = {s["index"]: s for s in g["spectra"]}
        count = total = 0
        import hashlib; roll = hashlib.sha1()
        for rec in parsers[g["format"]](path):
            count += 1; total += rec["n_peaks"]; roll.update(rec["peaks_sha1"].encode())
            if rec["index"] in stored:
                s = stored[rec["index"]]
                ok = (rec["n_peaks"]==s["n_peaks"] and rec["charge"]==s["charge"]
                      and rec["peaks_sha1"]==s["peaks_sha1"])
                check(f"{g['file']}#{rec['index']}", ok)
        check(f"{g['file']}.n_spectra", count==g["n_spectra"], f"{count} vs {g['n_spectra']}")
        check(f"{g['file']}.total_peaks", total==g["total_peaks"])
        check(f"{g['file']}.rolling_hash", roll.hexdigest()==g["rolling_peak_sha1"])


def test_param_inventory():
    p = os.path.join(GOLDEN, "models", "param_inventory.golden.json")
    if not os.path.exists(p): return
    import hashlib
    for m in load(p)["models"]:
        fp = os.path.join(DATA, "models", m["file"])
        if not os.path.exists(fp): print(f"  (skip {m['file']} — absent)"); continue
        b = open(fp, "rb").read()
        check(f"{m['file']}.size", len(b)==m["size"])
        check(f"{m['file']}.sha256", hashlib.sha256(b).hexdigest()==m["sha256"])
        off = 4; enums = []
        for _ in range(4):
            n = b[off]; off += 1; enums.append(b[off:off+2*n].decode("utf-16-be")); off += 2*n
        idn = m["identity"]
        check(f"{m['file']}.identity",
              enums[0]==idn["activation"] and enums[2]==idn["enzyme"])


def test_worked_example():
    p = os.path.join(GOLDEN, "worked_example.golden.json")
    if not os.path.exists(p): return
    g = load(p)
    check("worked_example.n_psms", g["n_psms"] == len(g["psms"]))
    for psm in g["psms"]:
        check(f"worked_example[{psm['peptide'][:12]}].bounds",
              psm["raw_score"] <= psm["denovo_score"] and 0 < psm["spec_evalue"] <= 1)
    check("worked_example.input_present", os.path.exists(os.path.join(DATA, "spectra", "test.mgf")),
          "test.mgf missing (fetch_reference_data.sh --full)")


def test_msgf_golden():
    """Authoritative MS-GF+ search goldens (e.g. iprg2013_F13) if generated."""
    for fn in sorted(os.listdir(GOLDEN)):
        if not fn.endswith(".golden.json") or fn == "worked_example.golden.json":
            continue
        g = load(os.path.join(GOLDEN, fn))
        if "psms" not in g: continue
        n_bad = 0
        for psm in g["psms"]:
            if not (psm["raw_score"] <= psm["denovo_score"] and 0 < psm["spec_evalue"] <= 1
                    and psm["evalue"] > 0):
                n_bad += 1
        check(f"{fn}.psm_invariants", n_bad == 0, f"{n_bad}/{len(g['psms'])} bad")
        check(f"{fn}.nonempty", g["n_psms"] > 0)
        print(f"  [{fn}: {g['n_psms']} PSMs]")


def _f32(x):
    """Round an f64 to the nearest f32 — the Java oracle computes q-values in `float`."""
    return struct.unpack("f", struct.pack("f", x))[0]


def _fdr_map(targets, decoys, pit=1.0):
    """Transcription of edu.ucsd.msjava.fdr.TargetDecoyAnalysis.getFDRMap, smaller-is-better.

    A run of equal decoy scores is charged at its first index; a threshold with no better target
    is skipped without stopping the sweep; both infinite sentinels are seeded (the upper one at 0
    when there are no decoys); FDRs become q-values by a running minimum from worst to best.
    """
    t, d = sorted(targets), sorted(decoys)
    swept, ti, prev = [], 0, None
    for di, key in enumerate(d):
        if key == prev: continue
        prev = key
        while ti < len(t) and t[ti] < key: ti += 1
        if ti == 0: continue                     # no entry, but the sweep carries on
        fdr = 1.0 if ti <= di else min(_f32(math.floor(di * pit + 0.5) / ti), 1.0)
        swept.append((key, fdr))
        if fdr >= 1.0: break
    pairs = [(-math.inf, 0.0)] + swept + [(math.inf, 1.0 if d else 0.0)]
    run = 1.0
    for i in range(len(pairs) - 1, -1, -1):
        run = min(run, pairs[i][1]); pairs[i] = (pairs[i][0], run)
    return pairs


def _q_value(pairs, score):
    """Java TreeMap.higherEntry: the least threshold strictly greater than `score`."""
    for k, q in pairs:
        if k > score: return q
    return pairs[-1][1]


def test_fdr_map():
    """Re-derive the Java-dumped target-decoy cases (PLAN2 TD-2 Gate 2) without a JVM."""
    p = os.path.join(GOLDEN, "fdr", "fdrmap_cases.golden.json")
    if not os.path.exists(p):
        print("  (fdr golden absent — run reference/make_fdr_golden.sh, JVM only)"); return
    g = load(p)
    for case in g["cases"]:
        name = case["name"]
        targets = [float(s) for s in case["targets"]]
        decoys = [float(s) for s in case["decoys"]]
        pairs = _fdr_map(targets, decoys)
        want = [(float(e["key"]), _f32(float(e["q"]))) for e in case["map"]]
        check(f"fdrmap[{name}].map", pairs == want, f"{pairs} vs {want}")
        n_bad = sum(1 for l in case["lookups"] if l["q"] is not None
                    and _q_value(pairs, float(l["score"])) != _f32(float(l["q"])))
        check(f"fdrmap[{name}].lookups", n_bad == 0, f"{n_bad}/{len(case['lookups'])} differ")
        qs = [q for _, q in pairs]
        check(f"fdrmap[{name}].monotone", all(a <= b for a, b in zip(qs, qs[1:])))
    print(f"  [{len(g['cases'])} target-decoy cases re-derived from MS-GF+'s inputs]")


def main():
    print("== chemistry =="); test_chemistry()
    print("== spectra =="); test_spectra()
    print("== param inventory =="); test_param_inventory()
    print("== worked example =="); test_worked_example()
    print("== msgf search goldens =="); test_msgf_golden()
    print("== target-decoy fdr =="); test_fdr_map()
    print(f"\n{'OK' if FAIL==0 else 'FAILED'}: {PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
