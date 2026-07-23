#!/usr/bin/env python3
"""
run_regression.py — the MSGF_Rust regression suite.

Runs TODAY with no Java and no Rust: it re-derives every golden fixture from the raw data and
authoritative constants and asserts they match. This catches (a) data drift, (b) golden
corruption, and (c) constants regressions. Once the Rust CLI exists, the same golden files
become its oracle (the compare semantics are recorded inside each golden json).

Exit code 0 = all checks pass; nonzero = at least one failure (CI-friendly).
"""
import json, os, sys, math, importlib.util

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


def main():
    print("== chemistry =="); test_chemistry()
    print("== spectra =="); test_spectra()
    print("== param inventory =="); test_param_inventory()
    print("== worked example =="); test_worked_example()
    print("== msgf search goldens =="); test_msgf_golden()
    print(f"\n{'OK' if FAIL==0 else 'FAILED'}: {PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
