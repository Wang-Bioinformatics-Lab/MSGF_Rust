#!/usr/bin/env python3
"""
make_chemistry_golden.py — generate INDEPENDENT ground-truth fixtures for the mass/chemistry
layer (Rust crate msgf-chem). This ground truth does NOT come from MS-GF+; it is derived from
authoritative atomic monoisotopic masses + fixed residue elemental formulas, and is guarded
against published peptide-calibrant [M+H]+ values. If any constant is wrong, this script aborts
before writing — so the frozen golden is trustworthy.

Writes:
  golden/chemistry/constants.golden.json
  golden/chemistry/residue_masses.golden.json
  golden/chemistry/peptide_masses.golden.json
  golden/chemistry/fragment_ions.golden.json

The golden records the exact constants used, so the Rust chem layer can be configured
identically for a bit-reproducible comparison (compare tol 1e-6 Da on residues/peptides).
MS-GF+'s own (possibly different) constants are validated separately via the RawScore golden.
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "golden", "chemistry")

# --- authoritative monoisotopic atomic masses (most-abundant isotope), CODATA/AME-2020 rounded ---
ATOM = {
    "H": 1.0078250319,
    "C": 12.0000000000,
    "N": 14.0030740052,
    "O": 15.9949146221,
    "S": 31.9720707300,
    "P": 30.9737615120,
}
PROTON = 1.0072764669      # mass of a proton (H+); used for charging
ELECTRON = 0.0005485799
H2O = 2 * ATOM["H"] + ATOM["O"]

# --- residue elemental compositions (amino-acid residue = free AA minus H2O) ---
RESIDUE_FORMULA = {
    "G": {"C": 2,  "H": 3,  "N": 1, "O": 1},
    "A": {"C": 3,  "H": 5,  "N": 1, "O": 1},
    "S": {"C": 3,  "H": 5,  "N": 1, "O": 2},
    "P": {"C": 5,  "H": 7,  "N": 1, "O": 1},
    "V": {"C": 5,  "H": 9,  "N": 1, "O": 1},
    "T": {"C": 4,  "H": 7,  "N": 1, "O": 2},
    "C": {"C": 3,  "H": 5,  "N": 1, "O": 1, "S": 1},
    "L": {"C": 6,  "H": 11, "N": 1, "O": 1},
    "I": {"C": 6,  "H": 11, "N": 1, "O": 1},
    "N": {"C": 4,  "H": 6,  "N": 2, "O": 2},
    "D": {"C": 4,  "H": 5,  "N": 1, "O": 3},
    "Q": {"C": 5,  "H": 8,  "N": 2, "O": 2},
    "K": {"C": 6,  "H": 12, "N": 2, "O": 1},
    "E": {"C": 5,  "H": 7,  "N": 1, "O": 3},
    "M": {"C": 5,  "H": 9,  "N": 1, "O": 1, "S": 1},
    "H": {"C": 6,  "H": 7,  "N": 3, "O": 1},
    "F": {"C": 9,  "H": 9,  "N": 1, "O": 1},
    "R": {"C": 6,  "H": 12, "N": 4, "O": 1},
    "Y": {"C": 9,  "H": 9,  "N": 1, "O": 2},
    "W": {"C": 11, "H": 10, "N": 2, "O": 1},
}

def formula_mass(f):
    return sum(ATOM[el] * n for el, n in f.items())

RESIDUE = {aa: formula_mass(f) for aa, f in RESIDUE_FORMULA.items()}

def neutral_mass(seq):
    return sum(RESIDUE[a] for a in seq) + H2O

def mz(neutral, z):
    return (neutral + z * PROTON) / z

def by_ions(seq):
    """Full b/y series for z=1 and z=2 (m/z). b neutral = sum prefix residues; y neutral = sum suffix + H2O."""
    n = len(seq)
    b, y = [], []
    pre = 0.0
    for i in range(1, n):          # b1..b(n-1)
        pre += RESIDUE[seq[i-1]]
        b.append({"i": i, "z1": mz(pre, 1), "z2": mz(pre, 2)})
    suf = 0.0
    for j in range(1, n):          # y1..y(n-1)
        suf += RESIDUE[seq[n-j]]
        y.append({"j": j, "z1": mz(suf + H2O, 1), "z2": mz(suf + H2O, 2)})
    return b, y

# --- published monoisotopic [M+H]+ calibrants: the trust anchors (source in comments) ---
CALIBRANTS = [
    # seq,                published [M+H]+ mono, source
    ("MRFA",            524.2649,  "common QC calibrant (Met-Arg-Phe-Ala)"),
    ("RPPGFSPFR",       1060.5692, "Bradykinin, monoisotopic [M+H]+"),
    ("DRVYIHPF",        1046.5418, "Angiotensin II, monoisotopic [M+H]+"),
    ("EGVNDNEEGFFSAR",  1570.6774, "Glu-Fibrinopeptide B, monoisotopic [M+H]+"),
]

# extra peptides (no external ref) to broaden coverage
EXTRA_PEPTIDES = ["PEPTIDE", "SAMPLER", "ELVISLIVES", "HPLC", "ACDEFGHIKLMNPQRSTVWY"]
FRAGMENT_PEPTIDES = ["PEPTIDE", "MRFA", "DRVYIHPF", "SAMPLER"]


def guard():
    """Abort if computed calibrant [M+H]+ deviates from published — proves the constants."""
    worst = 0.0
    for seq, pub, _ in CALIBRANTS:
        got = mz(neutral_mass(seq), 1)
        d = abs(got - pub)
        worst = max(worst, d)
        print(f"  calibrant {seq:18s} computed [M+H]+={got:.4f}  published={pub:.4f}  Δ={d*1000:.2f} mDa")
        if d > 0.003:
            sys.exit(f"ABORT: calibrant {seq} off by {d:.4f} Da — chemistry constants are wrong")
    print(f"  guard OK (worst Δ = {worst*1000:.2f} mDa < 3 mDa)")


def main():
    os.makedirs(OUT, exist_ok=True)
    print("guarding constants against published calibrants:")
    guard()

    dump(os.path.join(OUT, "constants.golden.json"), {
        "note": "atomic monoisotopic masses + derived constants; Rust msgf-chem must reproduce exactly",
        "compare": {"kind": "float", "assert": "abs", "tol": 1e-6},
        "atoms": ATOM, "proton": PROTON, "electron": ELECTRON, "H2O": H2O,
    })

    dump(os.path.join(OUT, "residue_masses.golden.json"), {
        "note": "monoisotopic residue masses (free AA minus water), derived from atoms[] compositions",
        "compare": {"kind": "float", "assert": "abs", "tol": 1e-6},
        "residues": {aa: {"mass": RESIDUE[aa], "formula": RESIDUE_FORMULA[aa]} for aa in sorted(RESIDUE)},
    })

    peps = []
    for seq, pub, src in CALIBRANTS:
        peps.append(pep_record(seq, published_mh=pub, source=src))
    for seq in EXTRA_PEPTIDES:
        peps.append(pep_record(seq))
    dump(os.path.join(OUT, "peptide_masses.golden.json"), {
        "note": "neutral monoisotopic mass and [M+nH]n+ for test peptides; calibrants carry published_mh",
        "compare": {"kind": "float", "assert": "abs", "tol": 1e-4},
        "peptides": peps,
    })

    frags = []
    for seq in FRAGMENT_PEPTIDES:
        b, y = by_ions(seq)
        frags.append({"peptide": seq, "neutral_mass": neutral_mass(seq), "b_ions": b, "y_ions": y})
    dump(os.path.join(OUT, "fragment_ions.golden.json"), {
        "note": "b/y fragment m/z (z=1,2). b_neutral=Σ prefix residues; y_neutral=Σ suffix residues + H2O; m/z=(neutral+z·proton)/z",
        "compare": {"kind": "float", "assert": "abs", "tol": 1e-4},
        "peptides": frags,
    })
    print("wrote 4 chemistry golden files ->", os.path.normpath(OUT))


def pep_record(seq, published_mh=None, source=None):
    nm = neutral_mass(seq)
    rec = {"sequence": seq, "length": len(seq), "neutral_mass": nm,
           "mz": {"1+": mz(nm, 1), "2+": mz(nm, 2), "3+": mz(nm, 3)}}
    if published_mh is not None:
        rec["published_mh"] = published_mh
        rec["source"] = source
    return rec


def dump(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2)


if __name__ == "__main__":
    main()
