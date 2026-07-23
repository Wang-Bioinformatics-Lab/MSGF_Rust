#!/usr/bin/env python3
"""
parse_msgf_tsv.py — normalize an MS-GF+ TSV result (from MzIDToTsv) into frozen golden JSON.

The golden JSON is the regression oracle the Rust implementation is checked against. We keep
only the fields we assert on, with their comparison semantics recorded alongside:

  RawScore     (MS:1002049, TSV "MSGFScore")   int   -> assert EXACT
  DeNovoScore  (MS:1002050)                     int   -> assert EXACT
  SpecEValue   (MS:1002052)                     float -> assert |log10 ratio| <= tol
  EValue       (MS:1002053)                     float -> assert |log10 ratio| <= tol

Each PSM is keyed by (SpecFile, ScanNum, Charge, Peptide) so Rust output can be joined to it.

Usage:
  parse_msgf_tsv.py <result.tsv> -o <golden.json> [--tol-log10 0.05]
"""
import argparse, json, sys, os

# TSV column name -> (golden field, kind)
FIELDS = {
    "MSGFScore":    ("raw_score",    "int"),
    "DeNovoScore":  ("denovo_score", "int"),
    "SpecEValue":   ("spec_evalue",  "float"),
    "EValue":       ("evalue",       "float"),
    "QValue":       ("qvalue",       "float"),
    "PepQValue":    ("pep_qvalue",   "float"),
}
KEYS = ["#SpecFile", "ScanNum", "Charge", "Peptide"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tsv")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--tol-log10", type=float, default=0.05,
                    help="allowed |log10(rust/java)| for float fields")
    args = ap.parse_args()

    with open(args.tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for req in KEYS + list(FIELDS):
            if req not in idx:
                sys.exit(f"missing expected column {req!r} in {args.tsv}\nheader: {header}")
        psms = []
        for line in fh:
            if not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            rec = {
                "spec_file": os.path.basename(c[idx["#SpecFile"]]),
                "scan":      c[idx["ScanNum"]],
                "charge":    int(c[idx["Charge"]]),
                "peptide":   c[idx["Peptide"]],
                "protein":   c[idx["Protein"]] if "Protein" in idx else None,
            }
            for col, (field, kind) in FIELDS.items():
                v = c[idx[col]]
                rec[field] = int(v) if kind == "int" else float(v)
            psms.append(rec)

    out = {
        "source_tsv": os.path.basename(args.tsv),
        "generator": "MS-GF+ (Java) via MzIDToTsv; see generate_golden.sh",
        "n_psms": len(psms),
        "compare": {
            "raw_score":    {"kind": "int",   "assert": "exact"},
            "denovo_score": {"kind": "int",   "assert": "exact"},
            "spec_evalue":  {"kind": "float", "assert": "abs_log10_ratio", "tol": args.tol_log10},
            "evalue":       {"kind": "float", "assert": "abs_log10_ratio", "tol": args.tol_log10},
        },
        "key_fields": ["spec_file", "scan", "charge", "peptide"],
        "psms": psms,
    }
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"wrote {args.out}: {len(psms)} PSMs")


if __name__ == "__main__":
    main()
