#!/usr/bin/env python3
"""
select_scored_spectra.py — pick the spectra whose scored-spectrum data the Java
ScoredSpectrumDumper should dump for the RawScore validation oracle.

Reads the F13 golden PSMs, dedupes to distinct scans (keeping the best raw_score
per scan), cross-checks each scan against the actual F13.mgf (scan must exist and
its MGF precursor charge must equal the golden charge), then selects a spread of
charge-2 and charge-3 spectra preferring higher raw_score.

Emits a TSV (scan<TAB>charge<TAB>peptide<TAB>raw_score) consumed by the Java
driver. Selection lives here (Python reads the big JSON trivially); all MS-GF+
numbers come from actually running the Java scorer, never from this file.

Usage: select_scored_spectra.py <golden.json> <F13.mgf> <out.tsv> [n_per_charge]
"""
import json
import sys


def mgf_scan_charge(mgf_path):
    """Map MGF scan number -> precursor charge (int)."""
    scan2charge = {}
    scan = None
    charge = None
    with open(mgf_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("SCANS="):
                scan = int(line[len("SCANS="):])
            elif line.startswith("CHARGE="):
                # e.g. "2+" -> 2
                charge = int(line[len("CHARGE="):].rstrip("+").strip())
            elif line == "END IONS":
                if scan is not None and charge is not None:
                    scan2charge[scan] = charge
                scan = charge = None
    return scan2charge


def main():
    golden_json, mgf_path, out_tsv = sys.argv[1], sys.argv[2], sys.argv[3]
    n_per_charge = int(sys.argv[4]) if len(sys.argv) > 4 else 15

    scan2charge = mgf_scan_charge(mgf_path)
    psms = json.load(open(golden_json))["psms"]

    # Dedupe to distinct scans, keeping the PSM with the max raw_score per scan.
    best = {}
    for p in psms:
        scan = int(p["scan"])
        rs = int(p["raw_score"])
        if scan not in best or rs > best[scan]["raw_score"]:
            best[scan] = {"scan": scan, "charge": int(p["charge"]),
                          "peptide": p["peptide"], "raw_score": rs}

    # Keep only scans present in the MGF whose MGF charge matches the golden charge.
    ok = [b for b in best.values()
          if b["scan"] in scan2charge and scan2charge[b["scan"]] == b["charge"]]

    picked = []
    for target_charge in (2, 3):
        pool = sorted((b for b in ok if b["charge"] == target_charge),
                      key=lambda b: (-b["raw_score"], b["scan"]))
        picked.extend(pool[:n_per_charge])

    picked.sort(key=lambda b: b["scan"])  # stable, deterministic order
    with open(out_tsv, "w") as fh:
        for b in picked:
            fh.write(f"{b['scan']}\t{b['charge']}\t{b['peptide']}\t{b['raw_score']}\n")

    n2 = sum(1 for b in picked if b["charge"] == 2)
    n3 = sum(1 for b in picked if b["charge"] == 3)
    print(f"selected {len(picked)} scans ({n2} charge-2, {n3} charge-3) -> {out_tsv}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
