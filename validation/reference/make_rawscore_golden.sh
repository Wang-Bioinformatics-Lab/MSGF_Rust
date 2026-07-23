#!/usr/bin/env bash
#
# make_rawscore_golden.sh -- build the FINAL RawScore validation oracle.
#
# For the 30 (scan, charge, golden_peptide, golden_raw_score) PSMs already frozen
# in golden/rawscore/f13_scored_spectrum.golden.json, RawScoreDumper.java (compiled
# against MSGFPlus.jar) reconstructs MS-GF+'s peptide prefix-mass arrays exactly as
# edu.ucsd.msjava.msdbsearch's CandidatePeptideGrid + DBScanner do, then runs the
# reference FastScorer (node summation) and DBScanScorer (node + edge) and dumps
# node_only / full / edge / cleavage scores so the Rust RawScore summation + edge
# scoring can be validated. Every number comes from running MS-GF+.
#
# Output: golden/rawscore/f13_rawscore.golden.json
#
# JVM-only; needs reference/MSGFPlus.jar, data/models/HCD_QExactive_Tryp.param,
# data/spectra/F13.mgf, data/config/iprg-2013_Mods.txt and
# golden/rawscore/f13_scored_spectrum.golden.json (source of the 30 PSMs).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JAR="$HERE/MSGFPlus.jar"
MODEL="$HERE/../data/models/HCD_QExactive_Tryp.param"
MGF="$HERE/../data/spectra/F13.mgf"
MODS="$HERE/../data/config/iprg-2013_Mods.txt"
SCORED="$HERE/../golden/rawscore/f13_scored_spectrum.golden.json"
OUTDIR="$HERE/../golden/rawscore"
OUT="$OUTDIR/f13_rawscore.golden.json"

# Use a local javac/java if present, else the msgfjava conda env.
command -v javac >/dev/null 2>&1 && JVM() { "$@"; } || JVM() { conda run -n msgfjava "$@"; }
[[ -f "$JAR" ]] || { echo "ERROR: $JAR missing (fetch_reference_data.sh --jar)"; exit 1; }
[[ -f "$SCORED" ]] || { echo "ERROR: $SCORED missing (run make_scored_spectrum_golden.sh first)"; exit 1; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
mkdir -p "$OUTDIR"

# Selection TSV (scan<TAB>charge<TAB>golden_peptide<TAB>golden_raw_score), in the
# order the PSMs appear in the scored-spectrum golden.
python3 - "$SCORED" "$WORK/selection.tsv" <<'PY'
import json, sys
scored, out = sys.argv[1], sys.argv[2]
d = json.load(open(scored))
with open(out, "w") as fh:
    for s in d["spectra"]:
        fh.write(f"{s['scan']}\t{s['charge']}\t{s['golden_peptide']}\t{s['golden_raw_score']}\n")
print(f"wrote {len(d['spectra'])} selection rows")
PY

JVM javac -cp "$JAR" -d "$WORK" "$HERE/java/RawScoreDumper.java"
JVM java -cp "$WORK:$JAR" RawScoreDumper "$MODEL" "$MGF" "$MODS" "$WORK/selection.tsv" "$OUT"

echo "wrote $OUT ($(wc -c < "$OUT") bytes)"
