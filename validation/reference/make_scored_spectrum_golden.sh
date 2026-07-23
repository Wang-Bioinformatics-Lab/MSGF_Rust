#!/usr/bin/env bash
#
# make_scored_spectrum_golden.sh — build the RawScore validation oracle.
#
# 1. select_scored_spectra.py picks ~30 distinct F13 scans (spread of charge 2/3,
#    higher raw_score preferred) from the golden PSMs, cross-checked against F13.mgf.
# 2. ScoredSpectrumDumper.java (compiled against MSGFPlus.jar) runs MS-GF+
#    NewRankScorer + NewScoredSpectrum for each and dumps the preprocessed peak
#    list plus per-nominal-mass prefix/suffix node scores.
# Output: golden/rawscore/f13_scored_spectrum.golden.json
#
# JVM-only; needs reference/MSGFPlus.jar, data/models/HCD_QExactive_Tryp.param,
# data/spectra/F13.mgf and golden/iprg2013_F13.golden.json.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JAR="$HERE/MSGFPlus.jar"
MODEL="$HERE/../data/models/HCD_QExactive_Tryp.param"
MGF="$HERE/../data/spectra/F13.mgf"
PSMS="$HERE/../golden/iprg2013_F13.golden.json"
OUTDIR="$HERE/../golden/rawscore"
OUT="$OUTDIR/f13_scored_spectrum.golden.json"
N_PER_CHARGE="${1:-15}"

# Use a local javac/java if present, else the msgfjava conda env.
command -v javac >/dev/null 2>&1 && JVM() { "$@"; } || JVM() { conda run -n msgfjava "$@"; }
[[ -f "$JAR" ]] || { echo "ERROR: $JAR missing (fetch_reference_data.sh --jar)"; exit 1; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
mkdir -p "$OUTDIR"

python3 "$HERE/select_scored_spectra.py" "$PSMS" "$MGF" "$WORK/selection.tsv" "$N_PER_CHARGE"

JVM javac -cp "$JAR" -d "$WORK" "$HERE/java/ScoredSpectrumDumper.java"
JVM java -cp "$WORK:$JAR" ScoredSpectrumDumper "$MODEL" "$MGF" "$WORK/selection.tsv" "$OUT"

echo "wrote $OUT ($(wc -c < "$OUT") bytes)"
