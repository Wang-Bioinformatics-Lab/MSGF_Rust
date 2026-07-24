#!/usr/bin/env bash
#
# make_highres_goldens.sh -- regenerate the two RawScore validation goldens on
# the HighRes MS-GF+ model (HCD_HighRes_Tryp.param), the model the iprg2013_F13
# search actually used (-inst 1 = HighRes / Orbitrap, NOT -inst 3 = QExactive).
#
# It reuses the EXACT 30-spectrum selection + golden peptides already frozen in
# golden/rawscore/f13_scored_spectrum.golden.json (the QExactive scored-spectrum
# golden): scan, charge, golden_peptide (with flanking context) and
# golden_raw_score. Only the scoring MODEL changes -- everything else is held
# fixed -- so the QExactive and HighRes goldens are directly comparable.
#
# 1. ScoredSpectrumDumper_HighRes.java runs MS-GF+ NewRankScorer + NewScoredSpectrum
#    with the HighRes model and dumps, per spectrum, the preprocessed peak list
#    plus per-nominal-mass prefix/suffix node scores.
#      -> golden/rawscore/f13_scored_spectrum_highres.golden.json
# 2. RawScoreDumper_HighRes.java reconstructs MS-GF+'s peptide prefix-mass arrays
#    exactly as CandidatePeptideGrid + DBScanner do, then runs the reference
#    FastScorer (node summation) and DBScanScorer (node + edge) with the HighRes
#    scored spectra and dumps node_only / full / edge / cleavage scores.
#      -> golden/rawscore/f13_rawscore_highres.golden.json
#
# The two *_HighRes.java drivers are verbatim copies of ScoredSpectrumDumper.java
# and RawScoreDumper.java; the ONLY difference is the MODEL_NAME label constant
# (set to "HCD_HighRes_Tryp.param") so the emitted JSON "model" field honestly
# names the HighRes model. Both original drivers already take the model path as
# args[0], so the loaded model is identical to what these copies load.
#
# Every number comes from actually running MSGFPlus.jar; nothing is fabricated.
#
# JVM-only; needs reference/MSGFPlus.jar, data/models/HCD_HighRes_Tryp.param,
# data/spectra/F13.mgf, data/config/iprg-2013_Mods.txt and
# golden/rawscore/f13_scored_spectrum.golden.json (source of the 30 PSMs).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JAR="$HERE/MSGFPlus.jar"
# The iprg2013_F13 search ran `-inst 1` = HighRes, which loads HCD_HighRes_Tryp.param.
MODEL="$HERE/../data/models/HCD_HighRes_Tryp.param"
MGF="$HERE/../data/spectra/F13.mgf"
MODS="$HERE/../data/config/iprg-2013_Mods.txt"
SCORED="$HERE/../golden/rawscore/f13_scored_spectrum.golden.json"
OUTDIR="$HERE/../golden/rawscore"
OUT_SS="$OUTDIR/f13_scored_spectrum_highres.golden.json"
OUT_RS="$OUTDIR/f13_rawscore_highres.golden.json"

# Use a local javac/java if present, else the msgfjava conda env.
command -v javac >/dev/null 2>&1 && JVM() { "$@"; } || JVM() { conda run -n msgfjava "$@"; }
[[ -f "$JAR" ]]    || { echo "ERROR: $JAR missing (fetch_reference_data.sh --jar)"; exit 1; }
[[ -f "$MODEL" ]]  || { echo "ERROR: $MODEL missing"; exit 1; }
[[ -f "$SCORED" ]] || { echo "ERROR: $SCORED missing (run make_scored_spectrum_golden.sh first)"; exit 1; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
mkdir -p "$OUTDIR"

# Selection TSV (scan<TAB>charge<TAB>golden_peptide<TAB>golden_raw_score), in the
# order the PSMs appear in the (QExactive) scored-spectrum golden. This is the
# EXACT 30-spectrum selection + golden peptides both HighRes drivers consume.
python3 - "$SCORED" "$WORK/selection.tsv" <<'PY'
import json, sys
scored, out = sys.argv[1], sys.argv[2]
d = json.load(open(scored))
with open(out, "w") as fh:
    for s in d["spectra"]:
        fh.write(f"{s['scan']}\t{s['charge']}\t{s['golden_peptide']}\t{s['golden_raw_score']}\n")
print(f"wrote {len(d['spectra'])} selection rows")
PY

# Compile both HighRes drivers into the same work dir.
JVM javac -cp "$JAR" -d "$WORK" \
    "$HERE/java/ScoredSpectrumDumper_HighRes.java" \
    "$HERE/java/RawScoreDumper_HighRes.java"

# 1. HighRes scored-spectrum golden (preprocessed peaks + prefix/suffix scores).
JVM java -cp "$WORK:$JAR" ScoredSpectrumDumper_HighRes "$MODEL" "$MGF" "$WORK/selection.tsv" "$OUT_SS"

# 2. HighRes RawScore golden (node_only / full / edge / cleavage).
JVM java -cp "$WORK:$JAR" RawScoreDumper_HighRes "$MODEL" "$MGF" "$MODS" "$WORK/selection.tsv" "$OUT_RS"

echo "wrote $OUT_SS ($(wc -c < "$OUT_SS") bytes)"
echo "wrote $OUT_RS ($(wc -c < "$OUT_RS") bytes)"
