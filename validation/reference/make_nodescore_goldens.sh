#!/usr/bin/env bash
#
# make_nodescore_goldens.sh — compile the ScoreDumper Java driver against the MS-GF+ jar, run it
# for every .param model to emit authoritative getNodeScore/getMissingIonScore values, and
# distill them into golden/models/node_scores.golden.json for the Rust scoring primitives.
# JVM-only; needs reference/MSGFPlus.jar and the .param data.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JAR="$HERE/MSGFPlus.jar"
DATA="$HERE/../data/models"
GOLD="$HERE/../golden/models"
CLASSES="$(mktemp -d)"

command -v javac >/dev/null 2>&1 && JVM() { "$@"; } || JVM() { conda run -n msgfjava "$@"; }
[[ -f "$JAR" ]] || { echo "ERROR: $JAR missing (fetch_reference_data.sh --jar)"; exit 1; }

JVM javac -cp "$JAR" -d "$CLASSES" "$HERE/java/ScoreDumper.java"

tsvs=()
for param in "$DATA"/*.param; do
  [[ -e "$param" ]] || { echo "no .param models in $DATA"; exit 1; }
  base="$(basename "$param" .param)"
  tsv="$CLASSES/$base.nodescore.tsv"
  JVM java -cp "$CLASSES:$JAR" ScoreDumper "$param" "$tsv"
  tsvs+=("$tsv")
done

python3 "$HERE/make_nodescore_golden.py" "${tsvs[@]}" -o "$GOLD/node_scores.golden.json"
rm -rf "$CLASSES"
