#!/usr/bin/env bash
#
# make_model_goldens.sh — dump each MS-GF+ .param model to authoritative plain text (via
# NewRankScorer.writeParametersPlainText through jshell) and distill it into a structured golden
# for the Rust .param reader. JVM-only; needs reference/MSGFPlus.jar and the .param data.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JAR="$HERE/MSGFPlus.jar"
DATA="$HERE/../data/models"
GOLD="$HERE/../golden/models"

command -v jshell >/dev/null 2>&1 || JSHELL_VIA_CONDA=1
run_jshell() { if [[ "${JSHELL_VIA_CONDA:-0}" == 1 ]]; then conda run -n msgfjava jshell "$@"; else jshell "$@"; fi; }

[[ -f "$JAR" ]] || { echo "ERROR: $JAR missing (fetch_reference_data.sh --jar)"; exit 1; }
mkdir -p "$GOLD"

for param in "$DATA"/*.param; do
  [[ -e "$param" ]] || { echo "no .param models in $DATA (fetch_reference_data.sh)"; exit 1; }
  base="$(basename "$param" .param)"
  txt="$GOLD/$base.model.txt"
  printf 'var s = new edu.ucsd.msjava.msscorer.NewRankScorer();\ns.readFromFile(new java.io.File("%s"));\ns.writeParametersPlainText(new java.io.File("%s"));\n/exit\n' \
    "$param" "$txt" > "/tmp/dump_$base.jsh"
  run_jshell --class-path "$JAR" "/tmp/dump_$base.jsh" >/dev/null 2>&1
  python3 "$HERE/make_model_golden.py" "$txt" "$param" -o "$GOLD/$base.model.golden.json"
done
