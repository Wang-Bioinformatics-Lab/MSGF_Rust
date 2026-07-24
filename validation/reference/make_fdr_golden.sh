#!/usr/bin/env bash
#
# make_fdr_golden.sh -- build the TARGET-DECOY FDR oracle (PLAN2.md TD-2 "Gate 2").
#
# The F13 search golden already pins MS-GF+'s QValue/PepQValue columns end-to-end, but it only
# produces two distinct q-values (PLAN2.md section 4), so it cannot separate implementations that
# differ in tie handling, in the `targetIndex > 0` guard, or in the map-lookup rule. DumpFdrMap.java
# drives edu.ucsd.msjava.fdr.TargetDecoyAnalysis directly on small synthetic score lists chosen to
# separate exactly those behaviours, and freezes the map plus the lookups (including each key's
# immediate float neighbours).
#
# Output: golden/fdr/fdrmap_cases.golden.json
#
# JVM-only; needs reference/MSGFPlus.jar (fetch_reference_data.sh --jar) and nothing else -- no
# spectra, no models, no database. The output is committed, so the Rust test needs no fetched data.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JAR="$HERE/MSGFPlus.jar"
OUTDIR="$HERE/../golden/fdr"
OUT="$OUTDIR/fdrmap_cases.golden.json"

# Use $JAVA_HOME/bin, else a local javac/java, else the msgfjava conda env.
if [[ -n "${JAVA_HOME:-}" && -x "$JAVA_HOME/bin/javac" ]]; then
  JVM() { "$JAVA_HOME/bin/$1" "${@:2}"; }
elif command -v javac >/dev/null 2>&1; then
  JVM() { "$@"; }
else
  JVM() { conda run -n msgfjava "$@"; }
fi

[[ -f "$JAR" ]] || { echo "ERROR: $JAR missing (fetch_reference_data.sh --jar)"; exit 1; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
mkdir -p "$OUTDIR"

JVM javac -cp "$JAR" -d "$WORK" "$HERE/java/DumpFdrMap.java"
JVM java -cp "$WORK:$JAR" DumpFdrMap "$OUT"

echo "wrote $OUT ($(wc -c < "$OUT") bytes)"
