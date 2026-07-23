#!/usr/bin/env bash
#
# generate_golden.sh — run the reference MS-GF+ (Java) implementation to produce the frozen
# golden outputs the Rust port is validated against. This is the ONLY step that requires a JVM;
# once the golden JSON in validation/golden/ is committed, day-to-day Rust tests need no Java.
#
# Prereqs:
#   - Java 11+ on PATH               (the reference jar targets Java 11+)
#   - validation/reference/MSGFPlus.jar   (run: ../fetch_reference_data.sh --jar)
#   - input spectra + FASTA in validation/data/  (run: ../fetch_reference_data.sh [--full])
#
# Default run: the high-res iPRG-2013 set (F13.mgf) vs the iPRG human FASTA (needs --full data).
# Override with env vars:
#   SPECTRA=../data/spectra/F13.mgf  DB=../data/fasta/iprg2013_human.fasta \
#   MODS=../data/config/iprg-2013_Mods.txt  INST=1 FRAG=3 ENZYME=1 TOL=10ppm \
#   TAG=iprg2013_F13  ./generate_golden.sh
#
# INST/FRAG codes match MS-GF+: INST 1=Orbitrap/high-res, FRAG 3=HCD, ENZYME 1=Trypsin.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${HERE}/../data"
GOLDEN="${HERE}/../golden"
JAR="${HERE}/MSGFPlus.jar"

SPECTRA="${SPECTRA:-${DATA}/spectra/F13.mgf}"
DB="${DB:-${DATA}/fasta/iprg2013_human.fasta}"
MODS="${MODS:-${DATA}/config/iprg-2013_Mods.txt}"
INST="${INST:-1}"          # 1 = high-res (Orbitrap/FTICR/Lumos)
FRAG="${FRAG:-3}"          # 3 = HCD
ENZYME="${ENZYME:-1}"      # 1 = Trypsin
TOL="${TOL:-10ppm}"
TDA="${TDA:-1}"            # 1 = build target-decoy internally
THREADS="${THREADS:-$(nproc)}"
TAG="${TAG:-$(basename "${SPECTRA%.*}")}"

command -v java >/dev/null || { echo "ERROR: java not on PATH (install JDK 11+)"; exit 1; }
[[ -f "$JAR" ]] || { echo "ERROR: $JAR missing (run ../fetch_reference_data.sh --jar)"; exit 1; }
[[ -f "$SPECTRA" ]] || { echo "ERROR: spectra $SPECTRA missing (run ../fetch_reference_data.sh [--full])"; exit 1; }
[[ -f "$DB" ]] || { echo "ERROR: db $DB missing (run ../fetch_reference_data.sh --full)"; exit 1; }

mkdir -p "$GOLDEN"
work="$(mktemp -d)"
mzid="${work}/${TAG}.mzid"
tsv="${work}/${TAG}.tsv"

echo "==> MS-GF+ search: $(basename "$SPECTRA") vs $(basename "$DB")  [inst=$INST frag=$FRAG enz=$ENZYME tol=$TOL]"
java -Xmx3500M -jar "$JAR" \
  -s "$SPECTRA" -d "$DB" ${MODS:+-mod "$MODS"} \
  -inst "$INST" -m "$FRAG" -e "$ENZYME" -t "$TOL" -tda "$TDA" \
  -thread "$THREADS" -o "$mzid"

echo "==> mzid -> tsv (MzIDToTsv, keep all scores)"
java -Xmx2000M -cp "$JAR" edu.ucsd.msjava.ui.MzIDToTsv -i "$mzid" -o "$tsv" -showQValue 1 -showDecoy 1 -unroll 1

echo "==> tsv -> frozen golden json"
python3 "${HERE}/parse_msgf_tsv.py" "$tsv" -o "${GOLDEN}/${TAG}.golden.json"

# keep the raw reference artifacts next to the json for provenance
cp "$mzid" "${GOLDEN}/${TAG}.mzid"
cp "$tsv"  "${GOLDEN}/${TAG}.tsv"
rm -rf "$work"

# record exactly how this golden was produced, for reproducibility
cat > "${GOLDEN}/${TAG}.provenance.txt" <<EOF
tag=${TAG}
jar=$(basename "$JAR")   jar_source=$(cat "${HERE}/MSGFPlus.jar.source" 2>/dev/null || echo unknown)
spectra=$(basename "$SPECTRA")
db=$(basename "$DB")
mods=$(basename "$MODS")
params: -inst ${INST} -m ${FRAG} -e ${ENZYME} -t ${TOL} -tda ${TDA}
java=$(java -version 2>&1 | head -1)
EOF
echo "done -> ${GOLDEN}/${TAG}.golden.json (+ .mzid .tsv .provenance.txt)"
