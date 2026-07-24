#!/usr/bin/env bash
#
# make_specprob_golden.sh -- build the GENERATING-FUNCTION validation oracle
# (MS-GF:DeNovoScore + spectral probability / MS-GF:SpecEValue) for the Rust port.
#
# For the 30 (scan, charge, peptide) PSMs already frozen in
# golden/rawscore/f13_scored_spectrum.golden.json, SpecProbDumper.java (compiled
# against MSGFPlus.jar) MIRRORS EXACTLY the construction that
# edu.ucsd.msjava.msdbsearch.DBScanner.computeSpecEValue uses to compute the
# reported MS-GF:SpecEValue / MS-GF:DeNovoScore in the real search:
#
#   scorer = NewScorerFactory.get(HCD,QEXACTIVE,TRYPSIN,STANDARD).doNotUseError();
#   ss     = scorer.getScoredSpectrum(spec);            (spec.setCharge(charge))
#   scoredSpec = new DBScanScorer(ss, precursorNominalMass);   // node+edge
#   graph  = FlexAminoAcidGraph(aaSet, massIndex, TRYPSIN, scoredSpec, false,false)
#            registered in a GeneratingFunctionGroup over the -ti 0,1 / 10ppm range;
#   denovo_score = gf.getMaxScore()-1;  spec_prob = gf.getSpectralProbability(raw_score).
#
# It dumps denovo_score / spec_prob, cross-checks them against the search TSV
# (golden/iprg2013_F13.tsv), and emits the full ScoreDist for the first 3 spectra
# so the Rust DP can be validated bin-by-bin. Every number comes from running MS-GF+.
#
# Output: golden/rawscore/f13_specprob.golden.json
#
# JVM-only; needs reference/MSGFPlus.jar, data/models/HCD_QExactive_Tryp.param,
# data/spectra/F13.mgf, data/config/iprg-2013_Mods.txt,
# golden/rawscore/f13_scored_spectrum.golden.json (source of the 30 PSMs) and
# golden/iprg2013_F13.tsv (the search's DeNovoScore/SpecEValue).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JAR="$HERE/MSGFPlus.jar"
# The iprg2013_F13 search used `-inst 1` = HighRes (Orbitrap), which loads HCD_HighRes_Tryp.param
# (NOT QExactive / -inst 3). SpecProbDumper builds the scorer from the HighRes resource; this path
# is passed only as a label / for provenance.
MODEL="$HERE/../data/models/HCD_HighRes_Tryp.param"
MGF="$HERE/../data/spectra/F13.mgf"
MODS="$HERE/../data/config/iprg-2013_Mods.txt"
SCORED="$HERE/../golden/rawscore/f13_scored_spectrum.golden.json"
TSV="$HERE/../golden/iprg2013_F13.tsv"
# Searched database (provenance: db=iprg2013_human.fasta, -tda 1). Its amino-acid COMPOSITION sets the
# GF edge probabilities (DBScanner.setAminoAcidProbabilities) and is REQUIRED to reproduce SpecEValue.
# The revCat (target+decoy) DB has identical composition; target fasta is used here (matches provenance).
DB="$HERE/../data/fasta/iprg2013_human.fasta"
OUTDIR="$HERE/../golden/rawscore"
OUT="$OUTDIR/f13_specprob.golden.json"

# Use a local javac/java if present, else the msgfjava conda env.
command -v javac >/dev/null 2>&1 && JVM() { "$@"; } || JVM() { conda run -n msgfjava "$@"; }
[[ -f "$JAR" ]]    || { echo "ERROR: $JAR missing (fetch_reference_data.sh --jar)"; exit 1; }
[[ -f "$SCORED" ]] || { echo "ERROR: $SCORED missing (run make_scored_spectrum_golden.sh first)"; exit 1; }
[[ -f "$TSV" ]]    || { echo "ERROR: $TSV missing (fetch_reference_data.sh)"; exit 1; }
[[ -f "$DB" ]]     || { echo "ERROR: $DB missing (fetch_reference_data.sh)"; exit 1; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
mkdir -p "$OUTDIR"

# Selection TSV (scan<TAB>charge<TAB>context_peptide<TAB>golden_raw_score), same
# order and flanking-context peptides as make_rawscore_golden.sh.
python3 - "$SCORED" "$WORK/selection.tsv" <<'PY'
import json, sys
scored, out = sys.argv[1], sys.argv[2]
d = json.load(open(scored))
with open(out, "w") as fh:
    for s in d["spectra"]:
        fh.write(f"{s['scan']}\t{s['charge']}\t{s['golden_peptide']}\t{s['golden_raw_score']}\n")
print(f"wrote {len(d['spectra'])} selection rows")
PY

JVM javac -cp "$JAR" -d "$WORK" "$HERE/java/SpecProbDumper.java"
JVM java -cp "$WORK:$JAR" SpecProbDumper "$MODEL" "$MGF" "$MODS" "$DB" "$WORK/selection.tsv" "$TSV" "$OUT"

echo "wrote $OUT ($(wc -c < "$OUT") bytes)"
