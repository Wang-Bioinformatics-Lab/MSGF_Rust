#!/usr/bin/env bash
#
# build_all_golden.sh — (re)generate every golden fixture, then run the regression suite.
#
# No-Java fixtures (chemistry, spectra, param inventory, worked example) always build. Everything
# derived by *running* MS-GF+ — the scoring-model dumps, node scores, target-decoy FDR map, the
# rawscore/ family (scored spectrum, RawScore, SpecEValue, HighRes) and the F13 search golden —
# needs a JVM + jar + input data; pass --with-java to run those through the isolated conda env
# (msgfjava).
#
# --with-java covers every family listed in golden/UC_DERIVED.sha256. If you add a golden, add it
# here too, or its test will skip silently on a fresh checkout.
#
# Usage:
#   ./build_all_golden.sh                 # no-Java fixtures + regression
#   ./build_all_golden.sh --with-java     # + every MS-GF+-derived golden (the full corpus)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "### chemistry golden";       python3 "$HERE/make_chemistry_golden.py"
echo "### spectra golden";         python3 "$HERE/make_spectra_golden.py"
echo "### param inventory golden"; python3 "$HERE/make_param_inventory.py"

# The worked example is a re-parse of MS-GF+'s own shipped example result, so it needs no JVM —
# only the fetched data/example/test.tsv.
echo "### worked-example golden"
if [[ -f "$HERE/../data/example/test.tsv" ]]; then
  python3 "$HERE/parse_msgf_tsv.py" "$HERE/../data/example/test.tsv" \
    -o "$HERE/../golden/worked_example.golden.json"
else
  echo "  skip (missing): data/example/test.tsv (fetch_reference_data.sh)"
fi

if [[ "${1:-}" == "--with-java" ]]; then
  echo "### scoring-model goldens (.param dumps)"
  bash "$HERE/make_model_goldens.sh"
  echo "### node-score goldens (getNodeScore / getMissingIonScore)"
  bash "$HERE/make_nodescore_goldens.sh"
  echo "### target-decoy FDR golden (TargetDecoyAnalysis on synthetic score lists)"
  bash "$HERE/make_fdr_golden.sh"
  # Order matters from here down — this is a dependency chain, not a list:
  #   F13 search golden  ->  select_scored_spectra.py picks its 30 spectra from those PSMs
  #     scored-spectrum golden  ->  freezes the selection the next three all reuse
  #       RawScore / SpecEValue / HighRes goldens
  echo "### MS-GF+ high-res search golden (F13)"
  if command -v conda >/dev/null && conda env list | grep -q msgfjava; then
    TAG=iprg2013_F13 conda run -n msgfjava bash "$HERE/generate_golden.sh"
  else
    TAG=iprg2013_F13 bash "$HERE/generate_golden.sh"
  fi
  echo "### scored-spectrum golden (preprocessing + per-node prefix/suffix scores)"
  bash "$HERE/make_scored_spectrum_golden.sh"
  echo "### RawScore golden (FastScorer node sum + DBScanScorer node+edge)"
  bash "$HERE/make_rawscore_golden.sh"
  echo "### SpecEValue golden (generating function: score distribution + tail)"
  bash "$HERE/make_specprob_golden.sh"
  echo "### HighRes scored-spectrum + RawScore goldens (HCD_HighRes_Tryp model)"
  bash "$HERE/make_highres_goldens.sh"
fi

echo "### regression suite"
python3 "$HERE/../regression/run_regression.py"
