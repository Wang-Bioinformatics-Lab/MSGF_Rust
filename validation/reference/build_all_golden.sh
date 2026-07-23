#!/usr/bin/env bash
#
# build_all_golden.sh — (re)generate every golden fixture, then run the regression suite.
#
# No-Java fixtures (chemistry, spectra, param inventory) always build. The authoritative MS-GF+
# search golden builds only if a JVM + jar + input data are available; pass --with-java to run it
# through the isolated conda env (msgfjava).
#
# Usage:
#   ./build_all_golden.sh                 # no-Java fixtures + regression
#   ./build_all_golden.sh --with-java     # also (re)run the MS-GF+ F13 high-res search golden
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "### chemistry golden";       python3 "$HERE/make_chemistry_golden.py"
echo "### spectra golden";         python3 "$HERE/make_spectra_golden.py"
echo "### param inventory golden"; python3 "$HERE/make_param_inventory.py"

if [[ "${1:-}" == "--with-java" ]]; then
  echo "### scoring-model goldens (.param dumps)"
  bash "$HERE/make_model_goldens.sh"
  echo "### node-score goldens (getNodeScore / getMissingIonScore)"
  bash "$HERE/make_nodescore_goldens.sh"
  echo "### MS-GF+ high-res search golden (F13)"
  if command -v conda >/dev/null && conda env list | grep -q msgfjava; then
    TAG=iprg2013_F13 conda run -n msgfjava bash "$HERE/generate_golden.sh"
  else
    TAG=iprg2013_F13 bash "$HERE/generate_golden.sh"
  fi
fi

echo "### regression suite"
python3 "$HERE/../regression/run_regression.py"
