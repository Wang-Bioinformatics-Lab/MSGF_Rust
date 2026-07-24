#!/usr/bin/env bash
#
# fetch_reference_data.sh — populate validation/data/ with reference inputs from the
# upstream MS-GF+ (Java) repository, plus (optionally) the MS-GF+ release jar.
#
# We do NOT vendor these bytes into git (see validation/.gitignore); this script makes the
# test corpus reproducible on demand. Provenance + license: everything downloaded here is
# Copyright UC Regents under the MS-GF+ non-commercial/academic license (see README.md).
#
# Usage:
#   ./fetch_reference_data.sh            # small high-res reference set (default, ~4 MB)
#   ./fetch_reference_data.sh --full     # also fetch large FASTAs (iPRG human ~50 MB, ecoli ~2 MB)
#   ./fetch_reference_data.sh --jar      # also fetch + unzip the MS-GF+ release jar (~25 MB, needs Java to run)
#   ./fetch_reference_data.sh --all      # small + full + jar
#   ./fetch_reference_data.sh --training [N]  # also fetch N MassIVE-KB library MGFs (~48 MB each,
#                                             # default 5) into data/training/ — the msgf-train corpus
#   ./fetch_reference_data.sh --force    # re-download even if the file already exists
#
set -euo pipefail

REPO="MSGFPlus/msgfplus"
REF="master"
RAW="https://raw.githubusercontent.com/${REPO}/${REF}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${HERE}/data"
REFDIR="${HERE}/reference"

WANT_FULL=0; WANT_JAR=0; FORCE=0; WANT_TRAINING=0; TRAINING_N=5
prev=""
for arg in "$@"; do
  if [[ "$prev" == "--training" && "$arg" =~ ^[0-9]+$ ]]; then TRAINING_N="$arg"; prev=""; continue; fi
  case "$arg" in
    --full) WANT_FULL=1 ;;
    --jar)  WANT_JAR=1 ;;
    --all)  WANT_FULL=1; WANT_JAR=1 ;;
    --training) WANT_TRAINING=1 ;;
    --force) FORCE=1 ;;
    *) echo "unknown arg: $arg" >&2; exit 2 ;;
  esac
  prev="$arg"
done

# dest_relative_to_data :: upstream_repo_path
SMALL=(
  "spectra/F13.mgf::src/test/resources/iprg-2013/F13.mgf"
  "spectra/tiny.pwiz.mzML::src/test/resources/tiny.pwiz.mzML"
  "fasta/BSA.fasta::src/test/resources/BSA.fasta"
  "fasta/Tryp_Pig_Bov.fasta::src/test/resources/Tryp_Pig_Bov.fasta"
  "fasta/Tryp_Pig_Bov.revCat.fasta::docs/examples/Tryp_Pig_Bov.revCat.fasta"
  "models/HCD_QExactive_Tryp.param::src/main/resources/ionstat/HCD_QExactive_Tryp.param"
  "models/HCD_HighRes_Tryp.param::src/main/resources/ionstat/HCD_HighRes_Tryp.param"
  "models/CID_HighRes_Tryp.param::src/main/resources/ionstat/CID_HighRes_Tryp.param"
  "models/ETD_HighRes_Tryp.param::src/main/resources/ionstat/ETD_HighRes_Tryp.param"
  "example/test.mzid::docs/examples/test.mzid"
  "example/test.tsv::docs/examples/test.tsv"
  "example/test_Unrolled.tsv::docs/examples/test_Unrolled.tsv"
  "config/Mods.txt::docs/examples/Mods.txt"
  "config/enzymes.txt::docs/examples/enzymes.txt"
  "config/activationMethods.txt::docs/examples/activationMethods.txt"
  "config/protocols.txt::docs/examples/protocols.txt"
  "config/MSGFPlus_Params.txt::docs/examples/MSGFPlus_Params.txt"
  "config/iprg-2013_Mods.txt::src/test/resources/iprg-2013/Mods.txt"
)

LARGE=(
  "fasta/iprg2013_human.fasta::src/test/resources/iprg-2013/Homo_sapiens_non-redundant.GRCh37.68.pep.all_FPKM-cRAP.fasta"
  "fasta/ecoli.fasta::src/test/resources/ecoli.fasta"
  "spectra/test.mgf::src/test/resources/test.mgf"          # pairs with example/ expected outputs
)

fetch() {
  local dest="$1" src="$2" out="${DATA}/$1"
  mkdir -p "$(dirname "$out")"
  if [[ -s "$out" && "$FORCE" -eq 0 ]]; then
    echo "  skip (exists): $dest"; return 0
  fi
  echo "  get: $dest"
  curl -fsSL "${RAW}/${src}" -o "$out"
}

echo "==> small reference set -> ${DATA}"
for e in "${SMALL[@]}"; do fetch "${e%%::*}" "${e##*::}"; done

# fast unit-test fixture: first 25 spectra of the high-res set (deterministic subset)
SUB="${DATA}/spectra/F13_subset25.mgf"
if [[ -s "${DATA}/spectra/F13.mgf" && ( ! -s "$SUB" || "$FORCE" -eq 1 ) ]]; then
  echo "==> making fixture: spectra/F13_subset25.mgf (first 25 spectra of F13.mgf)"
  python3 - "$DATA/spectra/F13.mgf" "$SUB" <<'PY'
import sys
src, dst, n = sys.argv[1], sys.argv[2], 25
out, count = [], 0
with open(src) as fh:
    for line in fh:
        out.append(line)
        if line.startswith("END IONS"):
            count += 1
            if count >= n:
                break
open(dst, "w").writelines(out)
print(f"  wrote {count} spectra -> {dst}")
PY
fi

if [[ "$WANT_FULL" -eq 1 ]]; then
  echo "==> large FASTAs -> ${DATA}/fasta"
  for e in "${LARGE[@]}"; do fetch "${e%%::*}" "${e##*::}"; done
fi

if [[ "$WANT_JAR" -eq 1 ]]; then
  echo "==> MS-GF+ release jar -> ${REFDIR}"
  mkdir -p "$REFDIR"
  if [[ -s "${REFDIR}/MSGFPlus.jar" && "$FORCE" -eq 0 ]]; then
    echo "  skip (exists): reference/MSGFPlus.jar"
  else
    url="$(curl -fsSL "https://api.github.com/repos/${REPO}/releases/latest" \
      | python3 -c 'import json,sys; r=json.load(sys.stdin); print(next(a["browser_download_url"] for a in r["assets"] if a["name"].endswith(".zip")))')"
    echo "  release zip: $url"
    tmp="$(mktemp -d)"
    curl -fsSL "$url" -o "${tmp}/msgfplus.zip"
    ( cd "$tmp" && python3 -c 'import zipfile,sys; zipfile.ZipFile("msgfplus.zip").extractall("x")' )
    jar="$(find "$tmp" -iname 'MSGFPlus.jar' -o -iname 'MSGFPlus*.jar' | head -1)"
    cp "$jar" "${REFDIR}/MSGFPlus.jar"
    echo "$url" > "${REFDIR}/MSGFPlus.jar.source"
    rm -rf "$tmp"
    echo "  installed: reference/MSGFPlus.jar"
  fi
fi

# ---------------------------------------------------------------------------------------------
# Training corpus — MassIVE-KB peptide spectral libraries (MassIVE MSV000081142).
#
# DIFFERENT PROVENANCE AND LICENSE from everything above: MassIVE-KB is CC0, not UC-licensed.
# That is the whole point — a model trained only on these bytes carries no upstream restriction
# (see docs/models.md D5). Each file is an annotated MGF (peptide in `SEQ=`), ~48 MB / ~14k
# spectra, and is the input `msgf-train` counts a `.param` from.
if [[ "$WANT_TRAINING" -eq 1 ]]; then
  echo "==> MassIVE-KB training corpus (CC0) -> data/training/  [${TRAINING_N} files]"
  mkdir -p "${DATA}/training"
  PROXY="https://massiveproxy.gnps2.org/massiveproxy/MSV000081142/peak/filtered_library_mgf_files"
  # The first N library shards, by name, so a run is reproducible.
  NAMES=$(curl -fsSL --get "https://datasetcache.gnps2.org/datasette/database.json" \
    --data-urlencode "sql=select filepath from filename where dataset='MSV000081142' and filepath like '%filtered_library_mgf_files%' order by filepath limit ${TRAINING_N}" \
    --data-urlencode "_shape=array" \
    | python3 -c 'import json,sys,os; print("\n".join(os.path.basename(r["filepath"]) for r in json.load(sys.stdin)))')
  for name in $NAMES; do
    dest="${DATA}/training/${name}"
    if [[ -s "$dest" && "$FORCE" -eq 0 ]]; then
      echo "  skip (exists): training/${name}"
    else
      echo "  fetch: training/${name}"
      curl -fsSL "${PROXY}/${name}" -o "$dest"
    fi
  done
  echo "  corpus: $(ls -1 "${DATA}/training" | wc -l) files, $(du -sh "${DATA}/training" | cut -f1)"
fi

echo "==> writing MANIFEST.sha256"
( cd "$DATA" && find . -type f ! -name 'MANIFEST.sha256' ! -path './training/*' -print0 | sort -z \
  | xargs -0 sha256sum > MANIFEST.sha256 )
echo "done. $(wc -l < "${DATA}/MANIFEST.sha256") files in data/ (see data/MANIFEST.sha256)"
