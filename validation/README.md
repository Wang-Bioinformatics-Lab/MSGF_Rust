# validation/ — MSGF_Rust test & regression corpus

This is **the test folder**. It exists to prove the Rust port reproduces reference **MS-GF+
(Java)** numbers, and to catch regressions. The Java implementation is the oracle (see
`../PLAN.md`, decision D1: exact-reproduction-first).

## Layout

```
validation/
├── fetch_reference_data.sh   # download inputs from upstream MS-GF+ (reproducible; not vendored)
├── data/                     # inputs (gitignored — recreate with fetch_reference_data.sh)
│   ├── spectra/  F13.mgf (high-res QE, 1406 spectra), F13_subset25.mgf (fast fixture),
│   │             tiny.pwiz.mzML, [test.mgf via --full]
│   ├── fasta/    BSA, Tryp_Pig_Bov (+ revCat), [iprg2013_human, ecoli via --full]
│   ├── models/   HCD_QExactive_Tryp.param + other high-res .param scoring models
│   ├── config/   Mods / enzymes / activationMethods / protocols / MSGFPlus_Params
│   ├── example/  test.mzid / test.tsv — an MS-GF+-authored result, used as GOLDEN SCHEMA ref
│   └── MANIFEST.sha256       # checksums of everything fetched
├── reference/                # golden generators
│   ├── build_all_golden.sh   # (re)build every golden (--with-java), then run the regression suite
│   ├── make_chemistry_golden.py   # physics-based masses (no Java) — guarded by published calibrants
│   ├── make_spectra_golden.py     # per-spectrum parse facts + peak hashes (no Java)
│   ├── make_param_inventory.py    # .param size/sha256/identity header (no Java)
│   ├── generate_golden.sh    # JVM-only: run MS-GF+ -> mzid -> tsv -> frozen golden json
│   ├── make_fdr_golden.sh    # JVM-only: TargetDecoyAnalysis q-value maps (needs only the jar)
│   ├── parse_msgf_tsv.py     # MS-GF+ tsv -> normalized golden json (with compare tolerances)
│   └── MSGFPlus.jar          # (gitignored) via fetch_reference_data.sh --jar
├── golden/                   # frozen reference outputs — golden/README.md says what is
│   │                         # committed (chemistry/, param_inventory) vs generated locally
│   ├── iprg2013_F13.golden.json   # authoritative MS-GF+ high-res search, 4,133 PSMs [generated]
│   ├── worked_example.golden.json # 2 MS-GF+-authored PSMs (no Java needed)         [generated]
│   ├── fdr/fdrmap_cases.golden.json    # MS-GF+ target-decoy q-value maps, 14 cases  [generated]
│   ├── chemistry/                      # physics-based, ours                        [COMMITTED]
│   ├── spectra/  models/               # UC-derived fixtures                        [generated]
├── regression/
│   ├── run_regression.py     # re-derives every golden from raw data; runs now (no Java/Rust)
│   └── README.md             # the corpus + compare semantics + how Rust plugs in
└── diff_harness/             # run Java+Rust on the same inputs, report drift (later phase)
```

**Regression suite:** `python3 regression/run_regression.py` → currently **1,762 checks, 0 fail**.
See `regression/README.md` for the full corpus and per-family compare semantics.

## Quick start

```bash
cd validation
./fetch_reference_data.sh            # ~4 MB: high-res spectra, small FASTAs, .param models
./fetch_reference_data.sh --full     # + iPRG human FASTA (~50 MB) needed for the F13 golden
./fetch_reference_data.sh --jar      # + MS-GF+ jar (needs Java 11+ to run)

# generate the frozen golden (JVM-only; commit the json it writes)
cd reference && ./generate_golden.sh          # F13.mgf vs iPRG human FASTA (high-res HCD)
```

## Golden schema

`golden/*.golden.json` records, per PSM keyed by `(spec_file, scan, charge, peptide)`:

| field | mzid CV term | TSV column | Rust must match |
|---|---|---|---|
| `raw_score` | MS:1002049 MS-GF:RawScore | MSGFScore | **exact** (integer) |
| `denovo_score` | MS:1002050 MS-GF:DeNovoScore | DeNovoScore | **exact** (integer) |
| `spec_evalue` | MS:1002052 MS-GF:SpecEValue | SpecEValue | `|log10(rust/java)| ≤ 0.05` |
| `evalue` | MS:1002053 MS-GF:EValue | EValue | `|log10(rust/java)| ≤ 0.05` |

Rationale for tolerance on the E-values: exact IEEE reproduction across JVM/Rust and summation
order is not guaranteed; the integer scores, however, should be bit-exact.

## Provenance & license — READ THIS

Every file `fetch_reference_data.sh` downloads originates from the MS-GF+ repository
(`github.com/MSGFPlus/msgfplus`) and is **Copyright The Regents of the University of California**,
licensed for **educational / research / non-profit use** with attribution; **commercial use
requires a UCSD Technology Transfer agreement**. It is *not* an OSI/permissive license.

Consequences, enforced by `.gitignore`:
- These bytes are **not committed** to this repo — the script re-fetches them on demand.
- The `.param` scoring models are likewise not vendored.
- **Golden JSON derived by running MS-GF+ (or embedding its test data) is no longer committed
  either** — it is regenerated locally by `reference/build_all_golden.sh --with-java` (the flag is
  required; without it only the no-Java families build). Only the families that owe nothing to UC
  (`golden/chemistry/`, `golden/models/param_inventory.golden.json`) stay in git.
  See `golden/README.md` and `../LICENSING.md` for the split and the 2026-07-24 cleanup.
- Golden-backed tests therefore **skip** on a checkout that has not fetched/generated them;
  `cargo test` still passes, it just validates less.

The corpus this project *trains* on is a separate matter: `data/training/` holds MassIVE-KB
peptide libraries (`fetch_reference_data.sh --training`), which are **CC0**, not UC-licensed —
that is what makes the shipped model, and the repository, MIT.
