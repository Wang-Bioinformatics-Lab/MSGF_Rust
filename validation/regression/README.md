# validation/regression — the regression suite

`run_regression.py` re-derives every golden fixture from the raw data + authoritative constants
and asserts they match. It runs **today**, with no Java and no Rust, and is the same oracle the
Rust implementation will be checked against once it exists.

```bash
python3 regression/run_regression.py        # -> "OK: N passed, 0 failed" (exit 0)
```

## What it covers (current corpus)

| Family | Golden | Ground-truth source | What Rust must reproduce |
|---|---|---|---|
| **Chemistry** | `golden/chemistry/*.json` | Atomic monoisotopic masses + fixed residue formulas, **guarded against published calibrants** (MRFA, Bradykinin, Angiotensin II, Glu-Fib) | residue masses (1e-6), peptide neutral & `[M+nH]` (1e-4), b/y ion m/z (1e-4) |
| **Spectrum I/O** | `golden/spectra/*.json` | Direct parse of the MGF/mzML bytes | per-spectrum charge, precursor m/z, peak count, canonical peak-list hash; whole-file rolling hash |
| **Scoring models** | `golden/models/param_inventory.golden.json` | The `.param` bytes | file size, sha256, and parsed identity `(activation, resolution, enzyme, protocol)` |
| **Worked example** | `golden/worked_example.golden.json` | MS-GF+-authored `test.tsv` (2 PSMs, no Java needed) | RawScore/DeNovoScore exact; SpecEValue/EValue within tol |
| **MS-GF+ high-res search** | `golden/iprg2013_F13.golden.json` | **Authoritative** MS-GF+ search of `F13.mgf` (1,406 high-res spectra) vs the iPRG human DB — **4,133 PSMs** | RawScore/DeNovoScore exact; SpecEValue/EValue within `|log10 ratio| ≤ 0.05` |

Every golden json carries its own `compare` block (tolerance + assertion kind), so the Rust
harness reads the semantics from the data rather than hard-coding them.

## Regenerating

```bash
cd validation
./fetch_reference_data.sh --all              # data + iPRG FASTA + MS-GF+ jar
reference/build_all_golden.sh --with-java    # rebuild ALL goldens (incl. the F13 search) + run suite
```

The no-Java fixtures are deterministic. The MS-GF+ search golden is pinned by
`golden/iprg2013_F13.provenance.txt` (jar version, params, DB, JVM). Regenerating it is a
deliberate, reviewed action — commit the resulting `*.golden.json` only when the change is
understood.

## How Rust will plug in (next phase)

Each golden becomes a `cargo test` fixture: the Rust crate reads the same input from `data/`,
produces its output, and compares to the golden using the recorded `compare` semantics. The
integer scores (RawScore/DeNovoScore) are asserted **exact**; the E-values within log tolerance.
`run_regression.py` stays as the language-agnostic gate and the differential-testing entry point.
