# Performance — MSGF_Rust vs. MS-GF+ (Java)

Measured speed of the **MS-GF+ significance scoring** (the generating-function spectral E-value —
the "MSGF scoring" this project reimplements) in Rust versus the reference Java implementation, on
real high-resolution data. **The Rust output is bit-exact** to MS-GF+ (DeNovoScore + SpecEValue
30/30; distributions to ~2e-8), so these are apples-to-apples timings of identical results, not an
approximation trade-off.

## TL;DR

| Workload (per spectrum: preprocess → scored spectrum → generating function → SpecEValue) | Throughput | Time / 1,406 spectra |
|---|--:|--:|
| **MS-GF+ (Java), 1 core** (JIT-warm) | 72 spectra/s | 19.4 s |
| **MSGF_Rust, 1 core** | **89 spectra/s** | 15.8 s |
| **MSGF_Rust, 32 cores** (rayon) | **1,413 spectra/s** | 1.0 s |

Rust is **~1.23× faster single-threaded** and scales near-linearly across cores (this stage is
embarrassingly parallel). This is an **unoptimized** Rust generating function — no SIMD, no
score-distribution reuse yet — so there is substantial headroom (see *Caveats*).

## Environment

- **CPU:** Intel Xeon E5-2667 v3 @ 3.20 GHz, 32 logical cores (2× 8-core, HT).
- **Rust:** 1.94.0, `--release` / `cargo bench` (criterion), `opt-level=3`, thin LTO.
- **Java:** OpenJDK 17, MS-GF+ v2024.03.26 (`MSGFPlus.jar`), default heap up to `-Xmx8000M`.
- **Model:** `HCD_HighRes_Tryp.param` (the model the F13 search uses — `-inst 1` = HighRes).

## Data

- **`F13.mgf`** — the iPRG-2013 high-resolution Q-Exactive set, **1,406 MS/MS spectra**
  (from the MS-GF+ test resources).
- For the full-search timing: searched against the **iPRG-2013 human FASTA** (~40k proteins),
  trypsin, 10 ppm precursor tolerance, isotope error `-ti 0,1`, target-decoy.

## What is measured

The "MSGF scoring" a database search performs **once per spectrum**: preprocess the raw peaks
(precursor filtering + deconvolution + ranking) → build the scored spectrum (per-node scores) →
build the de-novo graph over the `-ti 0,1` isotope mass group → run the generating-function DP →
obtain the score distribution and the SpecEValue (p-value). Both implementations do exactly this,
over the same 1,406 spectra, with edge scoring on and DB-composition amino-acid probabilities.
Java uses `FlexAminoAcidGraph` + `GeneratingFunctionGroup`; Rust uses `msgf-genfunc`.

## Results

### Full SpecEValue per spectrum

| | Rust 1-core | Java 1-core | Rust 32-core |
|---|--:|--:|--:|
| median spectrum | 8.7 ms | ~13.8 ms | — |
| throughput | 89 spectra/s | 72 spectra/s | 1,413 spectra/s |
| 1,406 spectra | 15.8 s | 19.4 s | 1.0 s (wall) |

Java numbers are from JIT-warm passes (a fresh JVM's first pass was ~64 spectra/s before warm-up).
The generating-function DP dominates the per-spectrum cost in both; preprocessing is ~6 µs and the
scored-spectrum build ~0.8 ms (measured separately), i.e. the graph + DP is ~90% of the time.

### End-to-end MS-GF+ search (Java, for context)

A full MS-GF+ search of `F13.mgf` vs the 40k-protein human DB (32 threads, DB index pre-built):

| Phase | Wall time |
|---|--:|
| Load database index | 14.9 s |
| Read spectra | 15.3 s |
| **Search** (candidate scan + SpecEValue) | 50.1 s |
| q-values + write mzIdentML | 1.3 s |
| **Total** | **66.9 s** |

CPU time: **305.6 CPU-seconds** (User 299.5 + Sys 6.1). Note this whole-search number includes the
**database candidate scan**, which MSGF_Rust does **not** implement yet (that's the future
`msgf-search` engine) — so it is context, not a like-for-like comparison. The generating-function
portion of that search is what the table above isolates.

### CPU-hours

For the SpecEValue stage (the part MSGF_Rust implements), 1-core CPU-time:

| | per 1,406 spectra | projected per 100,000 spectra |
|---|--:|--:|
| MSGF_Rust (SpecEValue) | 0.0044 CPU-hr | **0.31 CPU-hr** |
| MS-GF+ generating function | 0.0054 CPU-hr | 0.38 CPU-hr |
| MS-GF+ full search (incl. DB scan) | 0.085 CPU-hr | — (dominated by DB scan) |

On the 32-core box, MSGF_Rust scores the SpecEValue for 100k spectra in **~71 s of wall-clock**.

## Caveats (read these)

- **Bit-exact output.** These timings compare implementations that produce identical SpecEValues
  (validated 30/30 against MS-GF+). Speed is not bought with accuracy.
- **Rust is unoptimized.** The generating function is a straightforward port of MS-GF+'s DP with a
  dense `ScoreDist` per node. No SIMD on the shift-add convolution, no sliding-window reuse of
  score distributions, no arena allocation. The ~1.2× single-core edge is essentially "Rust vs JVM
  with the same algorithm"; the algorithmic-optimization headroom is untapped.
- **No database search yet.** MSGF_Rust implements the per-spectrum *scoring* (through SpecEValue),
  not candidate generation / the fragment-index DB scan. The Java full-search number includes that
  scan (a large share of its 50 s search phase), so total-pipeline times are not comparable until
  the `msgf-search` engine exists.
- **Nominal mass grid.** This workload runs on the nominal (~1 bin/Da) grid, which is what MS-GF+
  uses for F13 and what our SpecEValue matches. The 274×-finer high-precision grid
  (`INTEGER_MASS_SCALER_HIGH_PRECISION`) is a separate high-res mode not exercised here; it is where
  the "super fast for high-res" optimization work will concentrate.

## Reproduce

```bash
# Rust (needs validation/data/ — run validation/fetch_reference_data.sh)
cd rust
cargo bench -p msgf-genfunc --bench genfunc      # SpecEValue: 1-core + 32-core (rayon)
cargo bench -p msgf-scorer  --bench scoring       # scored-spectrum sub-components

# Java generating-function timing (same workload, single-thread)
cd validation/reference/java
conda run -n msgfjava javac -cp ../MSGFPlus.jar -d /tmp/c TimeGenFunc.java
conda run -n msgfjava java -cp /tmp/c:../MSGFPlus.jar TimeGenFunc \
  ../../data/models/HCD_HighRes_Tryp.param ../../data/spectra/F13.mgf \
  ../../data/config/iprg-2013_Mods.txt ../../data/fasta/iprg2013_human.fasta

# Java full end-to-end search (needs --full FASTA + --jar)
conda run -n msgfjava java -Xmx8000M -jar ../MSGFPlus.jar \
  -s ../../data/spectra/F13.mgf -d ../../data/fasta/iprg2013_human.fasta \
  -mod ../../data/config/iprg-2013_Mods.txt -inst 1 -m 3 -e 1 -t 10ppm -ti 0,1 -tda 1 -thread 32 -o /tmp/f13.mzid
```

_Numbers above are from a single run on the machine described; expect ±10% run-to-run._
