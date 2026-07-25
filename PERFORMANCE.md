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
| **MSGF_Rust, 1 core** | **312 spectra/s** | 4.50 s |
| **MSGF_Rust, 32 cores** (rayon) | **4,209 spectra/s** | 0.33 s |

MSGF_Rust is **~4.3× faster single-threaded** than JIT-warm Java and scales near-linearly across
cores (this stage is embarrassingly parallel). These are the numbers **after optimization**; the
first working port was ~89 spectra/s (~1.23× Java), so the optimization work below is a further
**~3.5× single-thread** (and ~3.0× at 32 cores) on top of that — with **no loss of fidelity**.

## Platform

All numbers are from a single machine; expect ±10% run-to-run.

- **CPU:** Intel Xeon E5-2667 v3 @ 3.20 GHz (Haswell-EP). **Virtualized** (runs under a hypervisor);
  `lscpu` presents **32 logical CPUs** (1 thread/core in this VM) to the process. The 1-core numbers
  are single-threaded; the 32-core numbers use all 32 via rayon.
- **CPU ISA:** SSE4.2, AVX, **AVX2, FMA**; **no AVX-512**. This matters for the DP: its convolution
  kernel selects an **AVX** (256-bit, 4×f64) path at **runtime** (`is_x86_feature_detected!`), and
  deliberately does **not** use FMA (fused multiply-add would change rounding and break
  bit-exactness). Caches as reported by the VM: L1d 32 KiB/core, L3 16 MiB.
- **Rust:** 1.94.0. Built with `cargo build --release` / `cargo bench` (criterion). `profile.release`
  = `opt-level = 3`, `lto = "thin"`, `codegen-units = 1`. **No `-C target-cpu` flag** — because the
  AVX path is chosen at runtime, the *default* release binary is exactly what is measured (a stock
  `cargo build --release` gets the same speed on any AVX-capable x86-64; other targets fall back to a
  scalar kernel with identical results).
- **Java:** OpenJDK 17, MS-GF+ v2024.03.26 (`MSGFPlus.jar`), heap up to `-Xmx8000M`.
- **Model / grid:** `HCD_HighRes_Tryp.param` (the model the F13 search uses — `-inst 1` = HighRes),
  on the **nominal** (~1 bin/Da) mass grid. Edge scoring on, DB-composition amino-acid
  probabilities, `-ti 0,1` isotope-error group.

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
| throughput | 312 spectra/s | 72 spectra/s | 4,209 spectra/s |
| 1,406 spectra | 4.50 s | 19.4 s | 0.33 s (wall) |

Java numbers are from JIT-warm passes (a fresh JVM's first pass was ~64 spectra/s before warm-up)
and are the same reference figures as before optimization — only the Rust side changed.

### Where the Rust time goes now

Single-thread, over the 1,406-spectrum pass (from `examples/profile`), the SpecEValue pipeline
breaks down as:

| Stage | Share |
|---|--:|
| preprocess + scored-spectrum | ~1% |
| per-spectrum node tables (node mass, prefix/suffix node scores) | ~18% |
| graph build (CSR edges, once per spectrum) | ~11% |
| **generating-function DP (score-distribution convolution)** | **~70%** |

The DP now dominates; graph construction and the (previously dominant) scoring lookups were driven
down by the changes below. Allocation is effectively gone: ~19.8M allocator calls per run in the
first port → **~0** (a few dozen), and allocation traffic 9.0 GB → ~0.75 GB.

## How it got fast (all bit-exact)

Every change was verified against the frozen MS-GF+ golden set (DeNovoScore exact, SpecEValue within
`|log10 ratio| ≤ 0.05`; the checked-in run matches 30/30). Profiling — not intuition — chose each
target; the first assumption (that allocation was the bottleneck) was wrong and the profiler caught
it.

1. **Zero-allocation data layout.** The graph moved from `Vec<Node>` + a per-node `Vec<Edge>` (which
   reallocated on every push) to a flat **CSR** layout built in two counting passes, and the DP runs
   over a single **reusable arena** instead of a heap `ScoreDist` per node. This removed essentially
   all per-spectrum allocation but, tellingly, only ~10% of the time — allocation was never the
   bottleneck.
2. **Removed redundant scoring work** (the real win). Cache each node's main-ion mass once per
   spectrum instead of re-resolving it per incident edge; precompute the per-spectrum ion-existence
   and mass-error score tables (removing ~83M `ln()` calls per run); replace the per-lookup binary
   search in peak matching with an O(1) bucket index; and share the candidate-independent per-spectrum
   node tables **and one edge build** across both isotope-error candidate graphs (they differ only in
   node scores).
3. **Vectorized the DP convolution.** The shift-add kernel runs over contiguous slices through a
   runtime-selected **AVX** kernel (packed multiply + add, no FMA contraction), so the default release
   build is vectorized without a `target-cpu` flag and stays byte-identical to the scalar path.

### CPU-hours

For the SpecEValue stage (the part MSGF_Rust implements), 1-core CPU-time:

| | per 1,406 spectra | projected per 100,000 spectra |
|---|--:|--:|
| MSGF_Rust (SpecEValue) | 0.00125 CPU-hr | **0.089 CPU-hr** |
| MS-GF+ generating function | 0.0054 CPU-hr | 0.39 CPU-hr |

On the 32-core box, MSGF_Rust scores the SpecEValue for 100k spectra in **~24 s of wall-clock**.

### End-to-end database search (`msgf search` vs MS-GF+)

Since the `msgf-search` engine landed this *is* a like-for-like comparison: the same 1,406 F13
spectra against the same concatenated target-decoy human database (160,116 proteins = 80,058
targets + 80,058 decoys), same model (`HCD_HighRes_Tryp.param`), same parameters
(`-inst 1 -m 3 -e 1 -t 10ppm -ti 0,1`, `iprg-2013_Mods.txt`), both timed end-to-end including
index construction.

| | Wall | Threads | Peak RSS |
|---|---|---|---|
| **MS-GF+ (Java)** | 64.0 s | 6 | — |
| **MSGF_Rust** | 5.92 s | 6 | 3.6 GB |
| **MSGF_Rust** | **4.88 s** | 32 | 3.6 GB |

**~10.8× thread-matched**, ~13× at full width. MS-GF+ chose 6 threads itself here — it enforces a
250-spectra-per-thread minimum, so a 1,406-spectrum run caps at 6 no matter what `-thread` says;
the 6-thread Rust row is the fair comparison. The Java breakdown is load DB 13.2 s, read spectra
13.4 s, **search 48.8 s**, q-values 0.02 s, write 1.6 s. Rust builds its index in-process on every
run (30.8M peptides → 48.4M modified candidates), and that build — not the scoring — is most of its
5 s, which is why the 6→32 thread gain is modest.

Both find the same single target PSM at 1% FDR (F13 identifies essentially nothing; see PLAN2 §4).

## Caveats (read these)

- **Bit-exact output.** These timings compare implementations that produce identical SpecEValues
  (validated 30/30 against MS-GF+). Speed is not bought with accuracy — no `f32` probabilities, no
  FMA, no pruning beyond what MS-GF+ itself does.
- **Two different workloads below.** The tables above isolate the *generating function* (the DP
  this project exists to make fast); the end-to-end section adds candidate generation and the DB
  scan via `msgf-search`. Don't quote the 4.3× and the ~10.8× as if they measure the same thing —
  the first is a per-spectrum kernel, the second a whole search dominated by index build.
- **Nominal mass grid.** This workload runs on the nominal (~1 bin/Da) grid, which is what MS-GF+
  uses for F13 and what our SpecEValue matches. The 274×-finer high-precision grid
  (`INTEGER_MASS_SCALER_HIGH_PRECISION`) is a separate high-res mode not exercised here; the DP is
  now ~70% of the time and would grow ~274× on that grid, so it is where future work concentrates.
- **DP is near its single-thread floor for this algorithm.** It is memory/port-bound: an AVX2 build
  (`target-cpu=native`) beats the portable AVX kernel by only ~a few percent, and the working set is
  already cache-resident. Larger gains need the finer grid (above) or an algorithmic change that
  would forfeit bit-exactness (e.g. FFT convolution).

## Reproduce

```bash
# Rust (needs validation/data/ — run validation/fetch_reference_data.sh)
cd rust
cargo bench -p msgf-genfunc --bench genfunc         # SpecEValue: 1-core + 32-core (rayon)
cargo run  -p msgf-genfunc --example profile --release   # per-stage time + allocation breakdown

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

# Rust end-to-end search, same inputs (--threads 6 to match what MS-GF+ picks for this spectrum count)
cd rust && cargo build --release -p msgf-cli
/usr/bin/time -v ./target/release/msgf search \
  -s ../validation/data/spectra/F13.mgf -d ../validation/data/fasta/iprg2013_human.revCat.fasta \
  --mods ../validation/data/config/iprg-2013_Mods.txt \
  -p ../validation/data/models/HCD_HighRes_Tryp.param \
  -t 10ppm --ti 0,1 -e 1 --threads 6 -o /tmp/f13.rust.tsv
```

_Numbers above are from a single run on the machine described in **Platform**; expect ±10%
run-to-run._
