# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A Rust reimplementation of **MS-GF+ significance scoring** — the generating-function spectral
E-value (SpecEValue) for high-resolution tandem MS — validated to be **bit-exact** against the
reference Java MS-GF+ (`github.com/MSGFPlus/msgfplus`). The Java implementation is the numeric
**oracle**: every crate is checked against frozen "golden" outputs derived from it. A
Sage-inspired database search engine (`msgf-search`) is a later phase and does not exist yet.

`PLAN.md` is the authoritative design doc (phases, decisions, algorithm derivation). Read it before
substantial work. `PERFORMANCE.md` has the current Rust-vs-Java timings.

## Commands

All Rust work happens in the `rust/` workspace:

```bash
cd rust
cargo test --workspace                        # unit tests + golden validation
cargo test -p msgf-genfunc                    # one crate
cargo test -p msgf-genfunc --test golden_specprob   # one integration test file
cargo test -p msgf-scorer some_test_name      # one test by name
cargo clippy --workspace --all-targets
cargo fmt --all
cargo bench -p msgf-genfunc --bench genfunc   # SpecEValue DP, 1-core + 32-core (rayon)
cargo bench -p msgf-scorer  --bench scoring   # scored-spectrum sub-components
```

Reference data + regression (from `validation/`):

```bash
cd validation
./fetch_reference_data.sh                 # ~4 MB: high-res spectra, small FASTAs, .param models
./fetch_reference_data.sh --full          # + iPRG human FASTA (~50 MB), needed for the F13 golden
./fetch_reference_data.sh --jar           # + MS-GF+ jar (needs Java 11+ to run it)
python3 regression/run_regression.py      # re-derive every golden from raw data; no Java/Rust needed
```

## The data-absence contract (important)

The reference inputs under `validation/data/` are **UC-licensed and NOT committed** (see
`.gitignore` + `validation/README.md`); `fetch_reference_data.sh` recreates them on demand.
Consequently:

- **Golden integration tests skip gracefully when `validation/data/` is absent** — they locate the
  repo root via `env!("CARGO_MANIFEST_DIR")` joined with `../../..`, check for the model/data files,
  and `return` early with an `eprintln!("skip: ...")` if missing. A fresh clone therefore passes
  `cargo test` even with no data. **When adding a golden test, follow this skip pattern** or CI on a
  clean checkout will fail.
- Never `git add` anything under `validation/data/`, and never vendor `.param` models or the
  MS-GF+ jar. Only *derived numeric facts* (the `validation/golden/*.json`) are committed.

## Architecture

The scoring pipeline is a linear dependency chain of four crates (a spectrum + candidate peptide →
SpecEValue). Each crate is validated independently against its own golden family before the next
builds on it.

```
msgf-io  ──►  msgf-scorer  ──►  msgf-genfunc
(MGF read)   (.param model +    (de novo graph + score-distribution DP → SpecEValue)  ← hot core
              scored spectrum)
   msgf-chem  (masses, residues, fragment ions, tolerance, mass-grid scaling) — used by all
```

- **`msgf-chem`** — atomic/residue/peptide masses, b/y ions, tolerance, and the two mass-grid
  scalers in `scaling`: nominal (`0.999497`) and high-precision (`274.335215`). Java `Math.round`
  is reproduced as `round_half_up` (floor(x+0.5)); use it for **all** score/mass rounding to stay
  bit-exact. No dependencies.
- **`msgf-io`** — `Spectrum`/`Peak` types and a streaming `MgfReader`. Validated by byte-for-byte
  peak-list hashes. (mzML via the `mzdata` crate is planned, not present.)
- **`msgf-scorer`** — `read_param()` decodes the binary `.param` scoring model (`ScoringModel`,
  `Partition`, `FragOff`, `RankDist`, …); `preprocess()` does precursor filtering + deconvolution +
  intensity ranking; `ScoredSpectrum` produces per-node prefix/suffix scores and the full
  per-peptide **RawScore**. This is where MS-GF+'s preprocessing must be mirrored exactly.
- **`msgf-genfunc`** — the load-bearing core. `graph::build_reverse_graph` builds the de novo graph
  (nodes = scaled prefix masses, edges = amino acids weighted by background frequency); `compute()`
  runs the score-distribution DP producing a `GenFunc` with a `ScoreDist`; `merge_group()` combines
  the per-isotope-error graphs (the `-ti 0,1` mass group). DeNovoScore is the max-score path;
  SpecEValue is the upper-tail probability of the RawScore distribution.

### Why the generating function is the whole point

The DP `ScoreDist[m] = Σ_aa shift(ScoreDist[m − massBin(aa)], by s(m)) · freq(aa)` **depends only on
the spectrum + precursor mass, not on any candidate peptide** — so it is built **once per spectrum**
and every PSM becomes a cheap tail lookup. The high-res mass grid is ~274× finer than nominal
(~1.1M nodes vs ~4k), making this inner loop the dominant cost and the main optimization target
(flat `Vec` ScoreDist, sliding window, SIMD, rayon across spectra). See `PLAN.md` §4.

## Fidelity is the contract

The reason this project has value is that the output is **bit-exact to MS-GF+**, not an
approximation. Integer scores (RawScore, DeNovoScore) must match **exactly**; SpecEValue/EValue
within `|log10(rust/java)| ≤ 0.05`. When changing scoring, preprocessing, or the DP, the bar is
"golden tests still green," and reproducing Java's exact arithmetic (rounding mode, summation order,
`.param` decode) matters more than idiomatic Rust. If you must diverge from the Java algorithm,
that is a deliberate, reviewed change — flag it, don't silently "improve" it.

Golden generators live in `validation/reference/` (Python for no-Java fixtures; `*.java` dumpers +
`generate_golden.sh` for JVM-derived ones). Regenerating goldens is a deliberate action, never a
side effect of a code change.
