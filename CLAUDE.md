# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A Rust reimplementation of **MS-GF+ significance scoring** — the generating-function spectral
E-value (SpecEValue) for high-resolution tandem MS — validated to be **bit-exact** against the
reference Java MS-GF+ (`github.com/MSGFPlus/msgfplus`). The Java implementation is the numeric
**oracle**: every crate is checked against frozen "golden" outputs derived from it. A database
search engine (`msgf-search`) is built on top of that core and validated against MS-GF+ on F13;
candidates come from a mass-sorted peptide index, not the fragment index Sage uses.

Three workstreams are live: **fidelity/performance** of the scoring path (`plans/PLAN3.md` is the
current perf plan), **model ownership** — removing the last UC-licensed dependency (the trained
`.param` model) so the project can ship MIT — and **search/FDR** hardening (`plans/PLAN2.md` §4).

Doc map — read the relevant one before substantial work:

| Doc | What it is |
|---|---|
| `plans/PLAN.md` | authoritative design doc: phases, decisions (D1–D5), algorithm derivation |
| `plans/PLAN1.md` | the model-ownership execution plan (own the `.param` → train our own), with status |
| `plans/PLAN2.md` | target-decoy + FDR: decoy FASTA, q-values, `msgf-search` wiring; normative MS-GF+ TDA semantics |
| `plans/PLAN3.md` | spectral p-value acceleration (5–10× on the significance stage): exact DP work, gates, cascade |
| `plans/PLAN4.md` | desktop UI: a zero-dependency loopback web server (`msgf ui`) with an embedded front-end |
| `plans/PLAN5.md` | Nextflow scale-out: scatter on spectra, gather FDR; why the database must never be split |
| `plans/PLAN6.md` | timsTOF (Bruker `.d`) DDA support: direct reader, fragment tolerance, a TOF-trained model |
| `ALGORITHMIDEAS.md` + `research-trials/` | index and detailed reports for algorithm/perf experiments — PLAN3's evidence |
| `docs/models.md` | the two trained/implicit models, what taints what, plan to retrain (decision D5) |
| `docs/training.md` | how `msgf-train` counts a `.param` from a corpus + how the MassIVE-KB model compares to UC's |
| `LICENSING.md` | what the repo ships vs. fetches, why it is MIT-clean, and the clean-room boundary |
| `docs/param-format.md` | normative byte-level `.param` spec — the reference for any encoder/trainer |
| `PERFORMANCE.md` | current Rust-vs-Java timings and how the DP got fast |
| `validation/README.md` | golden corpus layout, provenance, license |

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

# The F13 end-to-end search oracle is #[ignore]d (~48M-candidate index over the human DB, ~2 GB):
cargo test -p msgf-search --release -- --ignored --nocapture
```

The `msgf` binary (crate `msgf-cli`) has four subcommands — `search`, `rescore`, `decoy`, `fdr`.
Full flags live in each subcommand's `USAGE` const (`msgf-cli/src/<subcommand>.rs`); `msgf
<cmd> --help` prints it. **Omitting `--param` silently uses the bundled MassIVE-KB model**, which is
a different scoring function from MS-GF+'s — always pass `--param` when comparing to Java:

```bash
cargo run -p msgf-cli --release -- search -s run.mgf -d human.revCat.fasta \
  --fixed-mod C+57.021464 --var-mod M+15.994915 -t 10ppm -o psms.tsv
cargo run -p msgf-cli --release -- decoy -d human.fasta -o human.revCat.fasta
cargo run -p msgf-cli --release -- rescore \
  --spectra spectra.mgf --param HCD_HighRes_Tryp.param --psms psms.tsv [--ti 0,1] [--aa-probs f.tsv]
cargo run -p msgf-cli --release -- fdr -i rescored.tsv -o rescored.q.tsv
```

Training a scoring model of our own (crate `msgf-train`; corpus = annotated MGF with `SEQ=`):

```bash
cargo run -p msgf-train --release -- --corpus lib.mgf --out MSGFRust_HCD_HighRes_Tryp.param \
  --report train.report.json
python3 validation/compare_models.py A.param B.param          # table-for-table diff
python3 validation/eval_trained_model.py library --mgf held_out.mgf --models A.param B.param --decoys 5
```

Reference data + regression (from `validation/`):

```bash
cd validation
./fetch_reference_data.sh                 # ~4 MB: high-res spectra, small FASTAs, .param models
./fetch_reference_data.sh --training 5    # + MassIVE-KB training corpus (CC0, ~240 MB, msgf-train)
./fetch_reference_data.sh --full          # + iPRG human FASTA (~50 MB), needed for the F13 golden
./fetch_reference_data.sh --jar           # + MS-GF+ jar (needs Java 11+ to run it)
python3 regression/run_regression.py      # re-derive every golden from raw data; no Java/Rust needed
```

## The data-absence contract (important)

The reference inputs under `validation/data/` are **UC-licensed and NOT committed** (see
`.gitignore` + `validation/README.md`); `fetch_reference_data.sh` recreates them on demand.
Consequently:

- **The UC-derived goldens are not committed either** (since 2026-07-24) — `validation/golden/`
  keeps only the families that owe nothing to MS-GF+ (`chemistry/`, `param_inventory`); the rest is
  regenerated by `validation/reference/build_all_golden.sh --with-java` (that flag is what builds
  the MS-GF+-derived families; without it they stay missing and their tests skip). See
  `LICENSING.md` and `validation/golden/README.md`. A new golden derived from MS-GF+ output must be
  gitignored too — and must be wired into `build_all_golden.sh`, or its test skips forever on a
  fresh checkout.
- **Golden integration tests skip gracefully when `validation/data/` *or the golden* is absent** — they locate the
  repo root via `env!("CARGO_MANIFEST_DIR")` joined with `../../..`, check for the model/data files,
  and `return` early with an `eprintln!("skip: ...")` if missing. A fresh clone therefore passes
  `cargo test` even with no data. **When adding a golden test, follow this skip pattern** or
  `cargo test` on a clean checkout will fail. (There is no test CI workflow — `.github/workflows/`
  holds only `release.yml`, which builds the `msgf` binaries on a `v*` tag. Run the suite locally.)
- Never `git add` anything under `validation/data/`, and never vendor `.param` models or the
  MS-GF+ jar. Only *derived numeric facts* (the `validation/golden/*.json`) are committed.

## The clean-room boundary (constrains how you write code)

Fidelity work reads MS-GF+'s Java freely — reproducing its arithmetic *is* the job for the scorer
and the DP. The **model-authoring path is different**: `write_param` (`msgf-scorer/src/write.rs`)
and any future trainer are deliberately **clean-room**, written from `docs/param-format.md` and this
repo's `read_param`, *not* transcribed from MS-GF+'s `NewRankScorer.writeParameters`. The defense is
that a file format is an interface; what is licensed is UC's trained numbers. Preserve that line —
if you extend the writer or start the trainer, work from the spec doc, and update
`docs/param-format.md` first if the spec is wrong. The `author_a_model_from_scratch` test enforces
the boundary by building and scoring a model with zero fetched bytes; keep it passing.

## Architecture

The scoring pipeline is a linear dependency chain (a spectrum + candidate peptide → SpecEValue).
Each crate is validated independently against its own golden family before the next builds on it.
The search engine is a second layer that drives that chain over database candidates.

```
scoring core:   msgf-io ──► msgf-scorer ──► msgf-genfunc ─┐
                (MGF read)  (.param model +  (de novo graph + DP  │
                             scored spectrum)  → SpecEValue)      │  ↑ hot core
                                  ▲                               │
                            msgf-train (counts a .param we own)   │
                                                                  ▼
search layer:   msgf-db ──────────────────► msgf-search ──► msgf-fdr
                (FASTA, decoys, digestion)   (index + engine)  (q-values)
                                                                  │
                msgf-chem (masses, ions, tolerance, mass-grid scaling) — used by all
                                                                  ▼
                    msgf (one-dependency facade)   msgf-cli (the `msgf` binary)
```

- **`msgf-chem`** — atomic/residue/peptide masses, b/y ions, tolerance, and the two mass-grid
  scalers in `scaling`: nominal (`0.999497`) and high-precision (`274.335215`). Java `Math.round`
  is reproduced as `round_half_up` (floor(x+0.5)); use it for **all** score/mass rounding to stay
  bit-exact. No dependencies.
- **`msgf-io`** — `Spectrum`/`Peak` types and a streaming `MgfReader`. Validated by byte-for-byte
  peak-list hashes. (mzML via the `mzdata` crate is planned, not present.)
- **`msgf-scorer`** — also **ships the default model**: `bundled::model()` decodes
  `models/MSGFRust_HCD_HighRes_Tryp_v1.param` (trained here from CC0 MassIVE-KB), which every CLI
  subcommand uses when `--param` is omitted; its bytes are SHA-pinned by a test, so retraining is a
  deliberate act (update `bundled.rs` + `models/README.md` together).
  `read_param()` decodes the binary `.param` scoring model (`ScoringModel`,
  `Partition`, `FragOff`, `RankDist`, …) and `write_param()` (`write.rs`) re-encodes it, round-trip
  byte-exact on all four high-res UC models; `preprocess()` does precursor filtering +
  deconvolution + intensity ranking; `ScoredSpectrum` produces per-node prefix/suffix scores and the
  full per-peptide **RawScore**. This is where MS-GF+'s preprocessing must be mirrored exactly.
- **`msgf-genfunc`** — the load-bearing core. `graph::build_reverse_graph` builds the de novo graph
  (nodes = scaled prefix masses, edges = amino acids weighted by background frequency) in a flat
  **CSR** layout; `compute()` / `compute_into()` run the score-distribution DP producing a `GenFunc`
  with a `ScoreDist`; `merge_group()` combines the per-isotope-error graphs (the `-ti 0,1` mass
  group). DeNovoScore is the max-score path; SpecEValue is the upper-tail probability of the
  RawScore distribution. On hot paths use `compute_into` with a reused `DpScratch` (one per thread)
  — the arena makes the whole spectrum allocation-free.
- **`msgf-train`** — the trainer: annotated spectra → a `.param` we own (`counts::train`, the
  `msgf-train` binary). Counting only — same corpus, byte-identical model. Clean-room like
  `write_param`: statistics defined from the scorer's consumption semantics and
  `docs/param-format.md`, never transcribed from MS-GF+'s `ScoringParameterGeneratorWithErrors`;
  `tests/train_smoke.rs` trains + scores from a synthetic corpus with zero fetched bytes. Read
  `docs/training.md` before touching the counting definitions.
- **`msgf-db`** — FASTA into a flat `ProteinDb` (one buffer, proteins as offset+length), the
  database's amino-acid composition (this is what weights the de novo graph's edges), enzymes and
  digestion, and MS-GF+-compatible decoy construction validated **byte-for-byte** against reference
  `.revCat.fasta`. Note `DigestParams` defaults to MS-GF+'s **unlimited** missed cleavages.
- **`msgf-fdr`** — target-decoy q-values mirroring `TargetDecoyAnalysis.java`. All arithmetic is
  **`f32`** because that is what Java emits, so q-values compare by *exact equality*. The lookup is
  Java's `TreeMap.higherEntry` (least key strictly greater). Its rules are subtle — the crate doc
  is normative; read it before touching the sweep.
- **`msgf-search`** — the engine. `PeptideIndex::build` generates mass-sorted modified candidates,
  `SearchEngine::run` scores in parallel (rayon), `assign_q_values` is a deliberately serial
  epilogue (FDR is a property of the whole result set). The generating function is built once per
  `(spectrum, charge)` and shared by every candidate in the precursor window.
- **`msgf`** — a facade re-exporting the workspace under short module names so downstreams take one
  dependency. `default-features = false` drops the search layer and its rayon dependency.
- **`msgf-cli`** — the `msgf` binary; one module per subcommand. `rescore` caches one generating
  function per `(scan, charge)` and turns each PSM into a RawScore + tail lookup;
  `tests/golden_rescore.rs` checks the binary end-to-end against MS-GF+ on F13.

### Why the generating function is the whole point

The DP `ScoreDist[m] = Σ_aa shift(ScoreDist[m − massBin(aa)], by s(m)) · freq(aa)` **depends only on
the spectrum + precursor mass, not on any candidate peptide** — so it is built **once per spectrum**
and every PSM becomes a cheap tail lookup. The high-res mass grid is ~274× finer than nominal
(~1.1M nodes vs ~4k), so this inner loop dominates runtime and is the main optimization target. See
`plans/PLAN.md` §4 for the derivation and `PERFORMANCE.md` for what has already been done (CSR graph +
arena, hoisted per-spectrum score tables, O(1) peak bucket index, shared node tables across the two
isotope-error graphs, AVX convolution kernel) — check there before "optimizing" something twice.

## Fidelity is the contract

The reason this project has value is that the output is **bit-exact to MS-GF+**, not an
approximation. Integer scores (RawScore, DeNovoScore) must match **exactly**; SpecEValue/EValue
within `|log10(rust/java)| ≤ 0.05`.

**That contract is conditional on MS-GF+'s own `.param`.** Every golden is generated with one
(usually `HCD_HighRes_Tryp.param`) and every fidelity test passes one explicitly. The *bundled*
default model (`msgf-scorer/models/`, trained here from MassIVE-KB) is a different scoring function
and is deliberately **not** held to bit-exactness. Don't "fix" that, don't write a golden test that
pins the bundled model against MS-GF+ output, and don't let a golden silently fall back to the
bundled default — assert the model you meant to load.

When judging the bundled model, use ground truth, not agreement with MS-GF+. Agreement measures
*sameness*; on a corpus like F13 it measures almost nothing, because MS-GF+'s own top hits there
are 50.0% decoy — chance. The real gate is `validation/eval_trained_model.py library` on a held-out
MassIVE-KB shard (annotated `SEQ=` peptides + mass-identical shuffled decoys), where the two models
are equivalent: true-peptide-above-decoy 0.9988 both. `msgf-scorer/models/README.md` has the
measured comparison and the reproduce command.

When changing scoring, preprocessing, or the DP, the bar is
"golden tests still green," and reproducing Java's exact arithmetic (rounding mode, summation order,
`.param` decode) matters more than idiomatic Rust. If you must diverge from the Java algorithm,
that is a deliberate, reviewed change — flag it, don't silently "improve" it.

**Two divergences are already deliberate** (documented at the top of `msgf-search/src/search.rs`) —
don't "fix" them without reading that doc: cleavage scoring is applied for C-terminal enzymes only
(disabled, with a warning, for Lys-N/Asp-N/unspecific, where MS-GF+ builds the graph in the other
direction), and `EValue = SpecEValue × database size` rather than MS-GF+'s internal candidate-count
estimate. Q-values come from SpecEValue, so the E-value scaling does not touch FDR.

**F13 is a scoring oracle, not an FDR oracle.** MS-GF+'s own F13 output has `QValue == 1.0` for
4132 of 4133 rows and its top hits are R/K-rich junk (`plans/PLAN2.md` §4 measures this). Per-PSM
score comparisons against it are valid; any gate phrased as "IDs at 1% FDR" is not — there is
nothing to compare. Don't build one on F13; the benchmark question is open.

This binds optimizations too: every perf change so far is bit-exact, and the vectorized DP kernel
must stay byte-identical to the scalar one. That is why `axpy_avx` does a packed multiply then a
packed add and **never** an FMA — contraction would change the result. Anything that reassociates
float summation, contracts operations, or reorders accumulation is a fidelity change, not a free
speedup.

Golden generators live in `validation/reference/` (Python for no-Java fixtures; `*.java` dumpers +
`generate_golden.sh` for JVM-derived ones). Regenerating goldens is a deliberate action, never a
side effect of a code change.

## Conventions

`AGENTS.md` carries the repo's commit/PR conventions: scoped Conventional Commit subjects
(`perf(genfunc): …`, `fix(search): …`), and PRs that name the affected crates, the validation
commands run, and any numeric/performance impact.
