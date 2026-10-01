# Provenance of the clean-room scorer and generating function

This records how the spectrum preparation, RawScore and generating-function code in
`rust/crates/msgf-scorer` and `rust/crates/msgf-genfunc` came to be, who saw what, and how the
result was checked. It combines three records: the spec writer's, the implementer's, and the
integration into this repository. All three steps happened on 2026-09-30.

## Why

Five pieces of MSGF_Rust up to commit `0fb0738` described themselves as written while reading
MS-GF+'s Java source. MS-GF+ is © The Regents of the University of California under a
non-commercial licence. The five pieces were:

- `rust/crates/msgf-scorer/src/preprocess.rs`
- `rust/crates/msgf-scorer/src/scored_spectrum.rs`
- the scoring-lookup half of `rust/crates/msgf-scorer/src/lib.rs`
- `rust/crates/msgf-genfunc/src/lib.rs`
- `rust/crates/msgf-genfunc/src/graph.rs`

A two-party clean-room process replaced them. The replacement is in this repository from the
branch `cleanroom-scorer` onward.

## Step 1: the specification (the "dirty" side)

A Claude agent (Opus) working for the user (mwang87) in the DIA_Proteomics_Rust project wrote
`SPEC.md` (in this directory), a test-vector package and a run-only oracle.

| saw | did not see |
|---|---|
| The published papers listed below. MSGF_Rust's documentation. The five MSGF_Rust files and their callers and tests, read **for behaviour**. The oracle's outputs. | MS-GF+ Java source, jars or decompiled code. They were not fetched. All MS-GF+ Java source had been removed from the machine beforehand. |

What `SPEC.md` contains:

- No code from MS-GF+ or from the five files. It is written in mathematics, prose and pseudocode
  of the spec writer's own structure.
- None of their function, type or variable names, none of their comments, and no Java class or
  method name.
- The one pseudocode block (§3.5, isotope-cluster reduction) fixes observable semantics only.
- Numeric constants (atomic masses, the 0.999497 nominal scaler, isotope spacings, the +2/−11
  cleavage terms, the 50…10000 nominal range). These are facts needed for numeric agreement, not
  expression.

A mechanical self-audit compared every delivered text file with the five source files. It looked
for identical lines of 20 or more characters, shared identifiers, and shared 7-word sequences.
Findings in earlier drafts were rewritten. The final result was 0 identical lines and 0 shared
7-word sequences. The only shared tokens were the published output names (RawScore, DeNovoScore,
SpecEValue), "NaN", one CLI output column name, and one file name in the provenance note.

The spec was validated by the spec writer with a throwaway Python implementation written from
`SPEC.md` alone, kept outside the package. It reproduced the oracle bit for bit on all
intermediate vectors (76 spectrum groups, 395 PSMs) and on two samples of search-mode PSMs
(167 and 246 PSMs). Several rules corrected in that loop are listed in `SPEC.md` §12 as "fixed
from observed behaviour".

**Sources.**

1. Kim S, Gupta N, Pevzner PA. *Spectral probabilities and generating functions of tandem mass
   spectra: a strike against decoy databases.* J Proteome Res 2008;7:3354–63.
2. Kim S, Mischerikow N, Bandeira N, et al. *The generating function of CID, ETD, and CID/ETD pairs
   of tandem mass spectra: applications to database search.* Mol Cell Proteomics 2010;9:2840–52.
3. Kim S, Pevzner PA. *MS-GF+ makes progress towards a universal database search tool for
   proteomics.* Nat Commun 2014;5:5277.
4. MSGF_Rust's MIT documentation (`docs/param-format.md`, `docs/training.md`, `docs/models.md`,
   `LICENSING.md`, the bundled model's README, the CLI help text).
5. Black-box runs of the MSGF_Rust `0fb0738` release binary (the oracle), and a small dump tool
   linking the same crates that printed intermediate quantities.

**Data.** 3,000 spectra from PRIDE PXD005573 (HeLa, Q Exactive HF), 15 synthetic edge-case
spectra, and UniProt UP000005640 (human, composition and search; not redistributed). The test
vectors' layout is described in `TEST_VECTORS.md`.

## Step 2: the implementation (the "clean" side)

A separate Claude agent (Opus) working for the user wrote `rust/src/dda/specprob/` in
DIA_Proteomics_Rust. It was branch `dda-specprob` from `main` `0395f88`, commit `b4485bb`,
merged into DIA_Proteomics_Rust `main` as `b930875`. Licence: MIT OR Apache-2.0, same copyright
holder as this repository.

**Inputs used, and only these:**

- the spec package: `SPEC.md`, the package's provenance note, the test vectors, the MIT model and
  its format document, the package's Python `.param` reader, and the oracle's usage note;
- the three papers, as background (no passage was needed beyond what `SPEC.md` states);
- general knowledge (IEEE-754 arithmetic, Rust, the MGF format);
- DIA_Proteomics_Rust's own code, for its conventions only.

**Not consulted:**

- any MSGF_Rust source;
- any MS-GF+ source or jar, or any decompiled code;
- the spec writer's working files and self-audit.

No web search or download was made.

**The oracle binaries** were only run as black boxes, inputs to outputs. They were never inspected
(no `strings`, disassembly or decompilation). Code structure, names and data layout are the
implementer's own. `SPEC.md` prescribes none of them.

Every stage matched the recorded vectors on the first complete build. The points where `SPEC.md`
needed a reading are in `SPEC_ISSUES.md` (in this directory). None of them changed a number on any
test vector.

## Step 3: integration into MSGF_Rust

A Claude agent (Opus) working for the user did the integration on 2026-09-30, on branch
`cleanroom-scorer` from MSGF_Rust `main` `0fb0738`.

1. Right after cloning, the five pieces were removed with `git rm` **without being opened**:
   `msgf-scorer/src/{preprocess,scored_spectrum,lib}.rs` and `msgf-genfunc/src/{lib,graph}.rs`.
   `msgf-scorer/src/lib.rs` was half encumbered. It was removed whole and rebuilt:
   - The `.param` record type and decoder (`msgf-scorer/src/param.rs`) were written fresh from
     `docs/param-format.md`.
   - Field names and types were taken from how the remaining clean code uses them (`write.rs`,
     `msgf-train`, the tests).
   - Scoring-table semantics were taken from the formula in `docs/param-format.md` and the
     `author_a_model_from_scratch` test.

   During the integration, neither the five files nor their history were read.
2. The clean-room modules were copied from DIA_Proteomics_Rust `b930875`:
   - `model.rs`, `spectrum.rs` and `candidate.rs` went into `msgf-scorer`.
   - `null.rs` became `msgf-genfunc/src/lib.rs`.
   - The constants and rounding helpers went into `msgf-scorer/src/lib.rs`.

   Changes:
   - module paths;
   - the `Cleavage` rule moved to `msgf-scorer/src/cleavage.rs`, since both crates need it;
   - spectrum preparation was factored into `preprocess` / `peak_by_mass` for the trainer, with
     the same code path that `prepare` uses.
3. The callers were rewired to the clean API: `msgf-cli` (`rescore`, `model`), `msgf-search`,
   `msgf-train` and the `msgf` facade. The CLI's commands, flags and output formats were not
   changed.
4. Other code from the old crates was dropped, because it was built on the removed internal API:
   - `msgf-genfunc/src/tilt.rs`;
   - the `profile` and `prunelab` examples;
   - the `graph_bitexact` test;
   - five scorer golden tests that read internal node-score arrays.

   The remaining MS-GF+ golden checks were re-expressed against the new API:
   `msgf-scorer/tests/golden_scoring.rs` and `msgf-genfunc/tests/golden_specprob.rs`.

**Validation.** The old release binary (built from `0fb0738`) was compared with the new one on the
same inputs:

- all seven `msgf rescore` test-vector sets;
- `msgf search` on the HeLa subset against UniProt UP000005640;
- the same rescore and search runs with MS-GF+'s own `.param` models on F13;
- the `decoy` and `fdr` subcommands;
- retraining a model with `msgf-train`.

The results are in the integration report and in `rust/crates/msgf-cli/tests/cleanroom_vectors.rs`,
which repeats the test-vector comparison whenever `MSGF_CLEANROOM_VECTORS` is set.

## Step 4: speed work (2026-10-01)

A Claude agent (Opus) working for the user made the clean-room scorer faster on branch
`cleanroom-speed`, from `cleanroom-scorer` `9a7068e`. It changed only clean-room code
(`msgf-genfunc/src/lib.rs`, `msgf-scorer/src/{spectrum,candidate,lib}.rs`).

What it read and did not read:

- It did **not** read MS-GF+ source (none was on the machine, none was fetched), the five removed
  files, any pre-`9a7068e` version of `msgf-scorer/src/{preprocess,scored_spectrum,lib}.rs` or
  `msgf-genfunc/src/{lib,graph}.rs`, or the dropped `tilt.rs`, in git history or anywhere else.
- It read the clean-room code, `docs/cleanroom/SPEC.md`, `LICENSING.md` and the integration
  report, and took *ideas only* from this repository's prose `PERFORMANCE.md` and
  `ALGORITHMIDEAS.md` (marked historical): sharing one edge build across isotope sinks, a
  runtime-selected AVX kernel without FMA, an ion-major node-table sweep, reusable buffers.
  Everything was implemented fresh against the clean-room code; the register-blocked
  convolution with zero-padded rows, the blocked integer passes and the segment-run node
  tables are new.
- For timing on the cluster, the `0fb0738` tree was built from `git archive` with compiler
  output discarded, as in step 3; only the resulting binary was run.

**Validation.** After every change: the integration harness (58 comparisons against the
`0fb0738` binary, plus two rescores with retrained models of other shapes) was byte-identical,
and `cargo test --workspace` with `MSGF_CLEANROOM_VECTORS` set passed. The rescore vectors and
`search` were also checked on a CPU without AVX2, which runs the baseline (non-dispatched) path.
Timings are in `PERFORMANCE.md` ("Current numbers").

## What remains of MS-GF+

**No MS-GF+ code remains in this repository.** What remains:

- MS-GF+'s published method;
- the `.param` file *format*, which is an interface and is documented in `docs/param-format.md`;
- numeric conventions needed to interoperate with MS-GF+ outputs and models;
- names of MS-GF+ concepts in documentation and in validation tooling. The `validation/reference/java/*.java`
  golden dumpers *call* the MS-GF+ jar to produce test oracles and are not shipped in any build.
