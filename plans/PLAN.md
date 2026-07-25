# MSGF_Rust — Plan

A high-performance Rust reimplementation of the **MS-GF+ significance scoring** algorithm
(spectral probability / spectral E-value via the generating function), targeting
**high-resolution** tandem MS data, with a validation harness that proves numerical fidelity
to the reference Java implementation. A Sage-inspired search engine is built around it later.

---

## 1. Goal & scope

**Primary deliverable:** a Rust function that, given a preprocessed MS/MS spectrum and a
candidate peptide, returns the MS-GF **RawScore**, **DeNovoScore**, and **SpecEValue**
(spectral E-value) — *fast*, especially in high-resolution mode — and reproduces the Java
MS-GF+ numbers within a defined tolerance.

**Explicit non-goals (for the first milestones):** protein inference, full FDR pipeline,
every enzyme/activation/mod combination, retraining scoring models. These come after the core
is validated. The search engine (fragment index, candidate generation, DB search) is Phase 6.

**Why this is worth doing:** MS-GF+'s spectral E-value is the gold-standard rigorous
significance measure in proteomics, but the Java implementation is slow on high-res data
because the generating-function DP runs on a ~274× finer mass grid (see §4). A cache-friendly,
SIMD-and-rayon Rust core can plausibly be 10–50× faster per spectrum and scale linearly across
cores — the thing a modern search engine needs.

---

## 2. Reference material

### 2.1 MS-GF+ (Java) — the algorithm oracle
Repo: `github.com/MSGFPlus/msgfplus` (package root `edu.ucsd.msjava`).
Papers: Kim, Gupta & Pevzner 2008 (*Spectral Probabilities and Generating Functions…*, JPR);
Kim & Pevzner 2014 (*MS-GF+…*, Nat. Commun.).

Java classes we will port / mine, by concern:

| Concern | Key Java classes (package) | Rust target |
|---|---|---|
| Masses, amino acids, mods, enzymes, tolerance | `msutil`: `AminoAcid(Set)`, `Composition`, `Modification`, `Enzyme`, `IonType`, `ActivationMethod`, `Constants` | `msgf-chem` |
| Mass discretization | `msgf`: `NominalMass(Factory)`, `IntMassFactory`, `MassFactory` | `msgf-chem` |
| Trained scoring model (load) | `msscorer`: `NewRankScorer`, `NewScorerFactory`, `FragmentOffsetFrequency`, `Partition` | `msgf-scorer` |
| Per-spectrum node scores | `msscorer`: `NewScoredSpectrum`, `SimpleDBSearchScorer`; `msgf`: `ScoredSpectrum(Sum)` | `msgf-scorer` |
| De novo graph | `msgf`: `DeNovoGraph`, `GenericDeNovoGraph`, `FlexAminoAcidGraph`, `AminoAcidGraph`, `DeNovoNodeFactory` | `msgf-genfunc` |
| **Generating function DP** | `msgf`: `GeneratingFunction`, `GeneratingFunctionGroup`, `GF`, `ScoreDist(Factory)`, `ProfileGF`, `ScoreBound` | `msgf-genfunc` |
| Spectrum I/O | `mzml`, `parser` packages | `msgf-io` (prefer the `mzdata` crate) |
| CLI / orchestration | `ui`: `MSGFPlus`, `ScoringParamGen` | `msgf-cli` |

Trained model files ship as binary resources under `src/main/resources/ionstat/*.param`,
one per (activation × resolution × enzyme × protocol). **High-res targets first:**
`HCD_QExactive_Tryp.param`, `HCD_HighRes_Tryp.param`, `CID_HighRes_Tryp.param`,
`ETD_HighRes_Tryp.param` (plus TMT/iTRAQ variants later).

### 2.2 Sage (Rust) — the engine blueprint
Repo: `github.com/lazear/sage`, MIT. We borrow its *architecture*, not its scoring: fragment
indexing for fast candidate generation, gzipped-mzML reading, all-cores parallelism, LDA
rescoring + target-decoy FDR. Note Sage scores with a hyperscore, **not** the generating
function — the MSGF spectral E-value is exactly what we add.

### 2.3 License constraint (decision-driving — see §3)
MS-GF+ is **Copyright UC Regents**, licensed for *educational/research/non-profit use with
attribution*; **commercial use requires a UCSD Technology Transfer agreement**. It is not an
OSI license. Directly porting code and redistributing the `.param` models inherits this
restriction. A clean-room reimplementation from the papers, with independently derived or
retrained models, is needed for a permissively licensed (e.g. MIT, like Sage) release.

---

## 3. Key decisions (need sign-off)

**D1 — Fidelity strategy (biggest fork).**
- **(A, recommended) Exact-reproduction-first.** Reuse the `.param` models and mirror the Java
  DP so we reproduce MS-GF+ numbers bit-for-bit (RawScore) / within float tolerance
  (SpecEValue). This gives a *rock-solid regression oracle* and immediate scientific
  credibility. Then, in a later phase, add a high-res-native mode that may deviate but is
  validated statistically. Caveat: inherits the UC license until we swap in own models.
- (B) Clean-room from papers with retrained/derived models. Permissively licensable from day
  one, but no exact oracle — validation becomes statistical (rank correlation, ID-count
  parity) instead of numeric, and it's more research than engineering.
- *Recommendation:* **A now, B later.** Build the validated engine against the reference,
  keep the model layer swappable, and re-license once own models are trained. Flag the license
  clearly in the README.

**D2 — First public API.** The "one super-fast function" =
`spec_evalue(scored_spectrum, peptide) -> {RawScore, DeNovoScore, SpecEValue, EValue}`, with a
batched fast path (build the generating function *once per spectrum*, then evaluate many
candidate peptides as O(1) tail-lookups — see §4). Confirm this is the target shape.

**D3 — mzML parsing.** Reuse the `mzdata` crate (handles mzML/mzMLb/gzip, high-res, used across
the Rust MS ecosystem) rather than porting the Java `mzml` package. Recommend reuse.

**D4 — High-res priority order.** Start with **HCD on Q-Exactive/Orbitrap, tryptic** (the most
common high-res case), then CID/ETD high-res, then labeled (TMT/iTRAQ). Confirm.

**D5 — Own-model / retraining path (the license unlock).** The trained `.param` fragment-scoring
model is the *only* UC-encumbered piece on the scoring path; replacing it with one we train from
**MassIVE-KB (CC0)** is what lets MSGF_Rust go MIT. Adopt the clean-room target: add a
`ScoringModel` seam behind `node_score()` now, build an `msgf-train` crate, ship a CC0-trained
`HCD_HighRes_Tryp` first, and keep the UC `.param` path as a permanent test-only oracle. **Full plan
in [`docs/models.md`](../docs/models.md).** *Recommendation: yes — the concrete execution of D1's "A
now, B later."*

---

## 4. How MSGF scoring works (the load-bearing core)

For a spectrum with precursor mass *M* and a candidate peptide, MS-GF+ produces:

1. **RawScore** — sum, over the peptide's theoretical cleavage sites (prefix masses), of a
   per-site score from the trained model (how strongly observed peaks support b/y/c/z ions at
   that site). Deterministic in (spectrum, peptide, model).
2. **DeNovoScore** — the *maximum* RawScore achievable by *any* peptide of mass *M* (best path
   through the de novo graph); a normalization / upper bound.
3. **SpecEValue** — the significance. Using the **generating function**, compute the
   distribution of RawScore over the ensemble of *all* peptides of mass *M* (each amino-acid
   step weighted by its background frequency). SpecEValue is the tail probability
   `P(score ≥ observed RawScore)`. `EValue ≈ SpecEValue × (#candidate peptides)`.

**The generating function DP.** Build a graph whose nodes are integer *scaled prefix masses*
`0 … round(M·scaler)`; edges `m → m + round(residueMass(aa)·scaler)` for each amino acid *a*,
weighted by `freq(a)`. Each node *m* carries a score `s(m)` from the scored spectrum. Dynamic
program in increasing mass order:

```
ScoreDist[m] = Σ over aa:  shift( ScoreDist[m − massBin(aa)], by s(m) ) · freq(aa)
```

`ScoreDist[m]` is the probability distribution of total scores of all peptides ending at mass
*m*. `ScoreDist[round(M·scaler)]` is the full RawScore distribution → SpecEValue is its upper
tail. **Crucially the DP depends only on the spectrum + M, not on any particular candidate** —
so we build it **once per spectrum** and every PSM against that spectrum is a cheap tail
lookup. This is the batched fast path in D2 and the whole reason MSGF is tractable in a search.

**Why high-res is slow — and our opportunity.** The mass grid is set by `Constants.java`:
`INTEGER_MASS_SCALER = 0.999497` (low-res → a ~4,000-Da peptide is ~4,000 nodes) vs
`INTEGER_MASS_SCALER_HIGH_PRECISION = 274.335215` (high-res → ~1.1M nodes, ~0.0036-Da bins).
The DP cost is `#nodes × ~19 aa × score-support`, so high-res is ~274× heavier per spectrum.
That inner loop, on flat arrays, is exactly where Rust wins.

**Performance levers (Phase 5):**
- Flat `Vec`-backed `ScoreDist` with `[minScore, maxScore]` bounds (Java already bounds
  support; per-node support is narrow and grows with peptide length) — no per-node heap objects.
- **Sliding window:** node *m* only needs predecessors within the max residue mass
  (~186 Da → ~51k high-res bins), so keep a bounded ring of live distributions — cache-resident,
  O(window × support) memory, not O(M × support).
- Integer/fixed-point score arithmetic; `f32` probabilities in log space to avoid underflow
  (spectral probs reach 1e-30+).
- SIMD on the shift-add convolution; `rayon` across spectra (embarrassingly parallel).
- Prune unreachable masses; precompute amino-acid mass bins once.

---

## 5. Repository structure

Honors the requested "rust folder" + "test folder"; Rust unit tests live inside each crate,
while cross-language golden/regression data lives in a top-level `validation/` tree.

```
MSGF_Rust/
├── plans/                     # PLAN.md (this file), PLAN1/2/3
├── rust/                      # THE RUST FOLDER — a Cargo workspace
│   ├── Cargo.toml
│   ├── crates/
│   │   ├── msgf-chem/         # masses, amino acids, mods, enzymes, tolerance, mass scaling
│   │   ├── msgf-scorer/       # .param loader + scored-spectrum (per-node scores)
│   │   ├── msgf-genfunc/      # de novo graph, ScoreDist, generating function, SpecEValue  ← hot core
│   │   ├── msgf-io/           # spectra via `mzdata`; precursor/charge/activation detection
│   │   ├── msgf/              # facade + public API (spec_evalue, SpectralProbabilityTable)
│   │   ├── msgf-cli/          # command line for end-to-end differential testing
│   │   └── msgf-search/       # (Phase 6) Sage-inspired engine: fragment index, FDR
│   └── benches/               # criterion benchmarks
├── validation/                # THE TEST FOLDER — cross-language oracle + regression
│   ├── reference/             # msgfplus.jar (fetch script) + generate_golden.sh
│   ├── data/                  # curated spectra (mzML/mgf), tiny FASTA, peptide lists
│   ├── golden/                # frozen Java outputs: RawScore/DeNovoScore/SpecEValue (+ params)
│   ├── regression/            # frozen expected Rust outputs + per-metric tolerances
│   └── diff_harness/          # run Java+Rust on same inputs, compare, report drift
└── docs/
    ├── porting-map.md         # Java class → Rust module + fidelity status (living doc)
    └── algorithm.md           # derivations, param-file format notes
```

---

## 6. Validation / testing / regression strategy

The whole project is anchored on **the Java implementation as oracle** (per D1-A).

- **Tier 0 — Oracle capture.** `validation/reference/generate_golden.sh` runs Java MS-GF+ on a
  small curated high-res dataset → mzIdentML, and we parse the CV terms
  `MS-GF:RawScore` (MS:1002049), `MS-GF:DeNovoScore` (MS:1002050),
  `MS-GF:SpecEValue` (MS:1002052), `MS-GF:EValue` (MS:1002053) into frozen JSON/TSV in
  `validation/golden/`. Where per-PSM values aren't enough, add a tiny Java probe (or MSGF's
  own `MSGF`/`ScoringParamGen` entry points) to dump intermediate node scores and ScoreDist.
- **Tier 1 — Unit tests** (in-crate `#[cfg(test)]`): amino-acid & fragment masses vs known
  values; tolerance matching; mass-bin rounding vs Java `NominalMass`; `ScoreDist` shift/add
  identities; DP invariants (distribution sums to 1, tail monotone non-increasing in score).
- **Tier 2 — Golden component tests:** fixed (spectrum, peptide, `.param`) → assert
  RawScore/DeNovoScore **exact**, SpecEValue within tolerance (e.g. |Δlog10| ≤ 0.05).
- **Tier 3 — Regression tests:** a frozen corpus in `validation/regression/`; CI fails on any
  drift beyond tolerance. Golden regeneration is a deliberate, reviewed action.
- **Tier 4 — Property tests** (`proptest`): random spectra/peptides → SpecEValue ∈ (0,1],
  monotonic in score, RawScore ≤ DeNovoScore ≤ max of the score distribution's support.
- **Tier 5 — Differential/fuzz harness:** `validation/diff_harness/` runs both implementations
  on randomized inputs and reports the worst discrepancies — catches edge cases the fixed
  corpus misses (unusual charges, missing peaks, terminal mods).

CI (GitHub Actions): fmt + clippy + `cargo test` on every push; a slower job runs Tier 2/3
against checked-in golden data (no JVM needed at PR time); a manual/nightly job regenerates
golden data with the JVM and flags diffs.

---

## 7. Milestones (each ends at a validation gate)

**Phase 0 — Scaffolding & oracle.**
Cargo workspace under `rust/`; CI; `validation/` tree; fetch `msgfplus.jar`; curate a tiny
high-res dataset; `generate_golden.sh` producing frozen golden outputs; decode the `.param`
binary format (read `NewRankScorer`/`ScoringParameterGeneratorWithErrors`) and document it in
`docs/algorithm.md`.
*Gate:* golden JSON exists for ≥1 high-res run; `.param` format documented.

**Phase 1 — `msgf-chem`.**
Amino acids, mono masses, mods, enzymes, tolerance, both mass scalers.
*Gate:* peptide/fragment masses match Java for a fixture set.

**Phase 2 — `msgf-scorer`.**
Parse `HCD_QExactive_Tryp.param`; port `NewScoredSpectrum` (spectrum preprocessing/filters +
per-node scores).
*Gate:* **RawScore exact-matches Java** for the golden PSMs.

**Phase 3 — `msgf-genfunc`.**
De novo graph, `ScoreDist`, `GeneratingFunction`(+`Group`), SpecEValue; the batched
build-once/query-many API.
*Gate:* **DeNovoScore exact**, **SpecEValue within tolerance** on the golden set.

**Phase 4 — `msgf-io` + `msgf-cli`.**
Spectra via `mzdata`; precursor/charge/activation detection; a CLI that scores a spectra file ×
peptide list. Wire up the Tier-5 differential harness end-to-end.
*Gate:* CLI reproduces the Java run over the whole curated dataset within tolerance.

**Phase 5 — Performance.**
Criterion benches; apply §4 levers (flat ScoreDist, sliding window, SIMD, rayon, log-space).
*Gate:* documented speedup vs Java per spectrum and linear multi-core scaling; **no fidelity
regression** (Tier 3 still green).

**Phase 6 — Search engine (`msgf-search`, Sage-inspired).**
Fragment index over a FASTA DB; candidate generation; top-N by RawScore → SpecEValue via the
per-spectrum generating function; target-decoy FDR; mzIdentML/TSV output; parallel over spectra.
*Gate:* ID counts at 1% FDR comparable to Java MS-GF+ on a benchmark dataset.

---

## 8. Risks & open questions

- **License (D1).** Must be resolved before any public/commercial release; keep the model layer
  swappable so we can move to retrained models and a permissive license.
- **`.param` binary format** is undocumented — reverse-engineering from the Java serializer is a
  Phase-0 unknown; budget time for it.
- **Exact float reproduction** of SpecEValue may be impossible bit-for-bit (JVM vs Rust FP,
  summation order); tolerance-based comparison is the pragmatic answer — define tolerances up
  front (§6 Tier 2).
- **Semantics of SpecEValue vs spectral probability vs EValue** must be pinned by reading
  `GeneratingFunction.java` + `MSGFPlus.java`, not assumed.
- **High-res memory** for the DP window needs measurement (Phase 5) — bounded ScoreDist support
  should keep it small, but verify.

---

## 9. Immediate next step

On sign-off of §3, execute **Phase 0**: scaffold `rust/` + `validation/`, fetch the reference
jar, curate a minimal high-res spectrum set, and capture the first frozen golden outputs — so
every subsequent phase is measured against a real oracle from day one.

---

## 10. The implicit scoring models (and how we replace them)

MS-GF+ scoring rests on **two trained/implicit models**, not just the algorithm — and one of them
is the sole reason MSGF_Rust currently inherits UC's non-commercial license:

1. **Fragment-scoring model** — the `.param` files (`msgf-scorer::ScoringModel`). Trained by MS-GF+
   from confident PSMs; **UC Regents, non-commercial.** This is the release blocker.
2. **Amino-acid background-frequency model** — the null P(amino acid) the generating function
   integrates over (`msgf-genfunc`: `AA_PROB = 0.05`, or DB-composition). **Already ours.**

Our loader is **read-only** — there is no trainer or `.param` writer today. The plan to make our own
model — what "training" concretely is (a counting pass), **MassIVE-KB (CC0)** as the corpus, the
`msgf-train` crate, the `ScoringModel` swap seam, validation by trainer-mechanics oracle +
ID-count parity, and the milestones to an MIT release — is the dedicated design doc:

**→ [`docs/models.md`](../docs/models.md)** (decision **D5**).
