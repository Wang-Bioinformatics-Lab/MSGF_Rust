# MS-GF+'s implicit scoring models — inventory, licensing, and the path to our own

**Status:** design doc. Companion to `PLAN.md` (decision **D1** — fidelity strategy — and the new
decision **D5** below). Read `PLAN.md` §2–4 first for context.

## TL;DR

MS-GF+ scoring rests on **two trained/implicit models**, not just the algorithm:

1. A large **fragment-scoring model** (the `.param` files) — the thing that turns observed peaks
   into per-cleavage-site scores. **Trained** by MS-GF+ from a corpus of confident PSMs.
2. A small **amino-acid background-frequency model** — the null distribution the generating
   function integrates over.

Everything about the *code* is already model-agnostic: the loader is parameterized, the DP takes
per-amino-acid probabilities as inputs. **The only thing forcing MSGF_Rust to inherit MS-GF+'s
non-commercial UC license is that we currently ship numbers derived from UC's `.param` files.**
Replace model #1 with one we train ourselves from permissively-licensed data (**MassIVE-KB is
CC0**) and MSGF_Rust becomes MIT-releasable — with a bit-exact regression oracle kept the whole way.

This doc inventories both models, states exactly what taints what, and lays out a concrete plan to
**train a new model from real data**.

---

## 1. Model inventory — the two implicit models

| # | Model | What it is | Where it lives (code) | Where the bytes come from | License |
|---|---|---|---|---|---|
| 1 | **Fragment-scoring model** (`.param`) | Trained rank-scoring model per *(activation, instrument, enzyme, protocol)* | `msgf-scorer`: `read_param()` → `ScoringModel` (`rust/crates/msgf-scorer/src/lib.rs`) | MS-GF+ repo `src/main/resources/ionstat/*.param`, fetched by `validation/fetch_reference_data.sh` into `validation/data/models/` | **UC Regents, non-commercial/academic** |
| 2 | **AA background-frequency model** | The null P(amino acid) the generating function weights edges by | `msgf-genfunc`: `AA_PROB = 0.05` and the per-edge `aa_prob` (`rust/crates/msgf-genfunc/src/graph.rs`); overridable via `msgf-cli --aa-probs <tsv>` (`load_aa_probs`, `msgf-cli/src/main.rs:460`) | Either uniform `1/20`, or counted from the searched FASTA (`count/total`) | **Ours already** (trivial arithmetic / user's data) |

### 1.1 What model #1 (`.param`) actually contains

Decoded by `read_param()` into `ScoringModel` (`msgf-scorer/src/lib.rs`). The trained payload:

- **`partitions`** — the data is sliced by *(precursor charge, parent-mass segment)*; each slice
  gets its own trained tables. (140 partitions in `HCD_QExactive_Tryp`.)
- **`frag_off`** (`FragOff`) — **which fragment ion types are worth scoring** in each partition and
  how often each is observed (`b`, `y`, charge-2 variants, neutral losses, …), as offset +
  frequency. This is the model deciding "in HCD-tryptic-charge-2, score b and y and b–H₂O but not
  a-ions."
- **`rank_dist`** (`RankDist`) — **the actual node scores.** For each scored ion type, a row of
  `max_rank + 1` frequencies indexed by the observed peak's intensity **rank** (rank 1 = most
  intense; the last bin = "ion absent"), plus a parallel `noise` row. The per-site score is
  `log( ionFreq[rank] / (noiseFreq[rank] · min(ionCharge, numSegments)) )` — see
  `ScoringModel::score_from_table` (`lib.rs:452`). **These two frequency tables are the heart of
  the model**; training is, in essence, the exercise of filling them in.
- **`precursor_off`** (`PrecursorOff`) — precursor m/z offset frequencies (charge-reduced species,
  isotopes).
- **`error_dist`** (`ErrorDist`) — per-partition high-res **mass-error** distributions (signal vs.
  noise, ± `error_scaling_factor` bins) + ion-existence priors, used by the high-res scoring layer.

MS-GF+ produces these with its own trainer, `ScoringParameterGeneratorWithErrors` (Java `ui`:
`ScoringParamGen`), run over a labeled PSM corpus. **We have no trainer today** — our loader is
strictly **read-only** (there is no `write_param`; grep confirms only `read_param`/`read_param_file`).

### 1.2 What model #2 (AA background) contains

Just P(amino acid) for the ~20 residues (plus any variable mods as distinct residues). MS-GF+'s
default is **uniform** (`AminoAcid.probability = 0.05`). Our F13 bit-exact result instead uses
**DB-composition** probabilities (`count/total` over the searched FASTA) — see the `pvalue-status`
memory and `graph.rs`. **This model is already ours** and swapping it is a non-event: it's counting
residues, not licensed data. It affects the SpecEValue *distribution*, not the integer RawScore.

---

## 2. Licensing — what taints what

- **Model #1 (`.param`) is the only encumbered piece.** It is Copyright The Regents of the
  University of California, licensed for educational/research/non-profit use with attribution;
  commercial use needs a UCSD Tech-Transfer agreement. Not OSI. (`validation/README.md`, `PLAN.md`
  §2.3.)
- **We already avoid vendoring it** — it is gitignored, never committed, re-fetched on demand
  (`fetch_reference_data.sh`). The committed golden JSON is *derived numeric facts* used as a test
  oracle, but it too descends from UC software, so it can't be the basis of a permissive release.
- **Consequence for a public release:** as long as the shipping scorer needs a UC `.param` to
  produce its numbers, MSGF_Rust cannot be released under MIT. **This was the single release
  blocker; it is now cleared** — `msgf-scorer` ships a MassIVE-KB-trained model as the default
  (2026-07-24), and the UC `.param` path is validation-only. See `LICENSING.md`.
- **The fix** (already the stated intent in `PLAN.md` D1/§8): *keep the model layer swappable and
  move to a model we train ourselves.* This doc makes that concrete.

---

## 3. Swap-ability — three levels, increasing independence

The consumption interface is tiny and already isolates the model: everything downstream needs from
model #1 is `ScoringModel::node_score()` / `missing_ion_score()` (rank → log-likelihood) plus the
`FragOff` ion set and the error tables. That is the seam to swap on.

- **Level 0 — config swap (works today).** `read_param()` already loads *any* `.param` for any
  *(activation, instrument, enzyme, protocol)*. Choosing `HCD_HighRes_Tryp` vs `HCD_QExactive_Tryp`
  is a runtime selection, no code change. Still UC-licensed.
- **Level 1 — format-compatible retrain.** Produce a **new `.param` in MS-GF+'s exact binary
  format** (write the format `read_param` already decodes) from our own training run. Zero loader
  change; drops straight into the existing pipeline and all golden tests still parse it. Good
  interim step — but the *format* is still MS-GF+'s.
- **Level 2 — clean-room permissive (the goal).** Define **our own** model container, add a
  `msgf-train` crate that writes it, and put a small `ScoringModel` trait behind `node_score()` so
  `msgf-genfunc` doesn't care whether the numbers came from a UC `.param` or our CC0-trained model.
  Ship the trained model with the repo under MIT/CC0. The UC `.param` path stays in the tree as a
  **validation-only oracle** (behind the same trait), never shipped.

**Design rule (adopt now):** introduce the `ScoringModel` trait / `NodeScorer` seam *before*
training, so the retrained model is a drop-in and the Java-derived model degrades cleanly to
"test oracle." This is cheap today and avoids a later refactor.

---

## 4. The plan: train a new model from real data (MassIVE-KB)

This is the new work. Goal: a **CC0/MIT-clean fragment-scoring model**, trained from public data,
that scores within ID-count parity of MS-GF+'s stock model — unlocking a permissive release.

### 4.1 What "training" actually is

MS-GF+'s `ScoringParameterGeneratorWithErrors` is a **counting/histogramming** pass over confident
PSMs — no gradient descent, fully reproducible. For each partition *(charge × parent-mass segment)*:

1. **Fragment offset frequencies** — for each candidate ion type (b, y, a, b²⁺, y²⁺, c, z, common
   neutral losses), over every cleavage site of every PSM in the partition, measure the fraction of
   sites with a matching peak within tolerance. **Ion types above a frequency threshold get
   scored** → populates `frag_off`.
2. **Rank distributions** — for each scored ion type, histogram the **intensity rank** of the
   matched peak (1..`max_rank`, plus an "absent" bin) over true cleavage sites → `ionFreq[]`; do the
   same for random/decoy positions → `noiseFreq[]`. These two rows are exactly what
   `score_from_table` consumes → populates `rank_dist`. **This is the load-bearing step.**
3. **Precursor offset frequencies** — histogram precursor m/z offsets → `precursor_off`.
4. **Mass-error distributions** — per-partition signal vs. noise mass-error histograms + ion
   existence → `error_dist` (the high-res term).

Inputs required: *(spectrum, confident peptide, charge)* triples + the theoretical ion model (from
`msgf-chem`) + a fragment tolerance. Nothing else. We already have the ion model and tolerance; the
missing ingredient is a large, permissively-licensed corpus of confident PSMs — which is exactly
what MassIVE-KB is.

### 4.2 Corpus: MassIVE-KB

[MassIVE-KB](https://massive.ucsd.edu/ProteoSAFe/static/massive-kb-libraries.jsp) is a
community-scale set of peptide spectral libraries distilled from ~31 TB of human HCD proteomics:
**~2.1 M precursors across 19,610 proteins**, FDR-controlled at spectrum/precursor/protein level,
each entry carrying **open provenance** back to the original spectra
([Wang et al., *Cell Systems* 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6279426/)).

- **License: CC0 1.0** ([confirmed on the MassIVE-KB docs / dataset
  MSV000081142](https://www.omicsdi.org/dataset/massive/MSV000081142)). A model trained purely on
  CC0 data carries no upstream restriction → **MSGF_Rust can then go MIT.** This is the crux.
- **Fit:** MassIVE-KB is human **HCD**-dominated, which lines up with the project's D4 priority
  (**HCD / high-res / tryptic first**). The first retrained model is therefore
  `HCD_HighRes_Tryp` — the same identity we already validate bit-exactly.
- **This is the user's own infrastructure.** Pulling it is well within reach: the repo already has
  MassIVE/USI download skills (`massive`, `massive-usi-download`, `dataset-filename-cache`) and the
  `massiveproxy`/GNPS2 routes.

**Two ways to get training PSMs from it — do (a) first, then (b):**

- **(a) v0 — train on the library reference spectra directly.** Each MassIVE-KB entry is a
  peptide + a reference (consensus) spectrum + charge. Already FDR-controlled and immediately
  usable. Fast path to a working model. Caveat: consensus spectra have cleaner noise statistics
  than raw scans, so `noiseFreq[]` will be optimistic — fine for a v0.
- **(b) v1 — follow provenance USIs to the original raw spectra.** Each entry links back to the
  raw acquisitions; training on those matches the real scoring-time peak/noise distribution and
  yields calibrated `noiseFreq[]`/`error_dist`. Heavier (bulk raw fetch) but the user's
  USI/provenance tooling makes it tractable.

### 4.3 New crate: `msgf-train` — **built**; see [`training.md`](training.md) for what it counts

A dedicated crate (sibling to `msgf-scorer`), read-corpus → write-model:

```
msgf-train/
  corpus.rs   # ingest MassIVE-KB (library .mgf/.json or provenance USIs) -> (spectrum, peptide, charge)
  partition.rs# assign each PSM to a (charge, parent-mass segment) partition (mirror MS-GF+ boundaries)
  counts.rs   # frag-offset freqs, rank histograms (signal + noise/decoy), precursor offsets, error dist
  emit.rs     # write model: Level-1 (.param binary, round-trips read_param) OR Level-2 (our format)
  bin/train.rs# CLI: msgf-train --corpus <...> --activation HCD --instrument HighRes --enzyme Tryp -o model
```

The counting logic re-uses `msgf-chem` (ion m/z, tolerance, `round_half_up`) and mirrors the exact
bin/threshold definitions in MS-GF+'s `ScoringParameterGeneratorWithErrors` (partition boundaries,
`max_rank`, `error_scaling_factor`, ion-type candidate list). Those constants must be read out of
the Java source — budget time for it, same as the `.param` format reverse-engineering in Phase 0.

### 4.4 Validation — how we trust a model we trained

Two independent gates (this is D1-B's "statistical validation," made concrete):

1. **Trainer-mechanics oracle (numeric).** Run Java `ScoringParamGen` and `msgf-train` on the
   **same small corpus** with the same settings; compare the resulting frequency tables. They
   should match within counting/float tolerance. This proves the *trainer* is faithful independent
   of corpus size — the analogue of our existing golden capture, but for training. (Add under
   `validation/reference/` as a `generate_training_golden.sh`.)
2. **Downstream ID-count parity (statistical).** Train on MassIVE-KB → score a **held-out**
   benchmark (iPRG2013 F13, or a MassIVE dataset absent from training) → compare **PSM/peptide IDs
   at 1 % FDR** against MS-GF+ with its stock `.param`. Target: comparable ID counts and high rank
   correlation of SpecEValue. This is the real-world "is the model good" test and reuses the
   Phase-6 search harness.

Note we **cannot** expect bit-exactness against MS-GF+'s shipped `.param` — we don't have their
exact training corpus or settings — so the shipped-model bar is *parity*, while the *trainer* bar
(gate 1) is numeric. The existing UC-`.param` golden tests stay green throughout as the
independent regression oracle.

### 4.5 The AA background model (model #2) — essentially free

Retraining model #2 is counting residues in a reference proteome (e.g. UniProt human) or, at search
time, in the searched FASTA — the DB-composition path already exists (`graph.rs`, F13 result). No
license issue, no new corpus. Provide a small default table (permissive) for the no-DB case and
keep the runtime DB-composition path. Done essentially for free alongside #1.

---

## 5. Milestones

Each ends at a gate; ordered so value lands early and the release blocker clears last.

- **M0 — Seam.** Introduce the `ScoringModel`/`NodeScorer` trait behind `node_score()` so genfunc
  is model-source-agnostic; UC `.param` becomes one impl. *Gate:* all golden tests still green
  through the trait. *(Small, do first.)*
- **M1 — Format writer (Level 1).** Implement `.param` **write** that round-trips `read_param`
  byte-for-byte on the four high-res models. `read_param` is the exact spec: the writer must
  reproduce the partition **TreeSet re-sort** (`lib.rs:261`) and the **empty-partition skip** in the
  rank-distribution section (`lib.rs:318`), plus big-endian / UTF-16BE strings / the `0x7FFFFFFF`
  terminator, or the round-trip won't re-parse. *Gate:* `read(write(m)) == m` and the re-emitted
  file reproduces the model golden. *(Unlocks emitting trained models in the existing format.)*
- **M2 — Trainer core.** ✅ `msgf-train` counting pipeline (`rust/crates/msgf-train`). The Java
  `ScoringParamGen` mechanics oracle was dropped by design — our statistics are defined from the
  scorer's consumption semantics, not transcribed, so bin-for-bin agreement is not the expectation.
  *Gate met instead by:* reproducible counting + a synthetic-corpus round-trip test with zero
  fetched bytes.
- **M3 — MassIVE-KB v0 model.** ✅ Trained `HCD_HighRes_Tryp` from 64,474 CC0 MassIVE-KB library
  PSMs in 3.4 s; it rediscovers the same ten ion types as the UC model. *Gate:* 99.5 % of the UC
  model's IDs on held-out MassIVE-KB (ground-truth peptides), 92 % on F13 raw spectra, ρ = 0.94–0.98
  on log10 SpecEValue, +7 % scoring time. Write-up: [`training.md`](training.md).
- **M4a — Ship the model (done, 2026-07-24).** ✅ The MassIVE-KB model is embedded in `msgf-scorer`
  and is the CLI default; UC-derived goldens are no longer committed; repo has an MIT `LICENSE` and
  a `LICENSING.md` accounting. **Gate met: the shipping path contains no UC-derived bytes.** What
  remains of M4 is quality, not licensing:
- **M4 — Provenance-raw v1 + permissive format (Level 2).** Train from raw provenance spectra
  (4.2b); define + ship our own model container under MIT/CC0; UC `.param` path demoted to
  test-only. *Gate:* v1 ≥ v0 on the benchmark; **repo builds and ships with no UC-derived bytes on
  the shipping path** → MIT release unblocked.
- **M5 — Coverage.** Additional identities as corpora allow (CID/ETD high-res, TMT/iTRAQ,
  non-tryptic). *Gate:* parity per identity.

---

## 6. Decision to record in PLAN.md

**D5 — Own-model / retraining path.** Adopt the **Level-2 clean-room** target with **MassIVE-KB
(CC0)** as the first corpus and **`HCD_HighRes_Tryp`** as the first retrained identity. Introduce
the `ScoringModel` seam now (M0); keep the UC `.param` path as a permanent validation oracle, never
on the shipping path. Re-license MSGF_Rust to MIT once M4's gate is met. *(Recommendation: yes —
this is the concrete execution of D1's "A now, B later.")*

## 7. Risks / open questions

- **MS-GF+ training constants are undocumented** — partition boundaries, `max_rank`,
  `error_scaling_factor`, the candidate ion-type list and frequency thresholds must be mined from
  `ScoringParameterGeneratorWithErrors`. Same reverse-engineering risk as the `.param` format.
- **Consensus-vs-raw noise mismatch** (4.2a vs 4.2b) — v0 `noiseFreq[]` will be optimistic; the
  parity gate must be judged on the raw-trained v1, not just v0.
- **Coverage skew** — MassIVE-KB is human/HCD/tryptic-heavy; CID/ETD/labeled/non-tryptic need other
  CC0 corpora (or are deferred).
- **Decoy/noise model for `noiseFreq[]`** — need to pin exactly how MS-GF+ defines the noise
  population (random positions vs. decoy peptides) to match gate 1; read the Java, don't assume.
- **Parity margin must be agreed up front** (as we did for the |Δlog10| ≤ 0.05 SpecEValue
  tolerance), so "good enough to ship" is objective.
