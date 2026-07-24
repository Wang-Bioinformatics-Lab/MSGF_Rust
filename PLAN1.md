# PLAN1 — Owning the fragment-scoring model

Execution plan for the workstream that removes the last licensed dependency from MSGF_Rust's
scoring path: the trained **fragment-scoring model** (the `.param` files). Strategy and the broader
retraining picture live in [`docs/models.md`](docs/models.md) (decision **D5**); this doc is the
concrete "understand it → isolate it → document it → write code to create it → make our own"
checklist, with current status.

The goal is a scoring model we can **produce ourselves**, so MSGF_Rust can ship under MIT/CC0 while
keeping the UC model as a bit-exact test oracle.

---

## Step 1 — Understand the model ✅

The fragment-scoring model is a per-*(activation, instrument, enzyme, protocol)* trained table that
turns observed peaks into per-cleavage-site scores. Decoded into `msgf_scorer::ScoringModel`
(`rust/crates/msgf-scorer/src/lib.rs`). The trained payload:

- **Fragment offset frequencies** (`frag_off`) — which ion types to score per partition + how often
  each is seen.
- **Rank distributions** (`rank_dist`) — the load-bearing scores: per ion type, a signal row and a
  noise row indexed by observed-peak intensity rank. The node score is
  `ln( ionFreq[rank] / (noiseFreq[rank] · min(ionCharge, num_segments)) )`
  (`ScoringModel::score_from_table`, `lib.rs:452`).
- **Precursor offset frequencies** (`precursor_off`) and **mass-error distributions** (`error_dist`,
  the high-res term).

Full field-by-field anatomy: `docs/models.md` §1.1.

## Step 2 — Isolate the clearly-licensed part ✅ (boundary drawn; enforced)

The point of the isolation is a clean line between **UC-encumbered artifacts** and **our
clean-room work**, so a permissive release is defensible.

| | Artifact | Status |
|---|---|---|
| **Encumbered (UC, non-commercial)** | The trained *bytes* in `validation/data/models/*.param`; MS-GF+ Java source; any golden JSON derived by running MS-GF+ | **Quarantined:** `.param`/data are gitignored, fetched on demand, used **only** as a test oracle — never vendored, never on the shipping path |
| **Ours (clean-room)** | The `.param` *format documentation*; the reader (`read_param`) and the new writer (`write_param`); any `ScoringModel` we construct; a model trained on openly-licensed data | Committed, MIT-intended |

**Why this line is defensible:** a file *format* is an interface (uncopyrightable); what's licensed
is the specific trained numbers. Our encoder/decoder are written from the documented format, not
transcribed from MS-GF+'s Java serializer — stated explicitly in the `src/write.rs` module header.

**Enforcement, not just intent:** the `author_a_model_from_scratch` test constructs and scores a
model with **zero fetched bytes** — proving the authoring path is independent of any UC artifact.
The remaining tie is that we don't yet *train* a good model; Steps 4–5 close that.

> Later hardening (from `docs/models.md` M0): put a small `ScoringModel`/`NodeScorer` trait behind
> `node_score()` so `msgf-genfunc` is agnostic to where the numbers came from, and the UC `.param`
> impl degrades cleanly to "oracle only." Not required to emit a model, so deferred — but it's the
> clean home for a future permissive-native format.

## Step 3 — Document the format ✅

Byte-level spec: **[`docs/param-format.md`](docs/param-format.md)** — endianness, string encoding,
the eight stream sections, the ordering constraints (TreeSet partition order; §6 skips empty
partitions, §7 doesn't), and the two read-side transforms (derived ion `name`; zero ion-existence
floored to `0.001`). Validated normative — the writer built from it reproduces all four UC models
byte-for-byte.

## Step 4 — Write code that creates a model ✅ (encoder)

`msgf_scorer::write_param(&ScoringModel) -> Vec<u8>` and `write_param_file(path, &m)` — the exact
inverse of `read_param`, in `rust/crates/msgf-scorer/src/write.rs`. Clean-room, ~150 lines,
mirrors the reader section-for-section.

**Evidence (`tests/roundtrip_write.rs`, both green):**

- `real_models_round_trip` — all four UC high-res models re-encode **byte-for-byte identical**
  (741,431 / 429,625 / 314,405 / 269,788 bytes). The guaranteed invariant is
  `read(write(read(f))) == read(f)`; byte-identity holds too because none of the shipped models has
  a zero ion-existence entry.
- `author_a_model_from_scratch` — a hand-built 1-partition model (b + y ions, 4-column rank table)
  is written (193 bytes), re-read to an identical struct, and produces the expected
  `node_score`/`missing_ion_score`.

This is milestone **M1** from `docs/models.md`. It unblocks emitting *trained* models in the
existing format, which the reader/scorer/generating-function already accept unchanged.

## Step 5 — Make our own model 🔜

### 5a — Author a model programmatically ✅ (mechanism proven)

`author_a_model_from_scratch` already assembles a `ScoringModel` from plain numbers and round-trips
it. That's the mechanism a trainer plugs into: produce the tables, hand them to `write_param`.

### 5b — Train from real data ✅ (`msgf-train`)

Training is a **counting pass** over confident PSMs — no learning-rate, fully reproducible (see
`docs/models.md` §4.1). New crate `msgf-train`:

1. **Ingest** a corpus of *(spectrum, confident peptide, charge)* — **MassIVE-KB (CC0)** first
   (human/HCD/tryptic ⇒ target identity `HCD_HighRes_Tryp`); library reference spectra for a v0,
   provenance-linked raw spectra for a calibrated v1. Pull via the repo's MassIVE/USI tooling.
2. **Partition** each PSM by *(charge, parent-mass segment)*, mirroring MS-GF+'s boundaries.
3. **Count** → fragment offset frequencies; signal vs. noise rank histograms (the §6 rows);
   precursor offsets; mass-error dists. Reuse `msgf-chem` for ion m/z, tolerance, `round_half_up`.
4. **Emit** with `write_param` (Step 4) → a `.param` that drops straight into the pipeline.

**Shipped** as `rust/crates/msgf-train` (`corpus.rs` / `partition.rs` / `ions.rs` / `counts.rs` +
the `msgf-train` binary). 64,474 MassIVE-KB PSMs → a 176-partition `.param` in **3.4 s**; counting
only, so the same corpus reproduces the model byte-for-byte.

The open unknowns were resolved **without** reading `ScoringParameterGeneratorWithErrors`: the
statistics are defined from how the scorer consumes each table, then checked against the numbers in
the shipped models (row sums, absent bins, the nominal-grid offsets). Full write-up, including the
per-section definitions and the normalisation evidence: **[`docs/training.md`](docs/training.md)**.

### 5c — Validate the model ◐ (downstream parity measured; mechanics oracle skipped)

Two gates (from `docs/models.md` §4.4):

1. **Trainer mechanics (numeric):** ~~run Java `ScoringParamGen` on the same corpus~~ — **dropped
   deliberately.** Our statistics are defined independently of theirs (that is the clean-room
   position), so bin-for-bin agreement is not the right expectation. Reproducibility is pinned
   instead by `tests/train_smoke.rs` (synthetic corpus, zero fetched bytes; byte-identical models
   across runs).
2. **Downstream parity (statistical):** ✅ measured with `validation/eval_trained_model.py` — both
   models rescore the same PSM list of identified peptides + 5 mass-identical shuffled decoys, so
   the gate needs no external FDR oracle (F13's own q-values are degenerate). Results:
   **held-out MassIVE-KB — 1475 vs 1482 IDs at a 1 %-decoy threshold (99.5 % of the UC model),
   ρ = 0.977 on log10 SpecEValue; F13 raw — 61 vs 66 IDs (92 %), ρ = 0.937.** Scoring costs +7 %
   time. The residual gap is corpus-domain, not code: consensus library spectra show a 0.88 y-ion
   hit rate vs 0.66 in raw data, so the model expects a more complete ion series than raw spectra
   deliver. See `docs/training.md`.

Throughout, the existing UC-`.param` golden tests stay green as the independent regression oracle.

---

## Status & next actions

| Step | State |
|---|---|
| 1 Understand | ✅ documented (`docs/models.md` §1) |
| 2 Isolate license boundary | ✅ boundary drawn + enforced by `author_a_model_from_scratch` |
| 3 Document format | ✅ `docs/param-format.md` |
| 4 Encoder (`write_param`) | ✅ shipped; 4/4 UC models byte-exact + synthetic model round-trips |
| 5a Author programmatically | ✅ proven |
| 5b `msgf-train` from MassIVE-KB | ✅ shipped; 64k CC0 PSMs → model in 3.4 s, reproducible |
| 5c Validation gates | ◐ parity measured (99.5 % held-out / 92 % F13); mechanics oracle dropped by design |
| (M0) `ScoringModel` trait seam | ⏳ deferred (permissive-native format) |

**Shipped (2026-07-24):** the trained model is now the **default** —
`msgf-scorer/models/MSGFRust_HCD_HighRes_Tryp_v1.param` (258,352 MassIVE-KB PSMs, 236 partitions)
is embedded in the crate and used by every CLI subcommand when `--param` is omitted. The UC-derived
goldens were untracked at the same time, so **nothing MSGF_Rust distributes is UC-licensed** and the
repo carries an MIT `LICENSE`; the full accounting is in `LICENSING.md`. That closes the release
blocker D5/M4 named for the *shipping path* — the remaining M4 item is model quality from raw
spectra, below.

**Immediate next action:** close the corpus gap (the quality half of milestone **M4**). The v0 model is trained on
MassIVE-KB *consensus* spectra, which are more ion-complete than real acquisitions and have their
precursor region stripped; both limits are measured in `docs/training.md`. Training from the raw
acquisitions MassIVE-KB's provenance links to needs an **mzML reader in `msgf-io`** (not present) —
that reader, not the trainer, is now the blocking work item.
