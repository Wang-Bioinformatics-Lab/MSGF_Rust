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

### 5b — Train from real data 🔜 (`msgf-train`)

Training is a **counting pass** over confident PSMs — no learning-rate, fully reproducible (see
`docs/models.md` §4.1). New crate `msgf-train`:

1. **Ingest** a corpus of *(spectrum, confident peptide, charge)* — **MassIVE-KB (CC0)** first
   (human/HCD/tryptic ⇒ target identity `HCD_HighRes_Tryp`); library reference spectra for a v0,
   provenance-linked raw spectra for a calibrated v1. Pull via the repo's MassIVE/USI tooling.
2. **Partition** each PSM by *(charge, parent-mass segment)*, mirroring MS-GF+'s boundaries.
3. **Count** → fragment offset frequencies; signal vs. noise rank histograms (the §6 rows);
   precursor offsets; mass-error dists. Reuse `msgf-chem` for ion m/z, tolerance, `round_half_up`.
4. **Emit** with `write_param` (Step 4) → a `.param` that drops straight into the pipeline.

Open unknown to resolve first: MS-GF+'s exact partition boundaries, `max_rank`,
`error_scaling_factor`, candidate ion-type list and frequency thresholds, and the **noise
population** definition — mine these from `ScoringParameterGeneratorWithErrors` (same class of
reverse-engineering as the format itself).

### 5c — Validate the model 🔜

Two gates (from `docs/models.md` §4.4):

1. **Trainer mechanics (numeric):** run Java `ScoringParamGen` and `msgf-train` on the **same small
   corpus**; the frequency tables should match within counting/float tolerance. Add
   `validation/reference/generate_training_golden.sh`. Proves the trainer independent of corpus.
2. **Downstream parity (statistical):** train on MassIVE-KB, score a **held-out** benchmark (F13, or
   a MassIVE set absent from training), compare **PSM/peptide IDs at 1 % FDR** and SpecEValue rank
   correlation against MS-GF+'s stock model. Bit-exactness is *not* expected (different corpus); the
   bar is parity. Reuses the Phase-6 search harness.

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
| 5b `msgf-train` from MassIVE-KB | 🔜 next |
| 5c Validation gates | 🔜 |
| (M0) `ScoringModel` trait seam | ⏳ deferred (permissive-native format) |

**Immediate next action:** scaffold `msgf-train` and, before writing counters, mine the training
constants from MS-GF+'s `ScoringParameterGeneratorWithErrors` (partition boundaries, `max_rank`,
`error_scaling_factor`, ion-type candidates, noise-population definition) — capture them as a small
golden so the trainer's mechanics gate (5c-1) has something to check against.
