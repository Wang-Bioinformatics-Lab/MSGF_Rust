# Training a fragment-scoring model (`msgf-train`)

How MSGF_Rust **produces** a `.param` scoring model of its own, what each trained number means,
and how the first model — counted from MassIVE-KB (CC0) — compares to the UC-licensed
`HCD_HighRes_Tryp.param` it is meant to replace.

Companions: `docs/models.md` (why we need our own model — decision **D5**), `plans/PLAN1.md` (the
execution plan; this doc closes step **5b** and reports step **5c**), `docs/param-format.md` (the
container).

## TL;DR

* `msgf-train` counts a complete, format-valid `.param` from annotated spectra. 258,352 MassIVE-KB
  PSMs → a 236-partition model in **12 s** (32 cores); counting only, so the same corpus reproduces
  the bytes exactly.
* **The result is what MSGF_Rust ships and uses by default** —
  `msgf-scorer/models/MSGFRust_HCD_HighRes_Tryp_v1.param`, embedded in the binary and used whenever
  `--param` is omitted. That is what makes the project MIT (`LICENSING.md`).
* It independently rediscovers **the same ion types** MS-GF+'s model scores, in the same order of
  prevalence (all ten at the default threshold; the shipped model's stricter 0.25 keeps the top
  nine), with table scales that line up (≈4.5 scored sites per spectrum per ion).
* On held-out MassIVE-KB spectra with ground-truth peptides it is **at parity**: 1480 vs 1485 IDs at
  a 1 %-decoy threshold, median target SpecEValue 10⁻¹⁴·⁴⁷ vs 10⁻¹⁴·⁶⁵.
* On F13 raw spectra (out-of-domain, and peptides originally *chosen* by the UC model) it is just
  behind: 63 vs 66 IDs (95 %), target-beats-decoy 0.746 vs 0.774, Spearman ρ = 0.94 on
  log₁₀ SpecEValue — though on the confident fifth of those PSMs it identifies *more* (24 vs 20).
* Scoring throughput is **the same** (1690 vs 1686 PSM/s; 3.69 ion types per partition vs 3.92).
* The one thing a library corpus **cannot** teach is precursor filtering — MassIVE-KB consensus
  spectra have the precursor region deleted (measured: 0.000 of charge-2 spectra retain it).

## Running it

```bash
cd validation && ./fetch_reference_data.sh --training 20    # ~950 MB of CC0 corpus
cd ../rust && cargo build --release -p msgf-train

# exactly how the shipped model was produced (see msgf-scorer/models/README.md)
./target/release/msgf-train \
  $(for f in ../validation/data/training/*.mgf; do echo -n "--corpus $f "; done) \
  --ion-threshold 0.25 \
  --out crates/msgf-scorer/models/MSGFRust_HCD_HighRes_Tryp_v1.param \
  --report ../validation/data/trained/train.report.json
```

The corpus is any MGF whose spectra carry `SEQ=` (MassIVE-KB peptide-library MGFs are exactly
this). `msgf-train --help` lists the knobs; the defaults are the HCD/HighRes/Tryptic identity.

Comparing and evaluating models:

```bash
python3 validation/compare_models.py A.param B.param        # table-for-table
python3 validation/eval_trained_model.py library --mgf held_out.mgf --models A.param B.param --decoys 5
python3 validation/eval_trained_model.py f13 --tsv validation/golden/iprg2013_F13.tsv \
        --mgf validation/data/spectra/F13.mgf --models A.param B.param --decoys 5 \
        --aa-probs iprg.tsv --ox-m
```

## What the trainer counts

Training is **histogramming, not optimisation** — no optimiser, no RNG, no learning rate. Same
corpus + same config ⇒ byte-identical model (verified: reordering and renaming the corpus files
reproduces the model byte-for-byte).

Each section of the format is filled from a definition derived from how the *scorer* consumes it
(`ScoringModel::score_from_table` computes `ln(ion[rank] / (noise[rank]·min(charge, segments)))`,
so the trainer's job is to produce that ratio's numerator and denominator):

| Section | Definition |
|---|---|
| §3 partitions | equal-count parent-mass quantiles **of our own corpus**, per charge (target ≈400 PSMs/partition, capped at 30), × `num_segments` m/z segments |
| §4 precursor offsets | fraction of spectra with a peak at each nominal-grid offset from the (charge-reduced) precursor m/z; kept if ≥0.15 **and** ≥2× the median of the scanned window |
| §5 fragment offsets | matched sites / all sites, per candidate ion type per partition; kept above a frequency threshold (default 0.15, shipped model 0.25), capped at 6 per partition; a populated partition always keeps its best ion |
| §6 rank distributions | `row[r]` = sites whose matched peak had intensity rank `r`, ÷ **spectra** in the partition; last bin = absent sites |
| §6 noise row | the same count at the node positions of a **decoy peptide**, averaged over the partition's singly-charged scored ion types |
| §7 signal error | main-ion edge mass error (`cur − prev − theoretical residue`) over true cleavage pairs, normalised |
| §7 noise error | the same over decoy-peptide edges |
| §7 ion existence | the four (cur present, prev present) combinations over true edges |

Three choices are worth calling out because they are where a trainer can go wrong:

**Row normalisation is per spectrum, not per site.** A rank row therefore sums to the average
number of scored sites per spectrum in that partition, and the ion/noise ratio is like-for-like.
This is corroborated by the UC model: its y row for charge 2 / mass 1200 / segment 0 sums to 4.66
with an "absent" bin of 1.60 — an absent *probability* above 1 is impossible, an absent *count per
spectrum* of 1.6 is exactly what a ~4.7-site partition should show.

**The noise population is decoy peptides, not random m/z.** The decoy is a deterministic shuffle of
the identified peptide with the C-terminal residue fixed, so it has the *same* parent mass and
partition; positions colliding with a true node of either series are skipped (a decoy prefix mass
can land on a real y-ion peak). This makes the denominator "what this spectrum offers a wrong
peptide at the same mass" — the question the score is actually asked at search time.

**Sparse high ranks are pooled.** Few spectra have a 120th peak, so a ratio of two such bins is
noise; neighbouring ranks are averaged with a window of ±10 % of the rank, identically for ion and
noise rows, and every bin is add-λ floored so `ln(ion/noise)` stays finite and an unobserved rank
scores ≈0 rather than ±∞.

### Candidate ion types

Textbook backbone fragments — b and y series, each with `−H₂O`, `−NH₃`, `−CO`, `+¹³C`, `+2¹³C`,
over fragment charges 1–2 (24 candidates). Which ones enter a model is decided by counting. Where
two candidates round to the same `.param` ion name (charge-2 losses do), the more frequent wins,
because the name is the model's lookup key.

## Clean-room boundary

`msgf-train` is written from the format spec and from the scorer's consumption semantics — **not**
transcribed from MS-GF+'s `ScoringParameterGeneratorWithErrors`. The statistical definitions above
were derived here and then checked *against the numbers* in the shipped models (row sums, absent
bins, offset grids), which is reverse-engineering of a data format, the same posture
`docs/param-format.md` already takes. `tests/train_smoke.rs` trains, writes, re-reads and scores a
model from a synthetic corpus with **zero fetched bytes**, the trainer-side counterpart of
`author_a_model_from_scratch`.

The corpus itself is **CC0** (MassIVE-KB, MassIVE `MSV000081142`), a different provenance from
everything else under `validation/data/` — which is the entire point of D5.

## How close is it to the UC model?

`HCD_HighRes_Tryp.param` (UC, 17,906 training spectra) vs the shipped model (258,352 MassIVE-KB
PSMs, `--ion-threshold 0.25`):

| | UC | ours |
|---|---|---|
| bytes / partitions | 429,625 / 92 | 1,067,562 / 236 |
| partitions per charge | 2:26, 3:17, 4:3 | 2:30, 3:30, 4:30, 5:23, 6:5 |
| ion types found | y, y+i, b, y−NH₃, y+2i, a, y−H₂O, b−H₂O, b+i, b−NH₃ | **the same set**, same prevalence order (b−NH₃ falls below threshold) |
| mean ion types / partition | 3.92 | 3.69 |
| precursor offsets | 14 | 1 counted (+ chemistry fallback) |

At the 0.15 default threshold the model keeps 5.3 ion types per partition and scores marginally
better on library-like data but worse on raw F13 (56 vs 64 IDs) while costing 7 % more time, which
is why the shipped model uses 0.25.

Scale of the load-bearing table (charge 2, mass 1200, segment 0):

| row | UC sum | UC hit rate | ours sum | ours hit rate |
|---|---|---|---|---|
| y (`S_1_19`) | 4.66 | 0.657 | 4.50 | **0.878** |
| b (`P_1_1`) | 5.11 | 0.473 | 4.91 | **0.600** |
| noise | 4.79 | 0.062 | 3.62 | **0.118** |

The row *scales* agree closely — independent evidence that the normalisation above is the one the
format intends. The **hit rates do not**, and that is the corpus talking: MassIVE-KB reference
spectra are consensus spectra, so a real ion series is far more complete there (0.88 of y sites
matched vs 0.66 in raw data) and a decoy position finds a peak twice as often. `plans/PLAN1.md` §4.2
predicted exactly this ("consensus spectra have cleaner noise statistics… `noiseFreq[]` will be
optimistic — fine for a v0"); the number to quote for it is **0.88 vs 0.66**.

The consequence is visible in the scores: our model charges −1.81 for a missing y ion where the UC
model charges −1.03, so on raw spectra — where y ions go missing more often than our corpus taught
us to expect — every correct PSM pays a larger penalty. Node scores for a matched y ion agree
within ≈1 log unit across ranks.

## Does it work? (step 5c evidence)

Both models rescore an identical PSM list: the identified peptide plus five mass-identical shuffled
decoys per spectrum. Metrics are decoy-referenced, so no external FDR oracle is needed — which
matters, because F13's own q-values are degenerate (see the `f13-degenerate-fdr-oracle` note).

**Held-out MassIVE-KB — 1,500 spectra, ground-truth peptides, trained on the other 19 shards:**

| | UC | ours |
|---|---|---|
| median target log₁₀ SpecEValue | −14.65 | −14.47 |
| target beats decoy | 0.999 | 0.999 |
| median log₁₀ gap target↔decoy | 10.04 | 9.92 |
| IDs at 1 %-decoy threshold | 1485 | **1480 (99.7 %)** |
| Spearman ρ (log₁₀ SpecEValue) | — | 0.977 |

**F13 raw spectra — 1,252 spectra, out-of-domain:**

| | UC | ours |
|---|---|---|
| median target log₁₀ SpecEValue | −5.50 | −5.35 |
| target beats decoy | 0.774 | 0.746 |
| IDs at 1 %-decoy threshold | 66 | **63 (95 %)** |
| Spearman ρ (log₁₀ SpecEValue / RawScore) | — | 0.937 / 0.925 |
| top-20 %-confidence stratum, target beats decoy | 0.936 | 0.928 |
| top-20 %-confidence stratum, IDs at 1 %-decoy | 20 | **24** |

Two caveats on the F13 column, both favouring the UC model: the peptides being rescored were
*selected* by MS-GF+ using that model, and F13's identifications are largely junk regardless
(4132/4133 rows have QValue 1). The stratified rows show the deficit is not concentrated in the
confident subset, i.e. it is a genuine but small model difference, not an artefact of the junk tail.

Calibration is equivalent: on decoys, both models' `P(SpecEValue < x)` deviate from the ideal `x`
by the same factor at every threshold — our null is neither more nor less optimistic than UC's.

**Cost:** none, at this configuration — rescoring 7,431 PSMs over 1,252 F13 spectra takes 4.41 s
with the UC model and 4.40 s with ours (1686 vs 1690 PSM/s), because both score ≈3.7–3.9 ion types
per partition. The 0.15-threshold variant, with 5.3, is the one that costs ~7 %.

## Known limitations / what a v1 needs

1. **Precursor filtering can't be learned from library spectra.** MassIVE-KB consensus spectra have
   had the precursor region removed: 0.000 of charge-2 spectra retain a peak within ±0.5 Da of the
   precursor m/z, so the counted §4 table is empty. The trainer falls back to a chemistry-derived
   table (precursor ±1 nominal bin and the water loss, per charge) — a filtering *rule*, not a
   trained number, since `preprocess` uses only `(reduced_charge, offset)` and ignores the
   frequency. Measured effect on F13 was ~1 ID, so it is not the source of the residual gap.
2. **The optimistic-completeness gap is the source**, and the fix is corpus, not code: train from
   the raw acquisitions MassIVE-KB's provenance links to (`docs/models.md` §4.2b, milestone M4).
   That needs an mzML reader in `msgf-io`, which does not exist yet — the corpus route is the
   blocker, not the trainer.
3. **Coverage.** Only the HCD/HighRes/Tryptic identity is trained — it is the bundled default, and
   any other acquisition still needs `--param` (MS-GF+'s models remain readable, just not shipped). CID/ETD need c/z candidates
   added to `ions.rs` and a matching corpus; TMT/iTRAQ and non-tryptic need corpora.
4. **No trainer-mechanics oracle.** `plans/PLAN1.md` 5c-1 wants Java `ScoringParamGen` and `msgf-train`
   run on the same small corpus and compared table-for-table. Not done: our statistics are defined
   independently (deliberately — see the clean-room boundary), so the two would not be expected to
   agree bin-for-bin; the meaningful gate is the downstream parity measured above.
