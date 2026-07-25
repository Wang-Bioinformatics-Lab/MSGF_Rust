# The bundled scoring model

`MSGFRust_HCD_HighRes_Tryp_v1.param` — the fragment-scoring model MSGF_Rust ships and uses when no
`--param` is given. It is the artifact that lets this project be MIT: **no byte of it comes from
MS-GF+**.

| | |
|---|---|
| identity | HCD / HighRes / Tryp / Automatic |
| size / SHA-256 | 1,067,562 B / `5f6ab76f5f849609f9901379536b4db98d93f355c5b428de108e3ad41e432d02` |
| partitions | 236 — precursor charges 2–6 × 2 m/z segments |
| ion types | y, y+¹³C, b, y−NH₃, y−H₂O, y+2¹³C, a, b−H₂O, b+¹³C (3.69 per partition) |
| trained from | 258,352 PSMs (281,855 spectra read) |
| corpus | MassIVE-KB peptide spectral libraries, MassIVE **MSV000081142** — 20 library shards |
| corpus license | **CC0 1.0** (public domain dedication) |
| trained by | `msgf-train` (this repo, `rust/crates/msgf-train`) |
| training time | ~12 s wall (32 cores), counting only |

## Reproducing it

```bash
cd validation && ./fetch_reference_data.sh --training 20     # ~950 MB, CC0
cd ../rust && cargo build --release -p msgf-train
./target/release/msgf-train \
  $(for f in ../validation/data/training/*.mgf; do echo -n "--corpus $f "; done) \
  --ion-threshold 0.25 \
  --out crates/msgf-scorer/models/MSGFRust_HCD_HighRes_Tryp_v1.param
```

Training is a counting pass — no optimiser, no RNG — so this reproduces the bytes exactly. The
corpus shards it was built from are pinned in `MSGFRust_HCD_HighRes_Tryp_v1.corpus.sha256`.

Everything else about the trainer (what each table counts, the clean-room position, the
measured quality) is in [`docs/training.md`](../../../../docs/training.md).

## How it compares to MS-GF+'s `HCD_HighRes_Tryp.param`

Both models rescoring the same PSM list (identified peptide + 5 mass-identical shuffled decoys per
spectrum):

| benchmark | MS-GF+ model | this model |
|---|---|---|
| held-out MassIVE-KB (1,500 spectra, ground-truth peptides) — IDs at a 1 %-decoy threshold | 1485 | **1480** |
| held-out MassIVE-KB — median target log₁₀ SpecEValue | −14.65 | −14.47 |
| F13 raw spectra (1,252 spectra, out-of-domain) — IDs at a 1 %-decoy threshold | 66 | **63** |
| F13, confident fifth of the PSMs — IDs at a 1 %-decoy threshold | 20 | **24** |
| F13 — Spearman ρ vs MS-GF+ on log₁₀ SpecEValue | — | 0.94 |
| rescoring throughput | 1686 PSM/s | 1690 PSM/s |

It independently rediscovers the same ion types MS-GF+ scores, in the same order of prevalence.
The residual difference is corpus domain, not code: MassIVE-KB reference spectra are *consensus*
spectra, so a real ion series is more complete there (y-ion hit rate 0.88) than in raw acquisitions
(0.66), which makes the model's missing-ion penalty slightly harsh on raw data. Training from
provenance-linked raw spectra is the planned v2 (`PLAN1.md`, milestone M4).

## What it is not: a way to reproduce MS-GF+'s output

The table above says this model is **about as good** as MS-GF+'s. It does not say it is the
**same**, and those are different claims. This repo's bit-exactness contract — RawScore and
DeNovoScore exact, SpecEValue within `|log10| ≤ 0.05` — is measured with **MS-GF+'s own
`HCD_HighRes_Tryp.param` loaded via `--param`**. It does not carry over to this model, and cannot:
a scoring model *is* the scoring function, so different tables mean different scores, hence
different top peptides and different SpecEValues.

Measured on the F13 iPRG-2013 set (1,406 spectra vs the concatenated human DB, `-t 10ppm -ti 0,1`,
oxidation on M — the same run in `PERFORMANCE.md`), top peptide per scan vs MS-GF+'s own output:

| Rust run | same top peptide as MS-GF+ | target PSMs |
|---|---|---|
| `--param HCD_HighRes_Tryp.param` (MS-GF+'s model) | **92.7 %** (1161/1253) | 624 |
| bundled model (no `--param`) | **66.6 %** (834/1253) | 609 |

Swapping only the model changes the reported peptide on ~1/3 of scans (68.1 % agreement between
the two Rust runs). Read that as a *divergence* measurement, not a quality verdict: F13 identifies
essentially nothing (MS-GF+ itself finds one target PSM at 1 % FDR), so most of those scans are
noise, and on noise the top hit is close to arbitrary — any two models will disagree. Agreement
does climb with confidence (on MS-GF+'s 1e-9–1e-6 tier, n=373: 95.2 % with its model, 75.3 % with
this one), but F13 has too few confident PSMs to resolve the high-confidence end. The quality
evidence is the held-out benchmark table above, not this one.

Practical rule:

- **Reproducing or diffing against MS-GF+ output** → pass `--param HCD_HighRes_Tryp.param`. The
  default will not match, and that is expected.
- **Running MSGF_Rust as its own search engine** → the bundled model is the right default: MIT-clean
  and benchmarked competitive on held-out ground truth.

Either way the model in use is printed on stderr at the start of every run, so a result is always
traceable to the tables that produced it.

## Changing it

Retraining is a deliberate, visible act: update the bytes, the `SHA256`/length constants in
`../src/bundled.rs`, this table, and `docs/training.md` together. `bundled_model_bytes_are_pinned`
fails until you do.

## Other identities

Only HCD / high-resolution / tryptic is trained. For CID, ETD, QExactive-tuned or non-tryptic runs,
pass `--param` with a model for that identity (including MS-GF+'s own, which remain readable — they
are simply not distributed here).
