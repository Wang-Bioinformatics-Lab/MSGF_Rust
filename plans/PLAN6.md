# PLAN6 — timsTOF (Bruker `.d`) DDA support

Execution plan for reading **Bruker timsTOF `.d` folders directly** and searching DDA-PASEF data
end-to-end: `msgf search -s run.d -d db.fasta`, with no msconvert step in between.

**Status: design doc, not started.**

*(Numbered PLAN6 because PLAN4 is the desktop-UI plan and PLAN5 is reserved by PLAN4 §5.1.)*

This is only half an I/O task. The other half is that MS-GF+'s scoring is conditioned on an
instrument through the `.param` model, and the fragment-matching tolerance the scorer uses comes
from **that model alone** (`msgf-scorer/src/scored_spectrum.rs:205` → `model.mme`), with no CLI
override. An Orbitrap-trained model applied to TOF fragments is the difference between "it runs"
and "it identifies peptides". §6 is not follow-up work; it is what makes §4 worth shipping.

---

## 1. Goal and scope

**Success criterion.** A public DDA-PASEF run searched directly from its `.d` folder produces a PSM
TSV whose identification rate at 1% FDR is competitive with the same run converted to mzML/MGF and
searched by MS-GF+ — measured on ground truth (§5, Layer 3), not on agreement with MS-GF+.

**In scope**

- a `.d` reader producing the existing `msgf_io::Spectrum` (`msgf-io/src/lib.rs:20`),
- DDA-PASEF precursor → merged-MS2 assembly, and the precursor metadata that comes with it,
- the spectrum-model fields TIMS needs and MGF has no place for (ion mobility, RT, precursor id),
- format dispatch in `msgf-io` and in the CLI's `-s` (`msgf-cli/src/search.rs:335`),
- fragment tolerance: a CLI override, and a timsTOF `.param` trained with `msgf-train`,
- validation that respects the data-absence contract and does not overclaim bit-exactness (§5).

**Out of scope** — diaPASEF/DIA, `.tsf`/MALDI/imaging, MS1 feature finding or 4D quant,
ion-mobility-aware rescoring or CCS prediction (emitting a mobility column is in scope; *scoring*
with it is not), chimeric-spectrum deconvolution (§7), mzML (orthogonal — but T1 keeps the door
open).

**Non-goals.** Changing scoring arithmetic. Making the `.d` path a fidelity oracle. Silently growing
the default dependency tree or breaking the single-static-binary release story (T1).

---

## 2. What TIMS data is, operationally

The facts below drive every decision in §3. None of them are true of MGF.

| Fact | Consequence here |
|---|---|
| A `.d` folder is `analysis.tdf` (SQLite) + `analysis.tdf_bin` (compressed frame blobs) | Reading is a dependency decision, not a parser we write (T1) |
| Frames are 2D: m/z × ion-mobility scan | There is no "spectrum" in the file; one has to be assembled |
| In DDA-PASEF a single precursor is fragmented across several frames / scan ranges | A usable MS2 is a **sum over windows** — so centroiding and merging become part of *our* numbers (T2) |
| `Precursors` gives `MonoisotopicMz` (nullable), `Charge` (nullable), intensity, mobility, parent frame | Charge/mono fallbacks matter far more than with MGF. Absent charge is already handled via `charge_range` (`msgf-search/src/search.rs:264-272`); absent precursor m/z is a hard skip |
| Mono-isotope picking is instrument-side and imperfect | `-ti 0,1` earns its keep — measure how often it is what saves the ID |
| Merged PASEF MS2 peak lists are dense; intensities are summed integer counts | `preprocess()` cost *and* the `probPeak` estimate (`scored_spectrum.rs:142`, `|peaks| / (peptideMass / 2·mme)`) are peak-density sensitive. Measure before assuming (§7) |
| TOF fragment accuracy is ppm-scale and wider than Orbitrap | The model's `mme` window is the whole ballgame (T4, §6) |
| MS-GF+ reads mzML/mzXML/MGF/ms2 — **not `.d`** | There is no direct Java oracle for this path. §5 splits the contract accordingly |

---

## 3. Decisions needed

**T1 — reader dependency, and the release-binary constraint.** Two credible routes, both Apache-2.0
(permissive, MIT-compatible; carry the NOTICE — add to `LICENSING.md`). Neither is UC-licensed, so
**neither touches the clean-room or model-ownership story.**

| Route | Gets us | Costs |
|---|---|---|
| `timsrust` directly (0.6.3, MannLabs) | current API, DDA precursor→spectrum assembly built in, leanest tree | `.d` only; mzML still unsolved |
| `mzdata` + `bruker_tdf` | mzML **and** `.d` in one dep; would settle `PLAN.md` D3 too | pins `timsrust` 0.4.1; pulls `rusqlite` + `parking_lot` + `mzsignal` |

**Recommend `timsrust` directly**, behind an **off-by-default `tims` feature on `msgf-io`**, with
both readers behind one `SpectrumSource` abstraction so revisiting this when mzML lands is a swap,
not a rewrite. Pin the version exactly — `timsrust` went 0.5.0 → 0.6.3 in about six weeks, and a
reader upgrade changes our peak lists.

**This collides with a stated repo posture and the collision has to be decided, not finessed.** The
workspace has exactly one runtime dependency (`rayon`) across ten crates, `release.yml` ships one
static `msgf` binary for five targets, and PLAN4 §2 commits to keeping both. Either route drags in
`rusqlite` — i.e. **C SQLite**, plus zstd — which needs a working C toolchain for every
cross-compiled target. So:

- keep `tims` **off by default even for `msgf-cli`** until cross-compilation is proven green on all
  five targets in a scratch branch;
- if a target cannot be made to build, ship `.d` support as a **separate artifact** rather than
  dropping a target or vendoring a pure-Rust SQLite reader on a whim;
- that experiment is **M1's gate**, not an afterthought — it can invalidate the whole approach, so
  run it first, before TIMS-1.

**T2 — assembly semantics.** Use `timsrust`'s DDA spectrum reader (one centroided spectrum per
precursor) rather than assembling from frames ourselves. Its merging/centroiding behavior is then
part of our output: record the crate version and settings in run provenance so a change in the
numbers is attributable to it.

**T3 — spectrum identity.** The report columns (`msgf-search/src/report.rs:11`) are MS-GF+'s and
assume scan numbers; `msgf fdr` and `rescore --psms` join on them. Decide once and lock with a test:
`ScanNum` = precursor index, `SpecID` = a string that round-trips to the raw data (precursor id +
parent frame). A user must be able to get from a PSM row back to the `.d`.

**T4 — fragment tolerance.** Today the only knob is the model's `mme`, and `msgf-train`'s CLI
exposes it as **Da only** (`msgf-train/src/bin/train.rs:38`, default 0.5 Da) even though the format
and the decoder both support ppm (`msgf-scorer/src/lib.rs:231`, `write.rs:90`). Do both:

1. add `msgf search --frag-tol <10ppm|0.02Da>` overriding `model.mme` — cheap, unblocks measurement,
   and useful well beyond TIMS;
2. accept ppm in `msgf-train --mme`, and verify the write/read round-trip preserves the unit.

Document that `--frag-tol` changes the scoring function and sits outside the bit-exactness contract.

**T5 — mobility in the output.** Emit `IonMobility` (1/K₀) and `RT`? Appending columns risks
strict-header consumers. Recommend appending at the end, only when the source carries them, and
updating `report.rs`'s header test in the same commit.

---

## 4. Deliverables

### TIMS-0 — spectrum model *(independent, ~0.5 day)*
Add `rt_seconds`, `ion_mobility`, `precursor_intensity`, `source_id: Option<String>` to `Spectrum`.
All `Option`, `Default` preserved, so the existing golden spectra JSON is unaffected. Fill
`rt_seconds` from MGF `RTINSECONDS=` while here — it is currently dropped
(`msgf-io/src/lib.rs:121`).

### TIMS-1 — `msgf-io` `tims` feature *(~2 days, gated on M1)*
`TimsReader` over a `.d` path → `Iterator<Item = io::Result<Spectrum>>`. Charge fallback (missing →
leave `None` and let `charge_range` drive), mono-m/z fallback, and peaks guaranteed ascending in m/z
so `preprocess()` keeps file order (`preprocess.rs:136`). Unit tests against the synthetic fixture
(§5) — zero fetched bytes.

### TIMS-2 — dispatch *(~0.5 day)*
`msgf_io::read_spectra(path)` dispatching on extension/dir shape; `-s` accepts `.d` in `search` and
`rescore`. A clear error when the binary was built without the feature.

### TIMS-3 — `--frag-tol` + ppm `--mme` *(~0.5 day)* — per T4.

### TIMS-4 — fixtures and goldens *(~2 days)* — see §5.

### TIMS-5 — a timsTOF `.param` *(~3 days + corpus time)* — see §6.

### TIMS-6 — docs *(~0.5 day)*
`README.md` quickstart, the `msgf-io` bullet in `CLAUDE.md`, the dependency note in `LICENSING.md`,
the `plans/README.md` row, and `docs/models.md` for the new model.

---

## 5. Validation — what the oracle can and cannot be

**MS-GF+ cannot read `.d`.** So the contract splits into three layers, and conflating them is the
main way this plan could produce a dishonest green suite.

**Layer 1 — unchanged, bit-exact.** Everything downstream of the reader is still held to the
existing contract. Test: dump our `.d`-assembled spectra to MGF, run the MS-GF+ jar on *that MGF*
with the same `.param`, and require exact RawScore/DeNovoScore plus `|log10(rust/java)| ≤ 0.05` on
SpecEValue. This is a real golden, and it validates the whole scoring path on TIMS-shaped input. It
is MS-GF+-derived → **gitignore it and wire it into `build_all_golden.sh --with-java`**, or its test
skips forever on a fresh checkout.

**Layer 2 — reader agreement, tolerance-based.** Our `.d` → spectra vs msconvert's `.d` →
mzML/MGF. **Do not assert bit-exactness here** — centroiding and PASEF merging legitimately differ.
Gate on: precursors matched, per-precursor peak-count ratio, m/z agreement within a stated ppm for
the top-N peaks, and — the one that actually matters — PSM overlap at 1% FDR.

**Layer 3 — quality, on ground truth.** Per `CLAUDE.md`'s standing rule: a held-out annotated
timsTOF shard (`SEQ=` peptides + mass-identical shuffled decoys) through
`validation/eval_trained_model.py library`. Agreement with MS-GF+ is not the gate.

**Fixtures.** Write a synthetic `.d` generator in `validation/reference/` (stdlib `sqlite3` plus the
frame encoding) so TIMS-1's unit tests need **zero fetched bytes** — the same discipline as
`author_a_model_from_scratch`. `timsrust`'s own `tests/test.d` (Apache-2.0, simulated, tiny) is the
fallback if writing the encoder proves disproportionate. A real public DDA-PASEF run goes behind
`fetch_reference_data.sh --tims` for the end-to-end.

**Data contract.** A golden derived from public open-licensed data plus *our own* reader owes
nothing to MS-GF+ and may be committed (like `chemistry/`). Anything produced by running the jar may
not. Every new golden test follows the skip-if-absent pattern.

---

## 6. The model problem

The `.param` is instrument-conditioned in three places that all matter here: the fragment tolerance
`mme`, the fragment-offset and rank distributions, and `ErrorDist` — which is literally a learned
**mass-error histogram** (`scored_spectrum.rs`'s `error_score`, quantized by
`error_scaling_factor`). An Orbitrap histogram is the wrong prior for TOF, and a too-tight `mme`
discards true fragment matches before scoring ever sees them.

Deliverable: **`MSGFRust_TIMS_DDA_Tryp_v1.param`**, trained by `msgf-train` from an annotated
DDA-PASEF corpus. Same machinery as the bundled MassIVE-KB model, same MIT-clean provenance — this
*extends* PLAN1's model-ownership story rather than competing with it.

Config deltas to determine empirically, not assume: `mme` (ppm, TOF-scale), `error_scaling_factor`,
segment/partition boundaries, and the charge histogram (PASEF skews high).

**Gate.** On a held-out timsTOF shard the TIMS model must beat the bundled HCD model. If it does
not, do not ship a second model — record the measurement and keep one.

---

## 7. Risks and open questions

- **Cross-compilation (T1).** The sharpest risk, and the reason M1 is ordered first. C SQLite in a
  five-target static release is where this plan is most likely to break.
- **Chimeric MS2.** PASEF co-isolates; MS-GF+'s model assumes one precursor per spectrum. Expect a
  lower ID rate than Bruker-native engines. Out of scope to fix — but measure it, so the README
  claim is honest.
- **Peak density.** `probPeak` and `preprocess()` were tuned on Orbitrap densities. Measure both
  before adding a pre-filter — a pre-filter is a scoring change, not a perf tweak.
- **Reader churn.** A `timsrust` bump can change peak lists. Treat it as a numbers-changing event
  requiring a Layer-2 re-run, never a routine `cargo update`.
- **Missing charge / mono rate.** Unknown until measured; determines how much work `-ti` and
  `charge_range` are really doing.
- **Open:** does a public DDA-PASEF dataset exist with annotations good enough for Layer 3, or must
  the corpus come from a MassIVE-KB timsTOF shard? Resolve before committing to TIMS-5.

---

## 8. Milestones

| M | Content | Done when |
|---|---|---|
| M1 | **T1 cross-compilation spike** | `timsrust` builds on all five `release.yml` targets, or the separate-artifact decision is made and written down |
| M2 | TIMS-0,1,2 + synthetic fixture | `msgf search -s run.d` reads and reports spectrum/peak counts; unit tests pass with zero fetched bytes |
| M3 | TIMS-3 + Layer-1 golden | a real DDA-PASEF run searches end-to-end; bit-exactness through the new front end is proven, not assumed |
| M4 | TIMS-5 + Layer-3 numbers | the TIMS model beats the HCD model on a timsTOF hold-out, or is dropped with the measurement written down |
| M5 | TIMS-6 + Layer-2 report | `.d` in the README quickstart; msconvert-comparison numbers published |

---

## 9. Rules this plan does not get to break

1. The reader must not touch the DP, the scorer's arithmetic, or preprocessing order.
2. No golden asserting bit-exactness against msconvert, or against MS-GF+ using a model trained here
   (`CLAUDE.md`'s standing rule — and TIMS-5 makes it newly tempting).
3. `--frag-tol` and the TIMS model sit explicitly outside the bit-exactness contract; the MS-GF+
   fidelity tests keep passing MS-GF+'s own `.param`.
4. Any MS-GF+-derived golden added here is gitignored **and** wired into `build_all_golden.sh`.
5. The default library build stays dependency-light, and the release story stays intact or changes
   deliberately (T1).
