# PLAN2 — Target–decoy and FDR

Execution plan for target-decoy analysis (TDA): decoy database construction, PSM/peptide-level
q-values, and how they wire into the future `msgf-search` engine. Strategy context is `PLAN.md`
(§7 Phase 6 lists "target-decoy FDR" as a single bullet); this doc is the concrete design.

**Status: TD-1, TD-2 (both gates) and TD-3 implemented** (`msgf-db`, `msgf-fdr`, `msgf-search`).
The decoy writer is byte-identical to both reference `.revCat.fasta` files; the q-value columns
reproduce the F13 golden exactly (1610/1610); and the Gate 2 Java probe now pins 110 thresholds
and 556 lookups across 14 synthetic cases, which **corrected three rules** F13 could not see (see
§1.4 and §3). **Still open:** the §4 oracle problem — there is no benchmark on which "ID counts at
1 % FDR" can be compared to Java at all. What we *do* already have is an oracle — the committed
`validation/golden/iprg2013_F13.golden.json` carries `protein`, `qvalue` and `pep_qvalue` per PSM,
and `fetch_reference_data.sh --full` puts Java-generated decoy FASTAs on disk
(`iprg2013_human.revCat.fasta`, `Tryp_Pig_Bov.revCat.fasta`).

Two of the three deliverables below need **no search engine** and can land now.

---

## 1. What MS-GF+ does (normative reference)

Line numbers are `github.com/MSGFPlus/msgfplus` @ `v2024.03.26`, the release our jar comes from.

### 1.1 Decoy database — `msdbsearch/ReverseDB.java:34`

`-tda 1` builds `<db>.revCat.fasta` (`MSGFPlus.java:216-230`) and searches **that** — one
concatenated pass, never a separate decoy search.

- Each protein is reversed **whole-sequence** (not peptide-wise), header `>XXX_<original header>`.
- Prefix default `XXX` (`MSGFPlus.java:29`); a user prefix has trailing `_` stripped, then `_` is
  re-appended, so `-decoy REV_` and `-decoy REV` both yield `REV_`.
- Byte details that matter for a byte-exact writer: the target block is copied **line-by-line,
  wrapping preserved**; each decoy protein is emitted as **one unwrapped line**, `.trim()`ed; the
  final record is flushed after the read loop.
- Load-time sanity gates (`MSGFPlus.java:238-252`): unique-protein ratio ≥ 0.5, decoy fraction in
  [0.4, 0.6], else hard exit.

### 1.2 Scoring is unchanged by TDA

Decoy peptides go through the identical scorer and generating function. One thing to note:
background amino-acid frequencies are counted over the **concatenated** FASTA
(`DBScanner.java:919`). Reversal preserves per-residue counts exactly, so every count and the total
exactly double and the `float` quotient is bit-identical to the target-only value — our existing
`--aa-probs` path needs no change.

### 1.3 What makes a PSM a decoy — `fdr/MSGFPlusPSMSet.java:61-77`

A match is decoy iff **every** protein occurrence's accession starts with the prefix. A peptide
shared between a target and a decoy protein counts as **target**.

### 1.4 q-values — `fdr/TargetDecoyAnalysis.java:160-261`, `fdr/ComputeFDR.java:247`

- Score is **SpecEValue**, smaller-is-better, `float` (**f32**) throughout.
- Two populations: PSM-level = one entry per reported match (`considerBestMatchOnly=false`, i.e.
  all `-n` matches per spectrum); peptide-level = peptide → **best (min)** SpecEValue. Peptide key
  is the mod-bearing sequence with flanks stripped and uppercased (`TSVPSMSet.java:220`).
- FDR sweep: sort both lists ascending; for each **distinct** decoy score advance a target cursor
  past all targets strictly better; if `targetIndex <= decoyIndex` → FDR 1, else
  `round(decoyIndex · pit) / targetIndex` (Käll et al. 2008, D/T; `pit` is always 1); clamp to 1;
  **break** once FDR ≥ 1; ±∞ sentinels are seeded first.
- q-value conversion: running minimum walking from worst key to best.
- Lookup is `higherEntry(score)` — the value at the **strictly next-larger** key — which is why a
  score exactly equal to a map key resolves to the next-worse threshold's q-value.
  - **Settled by the Gate 2 probe.** An earlier reading of F13 alone suggested a floor lookup
    (largest key ≤ score) instead; `DumpFdrMap.java` shows that is wrong. Probing each threshold's
    immediate float neighbours, MS-GF+ answers `q(1e-7) = 0.25` for the threshold `1e-7 → 0.0`,
    i.e. the **next** threshold's value — `higherEntry`, as the source says. On F13 the two rules
    happen to agree (§4: only two distinct q-values), which is exactly why the synthetic probe was
    needed. `msgf-fdr` now implements `higherEntry`, with **both** `±∞` sentinels seeded.
- Reported as `QValue` (MS:1002054) and `PepQValue` (MS:1002055).

**Verified:** a Python transcription of exactly the above reproduces both columns for all 1610
unique PSMs of the F13 golden, bit-for-bit in f32. The algorithm description here is not a guess.

---

## 2. Decisions needed

**D6 — decoy construction.** (a) MS-GF+-compatible whole-protein reversal; (b) peptide-level
reversal keeping the C-terminal residue fixed (Sage/Comet style, better-matched decoys for tryptic
search); (c) shuffling. *Recommend (a) as the default* — it is the oracle path and the only one
byte-comparable to Java — behind a `DecoyStrategy` trait so (b)/(c) can be added without churn.

**D7 — FDR formula.** MS-GF+ uses Käll `D/T`. Elias–Gygi `2D/(T+D)` is in the Java source as a
commented-out alternative. *Recommend* MS-GF+ default, `--fdr-formula` flag for the other.

**D8 — arithmetic width.** *Recommend f32* on the compatibility path: it is what the oracle emits
and q-values are then comparable by exact equality, not tolerance. f64 optional behind a flag.

**License note:** decoy generation and the FDR sweep are textbook published methods with no UC
trained data involved, so this whole workstream is MIT-safe by construction. Only *behavioral*
compatibility comes from reading the Java — the same posture as the rest of the fidelity work, and
unrelated to the `.param` model boundary in `PLAN1.md`.

---

## 3. Deliverables

### TD-1 — `msgf-db` crate + `msgf decoy` subcommand *(independent, ~1 day)*

New crate `msgf-db` (FASTA today; the natural home for digestion and the fragment index later):

- streaming FASTA reader/writer, no dependencies;
- `DecoyStrategy` trait + `ReverseProtein` impl reproducing §1.1 byte-for-byte, including the
  wrapping asymmetry between the target and decoy blocks;
- `msgf decoy --fasta in.fasta -o out.revCat.fasta [--prefix XXX] [--concat|--decoy-only]`;
- the two load-time sanity checks as a reusable `validate_concatenated()`.

*Gate:* output **byte-identical** to `validation/data/fasta/Tryp_Pig_Bov.revCat.fasta` and (with
`--full` data) `iprg2013_human.revCat.fasta`. Fetched data → the test must use the
`eprintln!("skip: ...")` pattern from CLAUDE.md.

### TD-2 — `msgf-fdr` crate + `msgf fdr` + `rescore --fdr` *(independent, ~1–2 days)*

- `PsmRecord { spec_key, peptide, proteins, score }`; `TargetDecoyAnalysis::new(psms, prefix, pit)`
  with `psm_qvalue(f32)` / `pep_qvalue(&str)`; f32 internals; the FDR map as a sorted
  `Vec<(f32, f32)>` with a binary search reproducing `higherEntry` semantics.
- `msgf fdr --psms x.tsv --decoy-prefix XXX -o y.tsv` appending `QValue`/`PepQValue`, and
  `rescore --fdr` doing the same over freshly computed SpecEValues. This is real user-facing value
  before the search engine exists: rescore a target+decoy PSM list and get q-values out.

*Gate 1 (data-free):* reproduce `qvalue` and `pep_qvalue` from `iprg2013_F13.golden.json` exactly.
The golden is committed, so this runs on a fresh clone. Requires rolling the `-unroll 1` TSV rows
(4133) back up into one record per match (1610) — MS-GF+ counts *matches*, not protein occurrences.

*Gate 2 (the real coverage): **done**.* `validation/reference/java/DumpFdrMap.java` (run by
`make_fdr_golden.sh`, JVM + jar only — no spectra, models or database) drives
`TargetDecoyAnalysis.getFDRMap` **and** `getPSMQValue` over 14 cases: ties within and across the
lists, an empty decoy list, an empty target list, all-decoys-better, single-element lists, the
`target_index == 0` skip, the early `break` at FDR ≥ 1, a non-monotone raw sweep, and three seeded
pseudo-random sets. Every threshold's immediate float neighbours are probed, which is what pins the
lookup rule. Frozen to `validation/golden/fdr/fdrmap_cases.golden.json` (committed) and checked
bit-for-bit by `msgf-fdr/tests/golden_fdrmap.rs` — 110 thresholds, 556 lookups — plus a no-JVM
re-derivation in `validation/regression/run_regression.py`.

Three rules the F13 gate could not see were **wrong** before this probe, each of which silently
mis-reports q-values on data where targets and decoys actually separate:

1. a run of equal decoy scores was charged at its *last* index instead of its first (Java takes
   `decoy_index` from the run's first member, so `tie_within_decoys` is FDR 0/2, not 1);
2. a threshold with no better target ended the sweep instead of being skipped, so
   `guard_no_target_better` reported q = 1 everywhere where MS-GF+ reports 1/3;
3. the lookup was a floor instead of `higherEntry` (§1.4).

### TD-3 — wire into `msgf-search` *(blocked on the engine)*

The engine must index the concatenated FASTA, keep each peptide's protein-index list, label decoys
by accession prefix per §1.3, retain top-N matches per spectrum, and hand a `Vec<PsmRecord>` to
`msgf-fdr` **once, after all spectra** — FDR is global, so it is a serial epilogue to the parallel
search. Output: MS-GF+'s TSV column set first, mzIdentML with MS:1002054/55 later.

*Gate:* PSM and peptide counts at 1% FDR within a stated margin of Java MS-GF+ — on a benchmark
that actually identifies peptides, which we do not currently have (§4).

Ordering: TD-1 is a prerequisite for the engine's index build anyway; TD-2 is fully independent and
can be done first or in parallel.

---

## 4. The oracle problem: F13 cannot validate FDR

Measured from the committed golden: 4133 unrolled rows → **1610 unique PSMs**, of which **765 are
decoy-only** vs 845 target. MS-GF+'s own `QValue` is **1.0 for 4132 of the 4133 rows** (exactly one
row at 0.0). Top-scoring hits are R/K-rich low-complexity junk (`R.RRSTRSEELTR.S`,
`R.RRKNKLKR.R`, …) and targets and decoys interleave from the second-best hit onward.

Consequences:

- F13 remains a perfectly good **scoring** oracle — per-PSM RawScore/DeNovoScore/SpecEValue match
  Java regardless of whether the identifications are biologically real. Nothing already validated
  is affected.
- It **cannot** support Phase 6's stated gate ("ID counts at 1% FDR comparable to Java MS-GF+",
  `PLAN.md:250`), nor any end-to-end TDA test: with q ≡ 1 there is nothing to compare.
- Both the spectra and the FASTA came from msgfplus' own `src/test/resources/iprg-2013/`
  (`fetch_reference_data.sh:39,60`) and the run used `-inst 1 -m 3 -e 1 -t 10ppm -tda 1`, so the
  pairing is presumably intended; why the search identifies essentially nothing is **open**.
  - **Partial answer (measured):** the run inherits MS-GF+'s default of **unlimited missed
    cleavages** (`-maxMissedCleavages` is not passed). Only 27% of its own top-hit sequences are
    reachable with ≤2 missed cleavages; the rest are long K/R-rich low-complexity peptides that a
    conventional `-c 2` search would never consider. Reproducing the golden requires the unlimited
    setting (`msgf-db`'s `UNLIMITED_MISSED_CLEAVAGES`, the default in `DigestParams`), and the
    permissiveness is a large part of why the top hits are junk. Whether the spectra are simply
    not human tryptic digest remains open.
  Resolve it before TD-3's gate — either fix the search configuration or adopt a different
  public spectra+DB benchmark for Phase 6.

---

## 5. Goldens and the data contract

- **Done:** `validation/golden/fdr/fdrmap_cases.golden.json` (Java probe, committed, JVM only at
  generation time), generated by `validation/reference/make_fdr_golden.sh` (also run by
  `build_all_golden.sh --with-java`) and re-derived without Java by `_fdr_map()` in
  `validation/regression/run_regression.py`.
- Edit: add `qvalue`/`pep_qvalue` to the `compare` block of `iprg2013_F13.golden.json` with
  `assert: exact` — the fields are already present in the file, so this is a `parse_msgf_tsv.py`
  change, not a golden regeneration.
- The decoy-FASTA byte gate depends on fetched data and must skip gracefully.

## 6. Non-goals

Protein inference and protein-level FDR; Percolator/LDA-style rescoring (Sage does this, MS-GF+
does not); entrapment or two-species FDR; π0/PIT estimation beyond MS-GF+'s fixed `pit = 1`;
group-specific or subset FDR.
