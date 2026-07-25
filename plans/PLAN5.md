# PLAN5 — Nextflow workflow for parallel scale-out

Execution plan for running one MSGF_Rust search **across many machines** with Nextflow:

```bash
nextflow run mwang87/MSGF_Rust --spectra 'data/*.mgf' --fasta human.fasta -profile slurm,singularity
```

**Status: design doc, not started.** Nothing Nextflow-related exists in the repo today (no `.nf`, no
`Dockerfile`, no `nextflow.config`).

**Standalone only.** This is a generic pipeline anyone can run — local, SLURM, or Kubernetes. It is
deliberately **not** coupled to GNPS2; no `workflowinput.yaml`, no GNPS2 result views. A GNPS2
wrapper over this core would be a separate, later plan.

---

## 1. Goal and success criteria

Search-per-spectrum is embarrassingly parallel; FDR is not. The plan is to exploit the first without
compromising the second.

**Primary success criterion (the identity gate):** for the same inputs, a distributed run with *N*
chunks produces **the same PSMs and the same q-values** as a single-process run, for every *N*.
Chunking is a scheduling decision and must have **zero** numerical consequence.

**Secondary:** near-linear speedup in chunk count until the per-task index build dominates (§5), and
one self-contained container image with no external model download.

**Non-goals.** Any GNPS2 coupling. mzML input (`msgf-io` is MGF-only). A distributed or shared
peptide index. Splitting the protein database — §2 explains why that is not a performance decision
but a correctness one. Changing any scoring or FDR code: this plan adds *plumbing* to the CLI and a
pipeline around it, nothing more.

---

## 2. The constraint that dictates the topology

Two facts about SpecEValue decide the entire shape of this pipeline.

1. **The null model is database-derived.** `msgf search --help` states it: "Amino-acid background
   frequencies for the generating function are computed from the database being searched, so
   SpecEValue reflects that database rather than a uniform alphabet."
2. **The E-value scale is database-sized.** `--db-size` defaults to the candidate-index size, and
   `EValue = SpecEValue × db_size`.

Therefore:

> **Scatter on spectra. Never on the database.**

Splitting the FASTA across tasks would give each shard a *different* background frequency vector and
a *different* E-value multiplier. The resulting SpecEValues would not be comparable, the merged score
ranking would be meaningless, and the target-decoy FDR computed over it would be silently wrong —
not noisy, wrong. This is the one thing in this plan that must never be "optimized."

The complementary fact is that **FDR is global by construction**. `assign_q_values` is documented as
requiring the whole result set at once, and is deliberately serial (`msgf-search/src/lib.rs`,
`PLAN2.md` §TD-3). That makes it the natural gather step — and `msgf fdr` already does exactly this
job over a merged PSM table, including the correct roll-up of `-unroll`-style repeated rows.

```
  FASTA ──► PREPARE_DB ──────────────┐  (decoys built ONCE)
                                     ▼
  MGF ────► SPLIT ──► chunk_1 ──► SEARCH ──┐
                      chunk_2 ──► SEARCH ──┤
                       ...      ──► ...    ├──► MERGE ──► FDR ──► REPORT
                      chunk_N ──► SEARCH ──┘   (concat)  (global   (IDs at
                         └── parallel ──┘                q-values)  1% FDR)
```

---

## 3. Blockers in the CLI today

The pipeline cannot be correct until these land. Each is small, additive, and independently testable
— and each is a genuine defect that a chunked run would otherwise hide.

### 3.1 `SpecID` collides across chunks — **must fix**

`report.rs` writes `index={p.spec_index}`, and `spec_index` comes from `enumerate()` over the
spectra handed to `run` (`msgf-search/src/search.rs:239`). It is **chunk-local**: every chunk restarts at
`index=0`. Concatenating chunk outputs therefore produces many rows claiming `index=0`, and any
downstream tool keyed on `SpecID` silently mis-joins.

Fix: `msgf search --spec-index-offset <N>`, added to `spec_index` at report time. The splitter emits
each chunk's starting offset in a manifest, and the pipeline passes it through.

### 3.2 `#SpecFile` loses provenance — **must fix**

Each chunk reports its own chunk filename, so the merged table names files that were never the
user's input. Fix: `msgf search --spec-file-name <NAME>` to override the reported name, set to the
original MGF.

### 3.3 `msgf fdr` takes one input — **should fix**

Merging is otherwise an `awk` step to strip repeated headers, which drags a shell-text dependency
into a pipeline whose entire appeal is one static binary. Fix: let `-i/--psms` repeat, concatenating
inputs (validating that headers agree) before the sweep. The MERGE process then disappears into FDR.

### 3.4 No MGF splitter — **must add**

Add `msgf split --spectra in.mgf --out-dir d/ --chunk-size <N|--chunks N>`, built on the existing
streaming `MgfReader` so it never loads the file. Emits `chunk_%04d.mgf` plus a manifest of
`(file, first_index, n_spectra)`. Doing this in Rust rather than Python/awk keeps the container a
single static binary.

### 3.5 Decoys must be generated once — **pipeline discipline, no code change**

`--tda` per chunk would rebuild identical reversed decoys N times. Reversal is deterministic so the
result would agree, but it is wasted work and an unnecessary divergence risk. The pipeline runs
`msgf decoy` once in PREPARE_DB and passes a concatenated FASTA to every SEARCH task, which then
runs **without** `--tda` (decoys are detected by prefix and never regenerated).

### 3.6 `--db-size` — pin it, don't trust it by default

The default is the per-run candidate-index size. Given the same FASTA and the same digest/mod
parameters, every chunk builds an identical index and gets an identical number, so the default is
*currently* safe. It is safe by coincidence of determinism, not by design. The pipeline computes the
size once in PREPARE_DB and passes `--db-size` explicitly to every SEARCH task, so E-values cannot
drift if index construction ever becomes non-deterministic (e.g. a parallel dedup with unstable
ordering). Needs a small `msgf index --fasta … --print-size`-style path, or PREPARE_DB reads it from
the `index: N peptides -> M candidates` line `msgf search` already prints to stderr.

---

## 4. Pipeline structure

`nextflow/main.nf`, DSL2, one process per stage.

| Process | Runs | Command | cpus / mem |
|---|---|---|---|
| `PREPARE_DB` | once | `msgf decoy` (if `--tda`) + report index size | 1–4 / index-sized |
| `SPLIT` | once per input MGF | `msgf split` | 1 / small |
| `SEARCH` | once per chunk (parallel) | `msgf search --threads ${task.cpus} --db-size … --spec-index-offset … --spec-file-name …` | `params.cpus` / index-sized |
| `FDR` | once | `msgf fdr -i chunk*.tsv -o final.tsv` | 1 / result-sized |
| `REPORT` | once | summary: PSMs, targets, IDs at 1% and 5% FDR | 1 / small |

Multiple input MGFs fan out naturally: `--spectra 'data/*.mgf'` splits each file, searches all
chunks from all files, and merges everything into **one** FDR calculation — which is the statistically
correct treatment of one experiment's runs, and matches what a single `msgf search` over a
concatenated MGF would do.

`--threads ${task.cpus}` is not optional. `msgf search` defaults to *all cores on the machine*
(`rayon` global pool), which on a shared SLURM node means oversubscribing every co-tenant job.

---

## 5. The repeated index build — the real efficiency limit

Every SEARCH task rebuilds the peptide index from the FASTA (`msgf-cli/src/search.rs:305` →
`msgf-search/src/index.rs:46`).
Nothing serializes it. With *N* chunks the pipeline pays *N* × index build, and that fixed cost per
task sets a **floor on useful chunk size**: past a point, adding chunks adds index builds faster than
it removes search work.

The sizing rule, to be stated in the README with measured numbers from §7's G4:

> Let `T_index` be the index build time and `T_spec` the mean per-spectrum search time. A chunk
> should hold **at least `10 × T_index / T_spec` spectra** so the index is ≤10% overhead.

Mitigations, in the order they should be considered:

1. **Bigger chunks** — free, always try first. `--chunk-size` default should follow the rule above.
2. **Index caching** (future): `msgf index build` → an mmap-able file, `msgf search --index`. A real
   feature with its own validation burden (an index file must produce results identical to an
   in-process build) and its own plan. Do not start it here.
3. **Never** split the database to avoid the cost (§2).

**Memory is the harder resource.** Peak RSS per SEARCH task is dominated by the peptide index, which
is a function of the FASTA, enzyme, missed cleavages and variable mods — **not** of chunk size.
Smaller chunks therefore do *not* reduce memory; they multiply it, because *N* tasks on one node each
hold a full index. This must be documented, because it is the opposite of the usual intuition and it
is how people OOM a cluster. `nextflow.config` gets a retry-with-escalating-memory directive, since
index size is not knowable a priori.

---

## 6. Container and profiles

**The image is trivially small and MIT-clean.** `release.yml` already builds an
`x86_64-unknown-linux-musl` static binary, and the default scoring model is *bundled inside it*
(`msgf-scorer/models/`, trained from CC0 MassIVE-KB). So the image is a base plus one file, with no
model download, no JVM, and **no UC-licensed bytes** — which matters, because a public container
image is a redistribution.

One gotcha worth stating so nobody wastes an afternoon: **`FROM scratch` will not work.** Nextflow
executes a generated `.command.run`/`.command.sh` inside the container and needs `/bin/sh` plus
basic coreutils. Use `alpine` (musl-native, ~8 MB) — not `distroless`, which also lacks a shell.

Extend `release.yml` to build and push the image on a `v*` tag, so the image version and the binary
version cannot diverge.

Profiles in `nextflow.config`:

| Profile | Target | Notes |
|---|---|---|
| `standard` | local | default; `params.cpus` bounded by the machine |
| `slurm` | HPCC | `process.executor = 'slurm'`, queue/partition params |
| `k8s` | NRP/Nautilus | storage class + PVC for inputs; see the `nrp-nautilus` conventions |
| `docker` / `singularity` | container engines | composable with any executor profile |

No `conda` profile — the artifact is a static binary; a Conda recipe would add a packaging surface
for no gain.

---

## 7. Validation gates

The repo's culture is that correctness claims are tested, not asserted. These gates need no Java and
no UC-licensed data — any MGF plus any FASTA exercises them.

| Gate | Check |
|---|---|
| **G1 identity** | `N=8` chunked run vs single `msgf search`: identical PSM rows (after sorting) and identical q-values |
| **G2 chunk invariance** | `N ∈ {1, 2, 8, 64}` all produce the identical result set |
| **G3 FDR** | IDs at 1% FDR from the pipeline == from the single run, exactly |
| **G4 scaling** | wall time vs N; publish the curve and the measured `T_index` floor from §5 |
| **G5 provenance** | every merged row has a unique `SpecID` and the original `#SpecFile` (§3.1, §3.2) |

G1 is the plan's contract. It is also a genuine test of §3.1–§3.6: if `SpecID` collides or `--db-size`
drifts, G1 fails. Ship `nextflow/test/identity_check.sh` that runs both paths and diffs them, and
run it in the same local-only way the rest of the suite runs (there is no test CI workflow —
`.github/workflows/` holds only `release.yml`).

A subtlety for G1/G2: row *order* may legitimately differ, since the merge concatenates per-chunk
outputs. Sort both sides on `(SpecFile, SpecID, SpecEValue, Peptide)` before diffing. Do not "fix"
an ordering difference by changing the search's sort.

---

## 8. Milestones

| # | Deliverable | Gate |
|---|---|---|
| N1 | `msgf split` + manifest | round-trip: split then concatenate == original spectrum set |
| N2 | `--spec-index-offset`, `--spec-file-name`, repeatable `-i` on `msgf fdr` | unit tests; G5 |
| N3 | `Dockerfile` (alpine + musl binary) + image build on tag | image runs `msgf --version` under Docker *and* Singularity |
| N4 | `main.nf` + `nextflow.config` with the `standard` profile | **G1**, G2, G3 locally |
| N5 | `slurm` + `k8s` profiles | a real multi-node run completes; G1 holds there too |
| N6 | Scaling measurement + README section | G4 curve published, chunk-size rule stated with numbers |

N1–N2 are CLI work in this repo and gate everything else; N4's G1 is the milestone that makes the
pipeline trustworthy.

---

## 9. Risks

**Silent FDR corruption is the severe failure mode.** Every other risk here costs time; getting §2
wrong costs correctness in a way that looks like a successful run — plausible peptides, plausible
q-values, wrong statistics. G1 exists specifically to make that failure loud, so G1 must be run on
every change to the pipeline, not just once.

**Per-task index memory** (§5) is the likely cause of a first failed cluster run. Mitigate with
retry-on-OOM and a documented sizing rule rather than an optimistic default.

**Many small chunks look attractive and are usually wrong** — they multiply both index CPU and index
RAM. The default `--chunk-size` should encode §5's rule so the naive invocation is already sensible.

**Version skew** between the container image and the pipeline. Pin the image tag in
`nextflow.config` to the release it was built from; never `:latest`.

---

## 10. Open questions

1. **Chunk by spectrum count or by peak volume?** Count is simple; MGF spectra vary enough in peak
   count that count-based chunks may straggle. Measure in G4 before adding complexity.
2. **Should `SPLIT` be avoidable?** For many input files with one chunk each, splitting is pure
   overhead — the pipeline could pass files straight through. Worth a fast path once G1 holds.
3. **Cloud object storage.** Nextflow gives S3/GS staging essentially for free, and the pipeline
   would inherit it; whether to test and document that in v1 or defer is unresolved.
4. **Index caching** (§5 mitigation 2) is the largest available win beyond linear scaling, and is
   plausibly its own plan rather than a milestone here.
