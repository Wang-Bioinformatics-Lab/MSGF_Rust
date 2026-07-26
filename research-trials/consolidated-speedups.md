# Consolidating the measured speedups — trial report

Branch `perf/consolidated-speedups`. Date 2026-07-25.

This is the report for the branch that takes every previously-measured optimization off its
experimental branch, puts them together, gives two of them the caller they never had, and adds three
new ones. It is the first time these results have been measured *in combination* rather than each
against its own baseline.

**Headline.** On the 1,406-spectrum F13 set, single thread, nominal grid, `HCD_HighRes_Tryp.param`,
`-ti 0,1`:

| | ms / 1,406 | spectra/s | vs `main` | CPU-hr / 100k | vs Java |
|---|--:|--:|--:|--:|--:|
| MS-GF+ Java (`TimeGenFunc`, JIT-warm) | 18,979 | 74 | — | 0.3750 | 1.0× |
| `main` @ `2341234` | 4,482 | 314 | 1.00× | 0.0886 | 4.2× |
| this branch, full distribution | 3,227 | 436 | 1.39× | 0.0638 | 5.9× |
| **this branch, tail-pruned (what `search` drives)** | **2,456** | **572** | **1.82×** | **0.0485** | **7.7×** |

Every number is the mean of three runs on an otherwise-idle machine; run-to-run spread was ≤1.5%.
The Java figure was re-measured on the same machine the same day (passes 1–4 of five, discarding the
cold pass: 19,056 / 18,818 / 19,027 / 19,015 ms).

## Method

```bash
cd rust
cargo run -p msgf-genfunc --example profile --release            # full distribution
cargo run -p msgf-genfunc --example profile --release -- thresh  # tail-pruned, search-like
```

`thresh` mode is new on this branch. It drives the DP the way `msgf search` now does: the tail
threshold for each scan is MS-GF+'s own observed top-hit RawScore, read from
`validation/golden/iprg2013_F13.tsv` (1,253 of 1,406 spectra have one; the rest run unpruned).
Without it the harness measures only the unpruned path, which no longer reflects the product.

Per-stage, mean of three:

| Stage | `main` | branch (full) | branch (`thresh`) |
|---|--:|--:|--:|
| preprocess + scored spectrum | 39 ms | 38 ms | 38 ms |
| per-spectrum node tables | 844 ms | 111 ms | 111 ms |
| graph build | 558 ms | 240 ms | 237 ms |
| generating-function DP | 3,041 ms | 2,818 ms | 2,052 ms |
| **total** | **4,482 ms** | **3,227 ms** | **2,456 ms** |

Allocation over the whole run fell from **750 MB to 72.5 MB**; the graph-build stage alone went from
682.5 MB to 3.4 MB.

## What landed, and what each part was worth

### 1. Consolidation (no new code)

- **Node-score cache + ion-major table sweep** (`spec-tables`), from `worktree-spec-tables-perf`:
  844 → 111 ms, a 7.6× on that stage. Detail in
  [spec-tables-node-scoring.md](spec-tables-node-scoring.md).
- **Exact tail-pruned DP** (`compute_tail_into`, `Prune`, `tilt.rs`, `examples/prunelab.rs`).
  This one was recovered from commit `b31c6d1`, whose branch had been **deleted** — the commit
  survived only as a dangling object. Detail in [dp-pruning-limits.md](dp-pruning-limits.md).

PR #9's sink-ancestor prune is **not** separately merged: it is subsumed by the tail prune, whose
`max_remaining[i] == i32::MIN` test performs the identical elimination.

### 2. Giving the tail prune a caller — the largest single win

The tail prune had been measured at 1.31× on the DP and then sat unused, because **nothing called
it**: `search_at_charge` built the generating function *before* scoring candidates, so it had no
threshold to offer.

Nothing in the generating function depends on the candidates, so the two halves can simply be
swapped. Scoring candidates first yields the RawScore of the worst PSM that will be reported, and
that is exactly the tail threshold. A spectrum with no candidates in its precursor window now skips
the DP entirely.

| | cells/graph | DP | TOTAL | throughput |
|---|--:|--:|--:|--:|
| full distribution | 136,351 | 2,818 ms | 3,227 ms | 436 spectra/s |
| tail-pruned | 67,739 | 2,052 ms | 2,456 ms | 572 spectra/s |
| | **2.01× fewer** | **1.37×** | **1.31×** | |

**F13 understates this by a wide margin.** It identifies essentially nothing, so its
DeNovoScore − RawScore gap has median 69 — the far right of the pruning curve
([measurement-traps.md](measurement-traps.md) §1). The same code measures 2.68× on the DP at
DeNovoScore − 20, where a corpus with real identifications sits.

### 3. Tail-pruned `rescore`

`msgf rescore` caches one generating function per `(scan, charge)` and reuses it across every PSM
sharing that key, so no single PSM's RawScore is a valid threshold — the **group minimum** is.
Restructured from one-pass-per-PSM into group-by-key, two-passes-per-group; rows and skip
diagnostics are buffered and replayed in input order.

Generating-function stage 3.42 → 2.71 s (1.26×), whole command 1.27× median over 20 interleaved
pairs, peak RSS 26.2 → 13.0 MB, output sha256-identical.

### 4. Fusing the DP's two per-edge passes

[dp-pruning-limits.md](dp-pruning-limits.md) §4 established a per-edge cost floor of ~2.9 ns/edge
(~16% of the unpruned DP) that caps *every* cell-pruning idea near 6×. The DP visited each incoming
edge twice — once for the destination's score range, once to convolve — gathering the same
`NodeDist` and re-deriving the same `node_score + edge_score` both times. The range pass now caches
`(src_start, src_len, src_min, score_diff, prob)` in a 32-entry stack array that the convolution
consumes.

| path | ratio | evidence |
|---|--:|---|
| unpruned | 1.037 | mean of 15 interleaved pairs, 14/15 favour, spread 0.997–1.078 |
| pruned | ~1.06 | 4 interleaved pairs, 1.031–1.080 |

The candidate flagged the pruned path as an unmeasured regression risk — a node whose row is later
emptied by the raised floor pays for descriptor writes it never uses. Measured, it is the opposite:
the fusion helps *more* on the pruned path, which is the one `search` drives.

Cell counts are unchanged (136,351/graph before and after) — this removes loads, not work. Nodes
with in-degree > 32 fall back to the original two-pass code; measured maximum in-degree on F13 is 21.

### 5. Graph build

Three changes, and an ablation that shows which one mattered:

1. `edge_prob: Vec<f64>` → `edge_aa: Vec<u16>` + an `aa_prob` table. The per-edge stream is 2 bytes
   instead of 8 and the DP reads the identical `f64` bit pattern. The index is explicit, never
   derived from nominal mass, because distinct amino acids share nominal masses (I/L at 113, K/Q at
   128) and two parallel edges can differ only by which one they are.
2. `build_reverse_graph_into` writes into a caller-owned `Graph`, so a whole run reuses one set of
   CSR buffers instead of allocating ~0.5 MB per spectrum.
3. The two-pass edge count collapsed: a node's incoming-edge count is `#{aa : nominal ≤ m}`, which
   saturates at `|aa|` above the heaviest residue, so only ~186 of ~1,528 nodes need a real scan.
   Edge-score constants hoisted out of the fill loop.

| | graph-build | ratio |
|---|--:|--:|
| baseline | 540 ms | 1.00× |
| + memory/allocation work only (1 and 2) | 476 ms | 1.14× |
| + edge-score constant hoist (3) | 238 ms | **2.27×** |

**The memory work was worth 1.14×; the constant hoist was worth 2.0× on top of it.** That ordering
was not what the candidate was set up to find, and it is the reusable lesson: the 682 MB of
allocation traffic was the conspicuous number, but arithmetic hoisted out of a hot loop was the
actual cost. This is [measurement-traps.md](measurement-traps.md) §3 again — allocation counters
point at anomalies, not bottlenecks.

The DP did **not** measurably benefit from the 4×-smaller probability stream (1.002–1.014, noise).

## Defects caught by adversarial review

Each candidate was reviewed by two independent agents, one hunting fidelity breaks and one hunting
correctness/edge-case bugs. Neither refuted the arithmetic of any candidate, but three real defects
surfaced — two of which would otherwise have shipped:

1. **Degenerate alphabet corrupted every CSR offset.** The fused edge count saturates above the
   heaviest residue, but an alphabet whose nominal masses are all ≤ 0 drives `saturated_from` to 0,
   so the saturating loop credited node 0 — the source, which has no incoming edges — with a full
   edge list and shifted every subsequent offset. Fixed by clamping to 1; the counting block is now
   `fill_edge_offsets`, unit-tested against a naive per-node count over 7 alphabet shapes × 5 graph
   sizes. The test was confirmed to fail without the clamp.
2. **A 256-entry alphabet cap could abort a legitimate search.** `msgf-search` builds one alphabet
   entry per (variable mod × target residue), and an unrestricted variable mod targets all 20 base
   residues, so ~12 such mods reach 260 entries and would have panicked inside a rayon worker. The
   index is now `u16`. Since the memory half of that optimization was worth only 1.14×, the second
   byte costs almost nothing.
3. **`msgf rescore --out` was opened after scoring.** An unwritable path failed after a full
   multi-minute run instead of immediately.

Separately, `GenFunc::spectral_probability` now asserts `raw_score >= valid_from` in **release**,
not only debug. Below the pruning threshold those cells were never computed, so the old behaviour
was a silently-too-small SpecEValue — far worse than a panic, and the check is one predictable
branch per PSM against a whole DP.

One reviewer returned `refuted: true` on the graph-build candidate. Its blocker was that the patch
did not apply to the branch and left three dangling `edge_prob` reads — exactly the three sites
already fixed by hand during the 3-way merge, so it was independent confirmation of the merge rather
than a new defect.

## Fidelity

Nothing here is an approximation. The tail prune is exact **at and above its declared threshold**;
everything else is exact everywhere.

| Gate | Result |
|---|---|
| `msgf search` on F13 vs pre-change `main` | **byte-identical**, same sha256, 1,256 PSMs |
| F13 end-to-end oracle (`golden_search`, `--ignored`) | 1161/1161 exact RawScore and DeNovoScore, SpecEValue in tolerance |
| `golden_specprob` | DeNovoScore 30/30, SpecEValue 30/30 |
| `golden_rescore` | raw 30/30, denovo 30/30, spec 30/30 |
| `cargo test --workspace --release` | green |
| `cargo clippy --workspace --all-targets` | 0 warnings |
| `cargo fmt --all --check` | clean |

New `f64::to_bits` tests, per [measurement-traps.md](measurement-traps.md) §7 (an exact change must
pass a bitwise test, not an epsilon):

- `fused_edge_pass_is_bit_exact` — fused vs. pre-fusion DP over 25 seeds × narrow/wide in-degree
  (straddling the 32 bound) × cleavage on/off × unpruned and four prune thresholds.
- `compact_graph_is_bit_exact_against_the_previous_builder` — the new builder vs. a verbatim copy of
  the old one on 40 real F13 spectra, with deliberately non-uniform amino-acid probabilities so an
  `edge_aa` mix-up would show, plus `reuse_survives_shrinking_and_growing` for buffer reuse across
  sizes 2500 → 400 → 1200 → 2500 → 300.
- `rescore::tests::pruned_matches_unpruned_bitwise` — pruned vs. unpruned SpecEValue over 250 golden
  PSMs plus injected bad rows.
- `fused_offsets_match_the_naive_count` — the CSR offset regression above.

## Where the remaining opportunity is

The stage split has changed shape, and that changes what is worth doing next:

| Stage | share on `main` | share now (`thresh`) |
|---|--:|--:|
| DP | 67.9% | **83.5%** |
| graph build | 12.5% | 9.6% |
| node tables | 18.8% | 4.5% |

Non-DP work can no longer fund much: even reducing graph build and node tables to zero is worth only
1.17× from here. And cell pruning is capped near 6× by the per-edge floor. So the next real step is
a gate that **does not run the DP at all** for spectra that cannot be significant — the two-sided
corridor and exact early-exit rejection in `plans/PLAN3.md` §5.2 — rather than a cheaper cell.

Not attempted here, and still open: sharing the reverse-reachability/`max_remaining` sweeps across
the two isotope candidates; sparse/dense score distributions; the saddlepoint estimator, which
remains on `worktree-genfunc-algo-speedups` and is deliberately not wired into anything.

## Reproduction

Requires the gitignored `validation/data/` and the MS-GF+-derived F13 golden.

```bash
cd rust
cargo run -p msgf-genfunc --example profile   --release            # full distribution
cargo run -p msgf-genfunc --example profile   --release -- thresh  # search-like path
cargo run -p msgf-genfunc --example prunelab  --release            # pruning curve vs threshold
cargo test --workspace --release
cargo test -p msgf-genfunc --release --test golden_specprob -- --nocapture
cargo test -p msgf-search  --release --test golden_search -- --ignored --nocapture

# byte-identical search output vs a pre-change build
./target/release/msgf search -s ../validation/data/spectra/F13.mgf \
  -d ../validation/data/fasta/iprg2013_human.revCat.fasta \
  --mods ../validation/data/config/iprg-2013_Mods.txt \
  -p ../validation/data/models/HCD_HighRes_Tryp.param \
  -t 10ppm --ti 0,1 -e 1 --threads 6 -o /tmp/f13.tsv
```
