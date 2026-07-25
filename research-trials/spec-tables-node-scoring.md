# Trial: per-spectrum node tables (`spec-tables`)

Raw measurement record for the `spec-tables` work summarized in
[generating-function-optimization.md](generating-function-optimization.md) §*Non-DP Stages*. This
file keeps the instrument output and the reasoning steps, including the wrong hypotheses, so the
investigation does not have to be redone to be trusted.

**Branch:** [`worktree-spec-tables-perf`](https://github.com/Wang-Bioinformatics-Lab/MSGF_Rust/tree/worktree-spec-tables-perf)
(`7047f3a` node-score cache, `c697864` ion-major sweep). Branched from `main` at `1b14c5f`, so it
does **not** contain draft PR #9 — its baseline is 314 spectra/s.

**Workload for every number below:** 1,406 F13 spectra, single thread, nominal grid,
`HCD_HighRes_Tryp.param`, `-ti 0,1`.

## Model shape

Measured, not assumed — the optimization only works because these numbers are small:

```
num_segments = 2
partitions   = 92
max_rank     = 150            (cache stride = 151 bins per ion)
fragment ions per partition: min 3, median 4, max 6
```

## Step 1 — inner-loop census

`node_score` forms a theoretical m/z for every (segment, ion) pair and then discards the iteration if
the ion has the wrong polarity or if `segment_num(theo) != seg`. How much survives:

```
node_score calls (nodes x 2 polarities): 4,295,210
inner iterations:            34,949,042   (8.1 per call)
  discarded, wrong polarity:   17,474,521   (50.0%)
  discarded, wrong segment:     8,687,844   (24.9%)
  reached the peak lookup:      8,786,677   (25.1%)
```

So a 4.0× ceiling on skipping dead iterations — but that says nothing about *time*, because the
surviving 25% carry all the peak lookups.

## Step 2 — locating the cost (the step that mattered)

The census invites the conclusion "remove the dead iterations." That would have been wrong as a first
move. Timing a hand-written replica of the same triple loop, with and without the peak lookup,
against the real `tables()`:

| | time |
|---|--:|
| real `tables()` | 791 ms |
| replica loop, no peak lookup | 113 ms |
| replica loop, with peak lookup | 234 ms |

The replica does the *same* iterations and the *same* lookups in 234 ms. **~70% of the real cost was
neither the loop nor the lookup**, and no amount of iteration-skipping would have touched it.

The gap is `ScoringModel::score_from_table`, called once per surviving ion per node — 8.8M times.
Each call:

1. linear-scans all 92 partitions' rank distributions to find this partition,
2. linear-scans that distribution's rows comparing the ion's **name `String`**,
3. takes an `ln`.

None of it depends on the node being scored. This is the same defect already fixed for *edge*
scoring (`ion_existence_cache`, `error_score_cache`); the node path never got the treatment.

## Step 3 — the node-score cache, and a wrong hypothesis

`NodeScoreCache` precomputes every distinct result per spectrum, indexed by
`(segment, ion, rank bin)`. First result:

| Stage | before | after |
|---|--:|--:|
| `spec-tables` | 820 ms | 341 ms |
| `scored-spec` (builds the cache) | 12 ms | **120 ms** |

The cache build cost 108 ms — a third of the win. The allocation counters showed reallocs going
0.0 → 8.9 per spectrum and allocation traffic 9.8 → 32 MB, which pointed at `Vec` growth.

**That hypothesis was wrong.** Adding `Vec::with_capacity` removed the reallocs (8.9 → 0.0, 32 →
16.6 MB) and moved the time by ~5 ms. The real cost was that the per-bin helper called
`score_from_table` once per bin, so the 92-partition scan and the name-string row lookup were being
re-paid on each of the **151 rank bins**. Resolving the rank distribution and the ion's row once per
ion (`extend_score_bins`) took the build from 115 ms to **36 ms**.

Lesson: allocation counters pointed at a real anomaly that was not the bottleneck. The same thing
happened in the original DP optimization (`PERFORMANCE.md` records that allocation "was never the
bottleneck"). Measure the hypothesis, not the counter.

## Step 4 — ion-major sweep

Only then was claiming the dead iterations worthwhile. Both filters are properties of
`(segment, ion)` and a node *range*, not of a node: `FragOff::mz` is monotone in node mass and
`segment_num` is monotone in m/z, so `segment_num(theo(k))` is a non-decreasing step function of the
node index and the surviving nodes are contiguous.

`tables()` now loops `(segment, ion)` outer, resolves the range once by binary search, and sweeps it.

Two design points worth preserving:

- **`segment_node_range` binary-searches the predicate itself**, evaluated with the same float
  arithmetic as the per-node path — deliberately *not* an algebraic inversion of `segment_num`. An
  inverted formula is a second expression that can round differently, and this code is held
  bit-exact.
- **The peak lookup became monotone for free.** Within one ion the match window only moves right, so
  a single cursor walks the peak list instead of re-entering through the bucket index and
  backtracking on every node. This was a consequence of the restructure, not a separate change.

## Result

| Stage | main | + node cache | + ion-major sweep |
|---|--:|--:|--:|
| preprocess | 24 ms | 24 ms | 23 ms |
| scored-spec | 12 ms | 36 ms | 35 ms |
| **spec-tables** | **820 ms** | **341 ms** | **111 ms** |
| graph-build | 533 ms | 546 ms | 540 ms |
| DP compute | 3085 ms | 3166 ms | 3082 ms |
| **pipeline** | **4475 ms** | 4114 ms | **3792 ms** |
| throughput | 314/s | 342/s | **371/s** |

`spec-tables` 18.3% → **2.9%** of the pipeline; 7.4× on the stage, 1.18× overall. The profile shape
changed as a result: DP is now 81.3% and graph build 14.2%, so DP work is worth proportionally more
than the original Amdahl estimate suggested.

## Fidelity

Bit-exact throughout — the cached value comes from the identical expression, and `prefix[k]` still
accumulates prefix ions in `(segment, ion)` order, which is exactly the order `node_score` adds them.

- `node_score_cache_is_bit_exact` — cached vs uncached, `to_bits`, every node, both polarities.
- `tables_match_node_score` — ion-major `tables()` vs per-node `node_score`, `to_bits`, every node,
  both polarities, across four charge/parent-mass combinations (the derived ranges move with parent
  mass, so one mass would not exercise them).

Both use the committed bundled model, so they run on a clean checkout with no fetched data.

A model whose rank distributions do not cover every scored ion declines the cache and takes the
uncached node-at-a-time path, preserving the original panic site rather than papering over it.

Full gate: `cargo test --workspace`, the `#[ignore]`d F13 search golden (1161/1161 exact RawScore and
DeNovoScore), `cargo clippy --workspace --all-targets`, and `msgf search` on F13 producing a
**byte-identical** output TSV.

## Not done

- **Cache per model rather than per spectrum** (~35 ms/run). Needs a field on `ScoringModel`, which
  breaks its `PartialEq`/`Clone` derives, or threading a table through `from_ranked_peaks` and its
  six call sites. Not obviously worth it at this size.
- **`node_masses()` is still a separate pass** inside `tables()`. Fusing it saves a pass over the
  node array, but its peak lookups remain, so the ceiling is small.

## Reproducing

```bash
git checkout worktree-spec-tables-perf
cd rust
cargo run -p msgf-genfunc --example profile --release      # the stage table above
cargo test -p msgf-scorer --release --lib                  # both bit-exactness tests
cargo test -p msgf-search --release --test golden_search -- --ignored
```

The census and replica-timing numbers came from a throwaway `msgf-scorer` example (`tableshape`)
that was deleted before commit along with the `debug_*` accessors it needed; the figures above are
its recorded output. Recreating it means re-adding accessors for `seg_partition`, `segment_num` and
the peak lookup.
