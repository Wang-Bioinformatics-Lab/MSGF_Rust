# Generating-Function and Scoring-Pipeline Optimization Trials

This report records measured experiments for reducing the generating-function dynamic program and
its surrounding scoring pipeline. It is evidence for a future consolidated implementation, not a
merge plan for any one experimental branch.

Numerical compatibility with MS-GF+ remains the primary constraint: exact optimizations must preserve
edge order, floating-point operations, the complete score distribution, DeNovoScore, and SpecEValue
unless an explicitly narrower API is introduced.

One explicitly *approximate* route is also recorded, in its own section, because it is the largest
measured reduction available. It is fenced off from the bit-exact path and nothing in the search or
CLI calls it.

## Experimental References

None of the references below are intended for direct merging; all are records of measured experiments.

### Draft PR #9 — prune dead DP subgraphs

Draft PR [#9 — prune dead DP subgraphs](https://github.com/Wang-Bioinformatics-Lab/MSGF_Rust/pull/9)
is a reference implementation. It explores two complementary optimizations:

1. **Remove redundant range work.** Each node's destination range is already constructed as the
   union of its shifted predecessor ranges. The convolution can therefore use the full predecessor
   slice without recomputing overlap clipping or bounds checks.
2. **Prune nodes outside sink paths.** A reverse CSR walk marks the union of all requested sinks'
   ancestors. The forward probability DP skips source-reachable nodes that cannot contribute to a
   sink.

On the 1,406-spectrum F13 benchmark, draft PR #9 produced these measured gains:

| Measurement | Before | PR #9 | Improvement |
|---|---:|---:|---:|
| DP compute | 3.09 s | 2.50 s | **19.1% faster** |
| Full pipeline | 4.48 s | 3.90 s | **12.9% faster** |
| Throughput | 314 spectra/s | 360 spectra/s | **14.6% higher** |
| Distribution cells per graph | 136,351 | 107,746 | **21.0% fewer** |

The sink-ancestor pruning step alone was approximately 15.6% faster in the DP and 10.6% faster
end-to-end. Golden DeNovoScore and SpecEValue remained exact for all 30 checked PSMs.

### Branch `worktree-genfunc-algo-speedups`

Branch [`worktree-genfunc-algo-speedups`](https://github.com/Wang-Bioinformatics-Lab/MSGF_Rust/tree/worktree-genfunc-algo-speedups)
(commit `7c8a314`) implements threshold-aware score pruning and a saddlepoint tail estimator, both
described below. It branches from `main` at `1b14c5f` and therefore does **not** contain PR #9, so
its measurements are relative to the unoptimized DP (the same 314 spectra/s baseline).

The two experiments overlap: the branch's `max_remaining[i] == i32::MIN` test is exactly PR #9's
sink-ancestor prune, arrived at as a special case of a more general bound, so **their gains are not
additive**. PR #9's first optimization also interacts — the pruned path raises per-node score floors,
so predecessor slices no longer fit whole and overlap clipping is needed again. The unpruned path
keeps it.

### Branch `worktree-genfunc-aggressive-prune`

Branch `worktree-genfunc-aggressive-prune` (commits `b31c6d1` and `76ad840`) ports exact threshold
pruning onto the newer mainline DP and measures its practical ceiling. It also retains reproducible
implementations of two rejected lossy policies and records certified error bounds for them. This
branch is experimental and was not pushed when this report was assembled.

### Branch `worktree-spec-tables-perf`

Branch [`worktree-spec-tables-perf`](https://github.com/Wang-Bioinformatics-Lab/MSGF_Rust/tree/worktree-spec-tables-perf)
(commits `7047f3a` and `c697864`) attacks the largest non-DP stage with a node-score cache followed by
an ion-major table sweep. Both steps are bit-exact and are described below.

## Where the DP's Time Goes

Single-thread, 1,406 F13 spectra, nominal grid, `HCD_HighRes_Tryp.param`, `-ti 0,1`, measured on
`main` at `1b14c5f` with `cargo run -p msgf-genfunc --example profile --release`:

| Stage | Share |
|---|--:|
| Preprocess and scored spectrum | ~0.8% |
| Per-spectrum node tables (`spec-tables`) | 18.3% |
| Graph build | 11.9% |
| Generating-function DP | 68.9% |

Per graph: 1,528 nodes (1,349 reachable), 29,576 edges, 136,351 distribution cells, and an average
score-support width of **101**. That width is the opportunity: every edge relaxation costs about 101
floating-point operations to convolve a whole distribution, when the quantity actually consumed is a
single tail probability.

This also bounds the exercise. With the DP at 68.9%, an infinitely fast DP still caps the
per-spectrum pipeline at approximately 3.2×; reductions past roughly 3× are not realizable until
`spec-tables` and graph build move as well. Separately, end-to-end `msgf search` on F13 is
index-build dominated (approximately 5 s of 6 s), so none of this is visible there — it applies to
the per-spectrum stage, which is what scales with spectrum count.

**Updated after the complete `spec-tables` trial below.** With that stage cut from 820 ms to 111 ms,
the breakdown is now `spec-tables` 2.9%, graph build 14.2%, and DP **81.3%**. Throughput moved from
314 to 371 spectra/s, making DP edge cost and graph construction the dominant remaining targets.

## Non-DP Stages

The DP is the interesting algorithm, but it is not the only place per-spectrum time goes, and the
other stages turned out to hide larger constant-factor mistakes than the DP does.

### Per-spectrum node tables — fixed

`spec-tables` was 18.3% of the pipeline, and almost none of it was the work it appears to do. An
inner-loop census over the F13 set found that only **25%** of `node_score`'s iterations reach a peak
lookup — 50% are discarded on a polarity mismatch and 25% because the theoretical m/z does not land
back in the segment being iterated. But timing a replica of that entire loop *including* the peak
lookups gave 234 ms against the real 791 ms, so roughly 70% of the cost was neither the loop nor the
lookup.

It was `ScoringModel::score_from_table`, called once per surviving ion per node — 8.8M times over the
F13 set. Each call linear-scans every partition's rank distribution to find this partition (92 of
them in `HCD_HighRes_Tryp`), then linear-scans that distribution's rows comparing the ion's *name
string*, then takes an `ln`. None of it depends on the node being scored.

The first fix on branch
[`worktree-spec-tables-perf`](https://github.com/Wang-Bioinformatics-Lab/MSGF_Rust/tree/worktree-spec-tables-perf)
by precomputing every distinct result per spectrum into a `NodeScoreCache` indexed by
`(segment, ion, rank bin)` — the same treatment `ion_existence_cache` and `error_score_cache` already
give edge scoring. A second step mattered as much as the first: filling a whole ion's row with the
rank distribution and the row resolved **once** (`extend_score_bins`), rather than re-paying the
92-partition scan on each of the 151 rank bins, took the cache-build cost from 115 ms to 36 ms.

The second fix claims the dead iterations. For a fixed `(segment, ion)`, theoretical m/z is monotone
in node mass and the surviving nodes form a contiguous range. `tables()` can therefore loop over
ions first, locate the valid range once by binary search, and sweep only those nodes. Within an ion
the peak window also moves monotonically, so one cursor replaces repeated bucket lookups and
backtracking.

| Stage | Main | + node cache | + ion-major sweep |
|---|---:|---:|---:|
| `spec-tables` | 820 ms | 341 ms | **111 ms** |
| Pipeline | 4,475 ms | 4,114 ms | **3,792 ms** |
| Throughput | 314 spectra/s | 342 spectra/s | **371 spectra/s** |

Bit-exact: the cached value comes from the identical expression and `node_score` adds the same terms
in the same order. `node_score_cache_is_bit_exact` compares cached against uncached with `to_bits`;
`tables_match_node_score` pins the ion-major path over every node, polarity, and several
charge/parent-mass combinations. A model whose rank distributions do not cover every scored ion
declines the cache and takes the uncached path, preserving the original behavior.

The branch reports a clean workspace suite, ignored F13 search golden (1,161/1,161 exact RawScore
and DeNovoScore), Clippy, and a byte-identical F13 search TSV.

### Remaining non-DP leads

- **Build the node-score cache once per model rather than per spectrum** (~24 ms/run). Needs either a
  field on `ScoringModel`, which breaks its derives, or threading a table through
  `from_ranked_peaks`; neither is obviously worth it at this size.
- **Graph build is now the second stage at 14.2%** and still allocates 682 MB per run, much of it the
  `edge_prob` array (see *Beyond the DP*).

## Further Exact-Pruning Ideas

### Threshold-aware score pruning — implemented and measured

Search previously computed the full generating function before candidate RawScores. The specialized
search path is:

1. score candidates and retain the top matches;
2. determine the lowest RawScore whose tail probability is needed;
3. compute an optimistic maximum remaining score from every node to a sink;
4. discard a state only when `current_score + max_remaining < required_score`.

Implemented on the branch as `compute_tail_into`, with `search_at_charge` reordered so candidates are
scored before the generating function is built. Nothing in the generating function depends on the
candidates, so the reordering is free, and spectra with no candidates now skip the DP entirely.

The retained cells are reached by the same multiply-adds in the same order as in the unpruned DP, so
they are bit-identical `f64` rather than approximately equal. `max_remaining` is one reverse integer
sweep of the existing CSR, approximately 1% of DP cost.

The two obstacles previously recorded here resolved as follows:

- **DeNovoScore from a truncated distribution.** No separate calculation is needed.
  `max_remaining[0]` is the source's best full-path score, so the sweep that supplies the bound also
  supplies the exact DeNovoScore. Clamping the effective threshold to it keeps the maximum-score cell
  alive even when the caller requests a threshold no peptide can reach.
- **Arbitrary `spectral_probability` queries.** `GenFunc` gained a `valid_from` field recording the
  threshold; querying below it is debug-asserted. Callers needing the whole distribution keep
  `compute_into`, which is unchanged.

Rescoring remains unaddressed: `msgf rescore` caches one generating function per `(scan, charge)`
across several PSMs, so it needs the per-key minimum RawScore as its threshold. It still uses the
unpruned `compute`.

Measured DP reduction over the 1,253 F13 spectra that have a golden MS-GF+ RawScore (2,506 graphs):

| Threshold | DP cells | DP time |
|---|--:|--:|
| DeNovoScore − 5 | 280× fewer | 6.5× |
| DeNovoScore − 10 | 75× fewer | 4.9× |
| DeNovoScore − 20 | 19× fewer | 3.0× |
| DeNovoScore − 40 | 5.6× fewer | 2.0× |
| DeNovoScore − 80 | 2.3× fewer | 1.5× |
| F13's own MS-GF+ RawScores | 2.1× fewer | 1.4× |

Two cautions on reading this. First, the cell reduction outruns wall-clock because narrower `axpy`
slices vectorize less effectively — a 2.1× cell reduction bought 1.4×. Second, **F13 is the worst
case and should not be used to judge the method**: it identifies essentially nothing, so its
DeNovoScore−RawScore gap has median 69, at the far right of the table. Even its best 50 PSMs by
SpecEValue have a median gap of 25. A corpus with real identifications sits in the 3–6× band. Cost
falls as the match improves, which is where the SpecEValue matters most.

Validation: bit-identical on synthetic graphs across thresholds, seeds, and cleavage settings
(`f64::to_bits`); bit-identical on 1,253 F13 spectra at MS-GF+'s own RawScores (`f64::to_bits`); and
the full F13 search produces a byte-identical output TSV. The `#[ignore]`d `golden_search` test is
unchanged at 1161/1161 exact RawScore and DeNovoScore.

### Reuse reachability across isotope candidates

The `-ti 0,1` graphs share one edge structure and have adjacent sinks. Investigate caching or
incrementally updating the reverse-reachability mask rather than rebuilding it for each candidate.
Measure the mask-building cost first; it is small relative to convolution. This now also covers the
`max_remaining` sweep, which has the same shape and the same per-candidate rebuild.

### Hybrid sparse/dense distributions

Early nodes may contain many zero score cells before distributions become dense. A sparse
representation could avoid those operations and switch permanently to the current contiguous arena
above a measured density threshold. Entries must retain score order, and conversion must not change
addition order.

## How Much Harder Can the DP Prune?

The aggressive-pruning branch measured the exact threshold prune and three attempts to push beyond
it. Reproduce these results with:

```bash
cargo run -p msgf-genfunc --example prunelab --release
```

All numbers use the 1,406-spectrum F13 set, single-threaded on the nominal grid, with
`HCD_HighRes_Tryp.param` and `-ti 0,1`.

### Exact-pruning ceiling

| Threshold | Cells per graph | Reduction | DP speedup |
|---|---:|---:|---:|
| MS-GF+ RawScore | 69,162 | 2.11× | **1.31×** |
| DeNovoScore − 40 | 22,137 | 6.2× | 1.76× |
| DeNovoScore − 20 | 5,973 | 22.8× | 2.68× |
| DeNovoScore − 10 | 1,414 | 96.4× | 4.31× |
| DeNovoScore − 5 | 350 | 390× | **6.17×** |

At DeNovoScore − 5, 390× fewer cells produce only a 6.2× speedup. Approximately 170 µs per graph
remains, or 16% of the unpruned DP, in the `max_remaining` and score-range edge passes. Those passes
are independent of score-support width. On this graph shape, no cell-pruning policy can beat roughly
6× without also making an edge visit cheaper.

### Rejected aggressive prunes

- **Top cap.** Retaining only `N` score units above the threshold is invalid because partial path
  scores are not monotone: signed node/edge scores and the cleavage penalty allow a path to rise far
  above its final score. At cap +30, the mean absolute log10 error was 0.50 and the worst was 27.2,
  despite a 1.55× DP speedup.
- **Chernoff/tilted low-end trim.** Bounding each cell's contribution with a backward tilted sweep
  removed only another 3% of cells at `ε = 1e-3`; spending 100× more error budget bought another
  1.8%. The sweep itself cost about 13% of the unpruned DP, resulting in 0.92×–1.03× overall.
- **Extra dead-node sweep.** A forward `max_achievable` pass can prove that some nodes lie on no path
  clearing the threshold, but measured 7–15% slower. The existing empty-predecessor range pass exits
  more cheaply than the extra arithmetic needed to classify every live edge.

The lossy experimental API removes probability only, so its result is a lower bound.
`GenFunc::err_bound` accumulates a certified upper bound on discarded mass and
`GenFunc::relative_error()` reports the interval width. Any future lossy experiment should expose
its uncertainty as explicitly.

The aggressive branch reports exact golden SpecEValue/rescore/search results, synthetic
`f64::to_bits` tail comparisons, and a bitwise scalar-versus-AVX kernel test.

## Approximate Routes (Outside the Bit-Exact Path)

### Saddlepoint inversion of the cumulant function

The score axis exists only so one number can be read off the end. Evaluating the generating function
at an exponential tilt `θ` collapses it to a scalar per node:

```
Z_θ(m) = Σ_aa e^{θ·(nodeScore(m) + edgeScore)} · p_aa · Z_θ(m − mass(aa))
```

This is the same recursion at width 1 instead of width ~101. Carrying `(Z, dZ/dθ, d²Z/dθ²)` yields
the cumulant function `K = ln Z` and its first two derivatives; Newton-solving `K'(θ̂) = T` and
applying the Lugannani–Rice formula with the lattice continuity correction yields the tail.
Implemented on the branch as `msgf_genfunc::saddle`.

Two properties made it practical:

1. **`Z(0)` is spectrum-independent.** At `θ = 0` every tilt weight is 1 and the recursion degenerates
   to `z[m] = Σ_aa p_aa · z[m − mass(aa)]` — the spectrum drops out entirely, as does the cleavage
   factor. The normalizing constant is therefore a mass-composition table built once per run rather
   than a sweep per graph. Omitting the term instead introduces a systematic +0.03 log10 bias.
2. **The tail is stationary in θ at the saddlepoint**, so solving `K'` to 0.25 score units rather than
   1e-9 costs nothing measurable in accuracy and saves approximately 20% of the sweeps.

Measured at **3.2×** the speed of the exact DP, at approximately 2.5 tilted sweeps per graph.
Accuracy against the exact DP on 1,253 F13 spectra at MS-GF+'s own RawScores: median log10 ratio
−0.002; 96.1% of spectra within the project's 0.05 SpecEValue tolerance, 98.6% within 0.10, and 100%
within 0.30.

The error is a smooth function of tail depth rather than noise. Per decade of *relative* depth
(`tail / Z(0)`, the tail as a fraction of the distribution's own total mass):

| Relative depth | Median \|log10\| | Worst |
|---|--:|--:|
| 1e0 – 1e-6 | 0.002 – 0.004 | 0.021 |
| 1e-6 – 1e-9 | 0.009 | 0.027 |
| 1e-9 – 1e-12 | 0.010 | 0.096 |
| 1e-12 – 1e-16 | 0.040 | 0.323 |
| Beyond 1e-16 | 0.229 | 1.013 |

Accuracy decays as the threshold approaches the maximum achievable score, where too few paths carry
the tail for a normal-based inversion. Real SpecEValues sit near 1e-6 relative depth, well inside the
dependable regime. Note the scale: `Z(0)` is the probability of reaching the precursor mass at all,
so an absolute SpecEValue of 1e-9 is nowhere near 1e-9 relative depth.

### Tiered exact/approximate evaluation

The two implemented reductions are complementary rather than additive, and the complementarity is
the more useful result. Exact pruning is cheap when the PSM is good, which is when the SpecEValue
must be right; it is expensive when the PSM is poor, which is when the SpecEValue is large and will
never approach an FDR cutoff.

A tiered path would run the saddlepoint estimator for every spectrum and escalate to the exact pruned
DP only for PSMs whose approximate E-value falls near the FDR decision boundary. Because the
approximation's measured error is bounded, a PSM more than about 1 log10 from the boundary cannot be
moved across it. Reported values in the region that decides anything stay bit-exact, and perturbing
the ordering far below the cutoff does not change decoy counts above it, so q-values are unaffected.

Not implemented; it requires a policy decision about what "near the boundary" means and whether it
belongs in the search engine or the FDR stage.

## Evaluated and Rejected

- **FFT convolution along the score axis.** Each edge applies a shift, not a general convolution, so
  there is no transform to amortize. The mass axis is not translation-invariant either, because node
  scores differ per node. This is a stronger objection than the fidelity one: the technique does not
  apply.
- **Score-lattice coarsening.** Binning scores by `b` shrinks the DP by `b×`, but SpecEValue moves by
  roughly `e^θ` per score unit, with `θ` typically 0.3–0.5 — approximately 0.15–0.2 log10 per unit of
  threshold misplacement. Already outside the 0.05 tolerance at `b = 2`. Analytical estimate; not
  measured.
- **Contracting score-inert graph regions.** Over a maximal run of nodes carrying zero node and edge
  scores the DP is a pure mass-mixing map with scalar coefficients, so the region could in principle
  be collapsed. The contraction costs approximately 186 entry-node terms per exit node against 21
  taps for the direct DP; the amino-acid kernel is already sparser than its own contraction.
  Analytical estimate.
- **Mass-axis reachability pruning on the nominal grid.** Every node at or above mass 57 is
  reachable, so there is nothing to prune. May be worth revisiting at the low-mass end of the
  high-precision grid, where achievable accurate masses are genuinely sparse. Untested.

## Beyond the DP

Given the Amdahl ceiling above, these now matter as much as further DP work:

- **`spec-tables`** was the largest non-DP stage; it is now 2.9% after the node-score cache and
  ion-major sweep.
- **`edge_prob` is effectively constant** (one value per amino acid) yet is stored as an `f64` per
  edge — approximately 236 KB per graph of redundancy, written during graph build and read on every
  DP edge. Replacing it with a per-amino-acid lookup removes a memory stream from both stages.
- **Per-edge fixed work now limits exact pruning.** Consider fusing reverse-bound and range metadata,
  reducing CSR memory streams, and retaining direct SIMD dispatch. Cell-only pruning cannot cross
  the measured ~6× ceiling.
- **Fewer saddlepoint sweeps.** Approximately 2.5 per graph is near the Newton floor. Fusing three or
  four θ values into one AVX sweep and recovering `K'` and `K''` by finite differences could approach
  one sweep-equivalent, roughly doubling the 3.2×.
- **The high-precision grid.** Both the exact DP and the tilted sweep scale linearly in node count
  (approximately 274×), so the ratio between them should hold. The saddlepoint's cost is independent
  of support width, which is the quantity most likely to change on the finer grid. Untested, and
  `PERFORMANCE.md` identifies this grid as where future work concentrates.

## Validation Gates

Every retained experiment should pass `cargo test --workspace`, the release
`golden_specprob` test, Clippy, and the full F13 profile. Do not use probability cutoffs, `f32`,
FMA, FFT convolution, or reordered summation in the bit-exact path.

An approximate estimator may violate those rules only in a separate, explicitly named API that the
bit-exact path does not call; it must document its measured error against the exact DP and assert
that error in a test.

Reproducing the trials requires the gitignored `validation/data/` and, where noted, the
MS-GF+-derived F13 golden:

```bash
cd rust
cargo run -p msgf-genfunc --example profile     --release   # stage breakdown
cargo run -p msgf-genfunc --example tailprune   --release   # pruning: bit-exactness and speedup curve
cargo run -p msgf-genfunc --example saddlepoint --release   # saddlepoint: accuracy and speedup
cargo run -p msgf-genfunc --example prunelab    --release   # exact ceiling + rejected policies
cargo test -p msgf-search --release --test golden_search -- --ignored
```

Run `tailprune` and `saddlepoint` on `worktree-genfunc-algo-speedups`, `prunelab` on
`worktree-genfunc-aggressive-prune`, and the profile comparison on `worktree-spec-tables-perf`.
