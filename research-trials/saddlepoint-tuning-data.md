# Trial data: saddlepoint tail approximation

Raw tuning and accuracy data behind
[generating-function-optimization.md](generating-function-optimization.md)
§*Approximate Routes*. The synthesis states the conclusions; this file keeps the numbers that
produced them, the two calibration mistakes made along the way, and what would be needed to go
faster.

**Branch:** [`worktree-genfunc-algo-speedups`](https://github.com/Wang-Bioinformatics-Lab/MSGF_Rust/tree/worktree-genfunc-algo-speedups)
(`7c8a314`), module `msgf_genfunc::saddle`, example `saddlepoint`.

**Workload:** the 1,253 F13 spectra that have a golden MS-GF+ RawScore (2,506 graphs), thresholds =
MS-GF+'s own observed top-hit RawScore. Exact DP over the same graphs = 2,946 ms.

## How it got from 1.26× to 3.2×

Each row is a real measurement, in the order taken. Accuracy is `log10(approx/exact)` against the
exact DP.

| # | change | sweeps/spectrum | time | vs exact | median | ≤0.05 |
|---|---|--:|--:|--:|--:|--:|
| 1 | first working version | 10.6 | 2340 ms | 1.26× | −0.002 | 96.1% |
| 2 | hoist per-sweep range scan; loosen solve; drop `θ=0` sweep | 6.0 | 1121 ms | 2.63× | **+0.030** | **79.8%** |
| 3 | restore `Z(0)` from a global mass table | 5.9 | 1090 ms | 2.75× | −0.002 | 96.1% |
| 4 | solve tolerance 0.02 → 0.25 | 5.0 | 934 ms | **3.17×** | −0.002 | 96.1% |

Row 2 is the instructive one: dropping the `θ = 0` sweep bought speed and cost a **systematic
+0.03 log10 bias** — 96.1% → 79.8% inside the project's 0.05 bar. The leading deep-tail term does not
need `K(0)` (it cancels: `Z(0)·φ(ŵ) = e^{K(θ̂)−θ̂t}/√(2π)` identically), but the Mills-ratio
correction does, and that correction is not negligible.

Row 3 recovers it for free via the observation that **`Z(0)` is spectrum-independent**: at `θ = 0`
every tilt weight is 1, the recursion degenerates to `z[m] = Σ_aa p_aa·z[m − mass(aa)]`, the spectrum
drops out entirely and so does the cleavage factor. One mass-composition table per run replaces a
sweep per graph. Accuracy returns to the row-1 figures at row-2 speed.

### Solve-tolerance sweep

The tail is stationary in θ at the saddlepoint, so `K'` precision is cheap to trade:

| tolerance on \|K′ − t\| | sweeps/spectrum | time | vs exact | ≤0.05 | ≤0.30 |
|---|--:|--:|--:|--:|--:|
| 0.02 | 5.9 | 1090 ms | 2.75× | 96.09% | 100% |
| **0.25** | **5.0** | **934 ms** | **3.17×** | **96.09%** | **100%** |
| 1.0 | 4.4 | 827 ms | 3.59× | 95.93% | 99.92% |

0.25 is free — identical accuracy to a 1e-9 solve. 1.0 starts to leak (the first point outside 0.30),
so 0.25 is the setting in `SOLVE_TOL`.

## Accuracy vs tail depth

Measured over 865 tail points on synthetic peptide-shaped graphs (`saddlepoint_tracks_the_exact_tail`),
binned by **relative** depth `tail / Z(0)`:

| relative depth | n | median \|log10\| | worst |
|---|--:|--:|--:|
| 1e0 – 1e-2 | 86 | 0.0012 | 0.0087 |
| 1e-2 – 1e-4 | 175 | 0.0018 | 0.0175 |
| 1e-4 – 1e-6 | 126 | 0.0039 | 0.0212 |
| 1e-6 – 1e-9 | 152 | 0.0094 | 0.0267 |
| 1e-9 – 1e-12 | 121 | 0.0103 | 0.0962 |
| 1e-12 – 1e-16 | 123 | 0.0397 | 0.3234 |
| beyond 1e-16 | 82 | 0.2294 | 1.0129 |

Overall median +0.0004 — no systematic bias. Degradation past 1e-12 is the threshold approaching the
maximum achievable score, where too few paths carry the tail for a normal-based inversion. Real
SpecEValues sit near 1e-6 relative depth.

`DEPTH_BOUNDS` in the test encodes this shape as the regression bound rather than a flat tolerance,
because a flat bound would either be a lie at depth or useless near the mean.

## Two calibration mistakes (both in the *test*, not the method)

Worth recording — both produced convincing-looking failures of correct code.

1. **Synthetic graphs with unit mass steps.** A 260-node graph with steps 1..6 gives source→sink
   paths of ~250 residues; path probabilities (`0.05^250`) underflow `f64` and `cumulants` returned
   `None`. Nothing was wrong with the estimator — the test graph was not peptide-shaped. Fixed by
   using steps `[9, 11, 14, 17, 22, 29]` over 260 nodes, giving ~10–30-residue paths. **Any synthetic
   de novo graph used for numerical testing must have a realistic mass-step-to-extent ratio.**
2. **Selecting tail points by absolute probability.** Filtering to `exact < 0.05` looked like "the
   tail," but these distributions are *sub*-probability measures whose total `Z(0)` can itself be
   ~1e-9 — so an absolute 2.3e-9 was the **70th percentile**, not a tail, and the estimator correctly
   returned ≈`Z(0)`. Tail depth must always be taken relative to `Z(0)`. This is the same trap for
   anyone reading a SpecEValue: an absolute 1e-9 is nowhere near 1e-9 relative depth.

## What would make it faster

Currently ~2.5 tilted sweeps per graph, near the Newton floor (one to evaluate, one to confirm). Per
sweep it is ~6.7× cheaper than the exact DP, so the ceiling at one sweep is ~6.7×.

- Fuse 3–4 θ values into one AVX sweep and recover `K'`/`K''` by finite differences instead of
  propagating `(Z, Z', Z'')` — fewer, wider sweeps, potentially ~1 sweep-equivalent.
- The per-sweep cost is ~4.5× less efficient per flop than the DP's AVX `axpy` (gather from the tilt
  table, three dependent accumulator streams, scalar reduction), so there is headroom in the sweep
  itself independent of sweep count.

## Reproducing

```bash
git checkout worktree-genfunc-algo-speedups
cd rust
cargo run -p msgf-genfunc --example saddlepoint --release   # accuracy + speed table
cargo test -p msgf-genfunc --release --lib -- --nocapture   # depth-banded bounds
```
