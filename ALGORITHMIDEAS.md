# Algorithm Ideas

This file is the index for algorithm and performance research. Detailed measurements, rejected
approaches, implementation notes, and reproduction commands belong under [`research-trials/`](research-trials/).

The **plan** these trials feed is [`plans/PLAN3.md`](plans/PLAN3.md) — spectral p-value acceleration
(5–10× on the significance stage): what to build, in what order, and the acceptance gates. This file
records what was measured; PLAN3 records what to do about it.

## Generating Function and Scoring Pipeline

See the [full optimization trial report](research-trials/generating-function-optimization.md), the
per-trial data behind it ([node tables](research-trials/spec-tables-node-scoring.md),
[saddlepoint tuning](research-trials/saddlepoint-tuning-data.md)), and
[measurement traps](research-trials/measurement-traps.md) — read that last one before designing a new
benchmark here.

**Consolidation status (2026-07-25).** Everything marked *merged* below lives on branch
`perf/consolidated-speedups`, together and green. On the 1,406-spectrum F13 profile that branch
measures **314 → 436 spectra/s** unpruned and **314 → 572 spectra/s** on the threshold-driven path
`msgf search` actually uses (**1.82×**). `msgf search` output on F13 is **byte-identical** to
pre-change `main` (same sha256), the F13 end-to-end oracle is unchanged at 1161/1161 exact
RawScore and DeNovoScore, and `golden_specprob` is 30/30.

| Idea or trial | Fidelity | Status | Key result |
|---|---|---|---|
| Remove redundant DP range and bounds work | Exact | **merged** (`c8400e3` on main) | ~3% of the DP on its own; PR #9's other half is subsumed below |
| Prune nodes outside all sink paths | Exact | **subsumed** by threshold pruning | `max_remaining == i32::MIN` performs the same sink-ancestor elimination |
| Threshold-aware tail pruning | Exact above declared threshold | **merged**, and now *wired into `search` and `rescore`* | 1.31× DP on F13's own RawScores, 2.68× at DeNovoScore − 20; 2.01× fewer cells |
| Node-score cache and ion-major table sweep | Exact | **merged** | `spec-tables` 844 ms → 111 ms (7.6×) |
| Fuse the DP's two per-edge passes | Exact | **merged** | DP 1.04× unpruned, ~1.06× pruned — attacks the 2.9 ns/edge floor |
| Reusable graph buffers + `edge_prob` as an aa index + hoisted edge-score constants | Exact | **merged** | `graph-build` 558 → 238 ms (2.3×); allocation 682 MB → 3.4 MB per run |
| Saddlepoint tail inversion | Approximate, opt-in | Branch `worktree-genfunc-algo-speedups` only | 3.2× DP speed; 96.1% within 0.05 log10 |
| Tiered saddlepoint then exact DP | Mixed | Proposed | Approximate screening with exact evaluation near decisions |
| Sparse/dense score distributions | Potentially exact | Proposed | Measure early-node sparsity before implementing |
| Reuse reachability across isotope sinks | Exact | Proposed | Avoid rebuilding nearly identical reverse bounds |
| Two-sided corridor + exact early-exit rejection | Exact | Proposed — see `plans/PLAN3.md` §5.2 | The only route past the ~6× cell-pruning ceiling |

## Important Conclusions

- PR #9 and threshold pruning overlap: `max_remaining == i32::MIN` performs the same sink-ancestor
  elimination, so their gains are not additive.
- Exact cell pruning reaches a fixed per-edge floor. At DeNovoScore − 5, 390× fewer cells produced
  only a 6.2× DP speedup; further work should reduce edge-visit cost or optimize graph construction.
- **A measured optimization is worth nothing until it has a caller.** The tail prune sat on a branch
  with no call site: `search` built the generating function *before* scoring candidates, so it had no
  threshold to offer. Reordering those two halves — which is free, since the DP depends on nothing
  the candidates produce — is what turned a 1.31× DP result into a 1.82× product.
- The DP is now **83.5%** of the per-spectrum stage and graph build 9.6%, so the Amdahl budget in
  `plans/PLAN3.md` §3.3 has tightened: non-DP work can no longer fund much. Past ~6× the DP needs a
  gate that skips it entirely, not a cheaper cell.
- Top caps, Chernoff trimming, extra dead-node sweeps, FFT convolution, and score-lattice coarsening
  were evaluated and rejected. The report records why and, where available, the measured error.
- Approximate algorithms must remain explicitly named and isolated from the bit-exact default path.

## Validation Standard

Retained exact changes must pass the full workspace suite, release golden SpecEValue tests, Clippy,
and the F13 profile. Approximate trials must publish accuracy against the exact DP and must never be
selected silently by search or CLI code.
