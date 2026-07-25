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

| Idea or trial | Fidelity | Status | Key result |
|---|---|---|---|
| Remove redundant DP range and bounds work | Exact | Draft PR [#9](https://github.com/Wang-Bioinformatics-Lab/MSGF_Rust/pull/9) | Combined with sink pruning: 19.1% faster DP |
| Prune nodes outside all sink paths | Exact | Draft PR #9 | 21% fewer cells; 12.9% faster pipeline with both PR changes |
| Threshold-aware tail pruning | Exact above declared threshold | Branch `worktree-genfunc-algo-speedups` / `worktree-genfunc-aggressive-prune` | 1.3× on poor F13 matches; up to 6.2× near DeNovoScore |
| Node-score cache and ion-major table sweep | Exact | Branch `worktree-spec-tables-perf` | `spec-tables` 820 ms → 111 ms; throughput 314/s → 371/s |
| Saddlepoint tail inversion | Approximate, opt-in | Experimental branch only | 3.2× DP speed; 96.1% within 0.05 log10 |
| Tiered saddlepoint then exact DP | Mixed | Proposed | Approximate screening with exact evaluation near decisions |
| Sparse/dense score distributions | Potentially exact | Proposed | Measure early-node sparsity before implementing |
| Reuse reachability across isotope sinks | Exact | Proposed | Avoid rebuilding nearly identical reverse bounds |

## Important Conclusions

- PR #9 and threshold pruning overlap: `max_remaining == i32::MIN` performs the same sink-ancestor
  elimination, so their gains are not additive.
- Exact cell pruning reaches a fixed per-edge floor. At DeNovoScore − 5, 390× fewer cells produced
  only a 6.2× DP speedup; further work should reduce edge-visit cost or optimize graph construction.
- Top caps, Chernoff trimming, extra dead-node sweeps, FFT convolution, and score-lattice coarsening
  were evaluated and rejected. The report records why and, where available, the measured error.
- Approximate algorithms must remain explicitly named and isolated from the bit-exact default path.

## Validation Standard

Retained exact changes must pass the full workspace suite, release golden SpecEValue tests, Clippy,
and the F13 profile. Approximate trials must publish accuracy against the exact DP and must never be
selected silently by search or CLI code.
