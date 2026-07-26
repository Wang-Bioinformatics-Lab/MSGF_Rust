# Research Trials

This directory holds detailed, reproducible reports for algorithm and performance experiments.
Reports should distinguish measured results from proposals, identify the branch or PR containing the
prototype, document fidelity tradeoffs, and include validation and reproduction commands.

## Reports

### Synthesis

- [Consolidating the measured speedups](consolidated-speedups.md) — **the current state.** What
  branch `perf/consolidated-speedups` merged, what each part was worth *in combination*, the three
  new optimizations, the defects adversarial review caught, and the fidelity evidence.
  314 → 572 spectra/s on the path `msgf search` drives (1.82×; 7.7× Java), byte-identical output.
- [Generating-function and scoring-pipeline optimization](generating-function-optimization.md) —
  DP profiling, exact sink and threshold pruning, aggressive-pruning limits, saddlepoint
  approximation, node-table optimization, rejected approaches, and remaining opportunities.
  Predates the consolidation; its per-experiment baselines differ from each other.

### Per-trial data

These keep the instrument output and the reasoning steps — including the hypotheses that turned out
wrong — behind the conclusions in the synthesis. The point is that a conclusion can be re-checked
without re-running the investigation.

- [Per-spectrum node tables (`spec-tables`)](spec-tables-node-scoring.md) — inner-loop census, the
  replica-timing step that located the real cost, the node-score cache, the ion-major sweep, and the
  allocation hypothesis that was wrong. 820 ms → 111 ms, bit-exact.
- [Saddlepoint tail approximation: tuning data](saddlepoint-tuning-data.md) — how it went from 1.26×
  to 3.2×, the `Z(0)` bias and its fix, the solve-tolerance sweep, accuracy per decade of tail
  depth, and two test-calibration mistakes.
- [How hard can the DP prune?](dp-pruning-limits.md) — the exact bound's derivation and why it is
  already optimal for exactness, the per-edge cost floor that caps every cell-pruning idea at ~6×,
  and full instrument output for three rejected aggressive prunes (top cap, Chernoff tilted trim,
  dead-node sweep) plus two changes that measured as exactly zero.

### Cross-cutting

- [Measurement traps in this repository](measurement-traps.md) — F13's degeneracy as a benchmark,
  cells ≠ time, allocation counters, relative vs absolute tail depth, peptide-shaped synthetic
  graphs, `to_bits` as the bit-exactness gate, and the worktree/gitignore hazard. Read before
  designing a new measurement.

`ALGORITHMIDEAS.md` at the repository root is the concise index; keep detailed trial narratives here.
