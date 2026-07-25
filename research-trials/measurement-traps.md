# Measurement traps in this repository

Cross-cutting lessons from the optimization trials. Each one produced a wrong conclusion, a wasted
experiment, or a convincing-looking failure of correct code. They are collected here because they are
properties of *this codebase and corpus*, so they will recur.

## 1. F13 is a worst case for anything threshold-driven

F13 identifies essentially nothing — MS-GF+'s own top hits there are ~50% decoy (see
`f13-degenerate-fdr-oracle`). Consequences for benchmarking:

| statistic | value |
|---|--:|
| DeNovoScore, median | 66 |
| MSGFScore (RawScore), median | 3 |
| DeNovoScore − RawScore gap, median | **69** |
| same gap, best 50 PSMs by SpecEValue | 25 |
| same gap, best 200 | 33 |

Any optimization whose payoff scales with how close the RawScore is to the DeNovoScore — tail
pruning above all — will look far worse on F13 than on a corpus with real identifications. Quote the
*curve* against the gap, not the single F13 aggregate, and say which regime a number came from.

F13 remains the right corpus for **fidelity** (it is the golden oracle); it is a poor one for
**judging threshold-sensitive speedups**.

## 2. Cell counts are not time

Tail-threshold pruning removed **2.11× the distribution cells** and bought **1.42× the time**.
Narrow `axpy` slices vectorize worse than wide ones and per-node overhead becomes a larger share, so
the AVX kernel's efficiency falls as the thing being pruned gets smaller. The same gap appears in the
aggressive-pruning trials, where 390× fewer cells produced only ~6× less time.

Always report both, and never convert one into the other.

## 3. Allocation counters point at anomalies, not bottlenecks

Twice now the allocator census flagged something real that was not the cost:

- The original DP port assumed allocation was the bottleneck; removing essentially all of it bought
  ~10% (`PERFORMANCE.md`).
- The node-score cache showed reallocs 0 → 8.9/spectrum and +22 MB traffic. `Vec::with_capacity`
  fixed the counters and moved the time by ~5 ms of 108 ms. The real cost was re-resolving a
  string-keyed table row on each of 151 rank bins.

Measure the hypothesis, not the counter.

## 4. A census of "wasted work" does not imply the waste is the cost

`node_score` discards 75% of its inner iterations, suggesting a 4× ceiling. But a replica of that
whole loop *including* peak lookups ran in 234 ms against the real 791 ms — ~70% of the cost was in
neither. Removing the dead iterations first would have chased a third of the available win.

Before optimizing a loop, time a replica of it. The gap between the replica and the real thing is
where the cost actually is.

## 5. Tail depth is relative to `Z(0)`, never absolute

The score distributions here are **sub-probability measures**: `Z(0)` is the probability that a
random peptide hits the precursor mass at all, and it can itself be ~1e-9. So an absolute
SpecEValue of 1e-9 may be the 70th percentile of its own distribution, not a deep tail.

This broke a saddlepoint test that filtered tail points by absolute probability, and it is the same
trap when reading a SpecEValue: `tail / Z(0)` is the only scale-free measure of "how far out."

## 6. Synthetic de novo graphs must be peptide-shaped

A 260-node test graph with mass steps 1..6 has source→sink paths of ~250 residues, whose path
probabilities (`0.05^250`) underflow `f64`. Correct code returns `None` and the test fails for
reasons unrelated to what it is testing.

Keep the mass-step-to-extent ratio realistic — steps like `[9, 11, 14, 17, 22, 29]` over 260 nodes
give ~10–30-residue paths, matching tryptic peptides. Real residue masses are 57–186 against
peptide masses of ~1,000–3,000.

## 7. Bit-exactness has a cheap, strong test — use it

Every exact optimization here was verified with `f64::to_bits` / `f32::to_bits` equality against the
unoptimized path, not with an epsilon. Three of them additionally produced a **byte-identical**
`msgf search` output TSV on F13. When a change is supposed to be exact, an approximate comparison
tests nothing:

- reorderings that preserve summation order are bit-identical and should be asserted as such;
- if a change cannot pass a `to_bits` test, it is a fidelity change and belongs in a separate,
  explicitly named API (see the saddlepoint module).

Keep the unoptimized path callable so the test has something to compare against. Both scorer tests
do this via a private `..._with(cache: Option<&_>)` seam.

## 8. Worktrees and the data-absence contract

`validation/.gitignore` uses directory-form rules (`data/`, `rawscore/`, `spectra/`). These do **not**
match symlinks. When working in a `.claude/worktrees/` checkout, the usual move is to symlink
`validation/data` and the UC-derived goldens from the main checkout — and those symlinks show up as
untracked, committable files.

Remove them before committing (`git status --short` should show only source files), and never
`git add -A` in such a worktree. Committed symlinks would not leak UC-licensed bytes, but they
violate the contract in `CLAUDE.md` and produce broken paths for everyone else.

Also note the main checkout may have a concurrent session editing the same files — re-read before
writing shared documents.
