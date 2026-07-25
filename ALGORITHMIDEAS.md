# Generating-Function Algorithm Ideas

This document collects exact or potentially exact ways to reduce the generating-function dynamic
program. Numerical compatibility with MS-GF+ remains the primary constraint: optimizations must
preserve edge order, floating-point operations, the complete score distribution, DeNovoScore, and
SpecEValue unless an explicitly narrower API is introduced.

## Experimental Reference

Draft PR [#9 — prune dead DP subgraphs](https://github.com/Wang-Bioinformatics-Lab/MSGF_Rust/pull/9)
is a reference implementation, not intended for direct merging. It explores two complementary
optimizations:

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

## Further Exact-Pruning Ideas

### Threshold-aware score pruning

Search currently computes the full generating function before candidate RawScores. A specialized
search path could instead:

1. score candidates and retain the top matches;
2. determine the lowest RawScore whose tail probability is needed;
3. compute an optimistic maximum remaining score from every node to a sink;
4. discard a state only when `current_score + max_remaining < required_score`.

This can preserve requested upper-tail probabilities, but it cannot implement today's arbitrary
`GenFunc::spectral_probability` queries or derive DeNovoScore from a truncated distribution.
DeNovoScore would need a separate exact maximum-path calculation, and rescoring would require
grouping PSMs by spectrum before building the bounded distribution.

### Reuse reachability across isotope candidates

The `-ti 0,1` graphs share one edge structure and have adjacent sinks. Investigate caching or
incrementally updating the reverse-reachability mask rather than rebuilding it for each candidate.
Measure the mask-building cost first; it is small relative to convolution.

### Hybrid sparse/dense distributions

Early nodes may contain many zero score cells before distributions become dense. A sparse
representation could avoid those operations and switch permanently to the current contiguous arena
above a measured density threshold. Entries must retain score order, and conversion must not change
addition order.

## Validation Gates

Every retained experiment should pass `cargo test --workspace`, the release
`golden_specprob` test, Clippy, and the full F13 profile. Do not use probability cutoffs, `f32`,
FMA, FFT convolution, or reordered summation in the bit-exact path.
