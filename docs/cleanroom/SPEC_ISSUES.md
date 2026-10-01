# SPEC_ISSUES: points in SPEC.md (v1, 2026-09-30) that needed a decision

None of these changed a number on any test vector. Every stage matched the recorded oracle
values on the first complete build, so no black-box oracle experiment was needed to resolve
behaviour. Each item records how the text was read.

1. **§8.7: where the pruning cut is clamped.** "Clamp the cut at the best full-path score" does
   not say whether this is per sink or over all sinks. Implemented per sink:
   `cut_p = min(θ − κ⁺, best_p)`. Either reading leaves every tail value at σ ≥ θ and the
   DeNovoScore unchanged. A sink whose best score is below θ − κ⁺ contributes nothing to any
   queried tail anyway. When cleavage is disabled the shift is 0 (the cut is θ).

2. **§3.5: `d` exactly equal to `±τ_d`.** The pseudocode tests `−τ_d < d < τ_d` (match) and
   then `d > τ_d` (stop). A value exactly at `+τ_d` or `−τ_d` hits neither branch. Implemented
   literally: the scan continues. This case did not occur in the vectors.

3. **§7.2 with no anchor type at all.** The edge direction is defined only for a suffix or
   prefix anchor type. With no anchor type (a model whose segment partitions have no ion types),
   every α(k > 0) is −1. Implemented as the suffix direction. The edge values do not depend on
   the direction in that case, because every edge takes a constant. This cannot happen with the
   bundled model.

4. **§8.2 and §8.3: alphabet entries with nominal mass ≤ 0.** These are not defined, and the
   DP would not terminate. Such entries are dropped from the graph. No real residue has one.

5. **Rescore peptide syntax (USAGE.md): how flanks are detected.** The text says flanks
   `X.SEQ.Y` are optional, but a modification value also contains `.`. Implemented: the string
   has flanks when it is at least 5 characters long, its second character is `.`, and its
   second-to-last character is `.`. A bare core sequence cannot trigger this, because a core
   never has `.` second (it starts with a residue letter), except in a pathological
   single-residue-plus-mod string. All 17,724 + 38 test PSMs parse as the oracle parsed them.

6. **SpecEValue text output: last printed digit.** `msgf rescore` prints 7 significant digits.
   The first build summed the tail from the top of the support downward. With that order, the
   uniform-alphabet hela sets had 3, 5 and 1 PSMs that differed from the recorded files by one
   unit in the 7th digit. These values are dyadic, so their decimal expansion ends in …5
   exactly, and a one-ulp difference in the f64 sum flips the printed rounding. The difference
   was already well inside §11's tolerance (max |Δlog10| 3.5e-7). SPEC does not fix the
   summation order. Summing from RawScore upward (ascending score) makes all seven
   `specprob-rescore` outputs **byte-identical** to the recorded oracle files, so that order is
   the one kept.

7. **§11.3 "ρ within 1 ulp" and §11.2 "≤ 1 ulp on rewritten m/z".** Both were met exactly
   (0 ulp) on every recorded spectrum. The harness still applies the 1-ulp tolerance.

8. **Search-mode ID check (§11.6).** The q-value procedure is not specified. The harness takes
   the top-1 PSM per spectrum (lowest SpecEValue, file order breaking ties), sorts ascending, and
   computes q = decoys/targets with a monotone running minimum. It applies the identical
   procedure to both score columns. The 25 `context_ambiguous` rows are excluded from both sides.
