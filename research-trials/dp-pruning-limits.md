# How Hard Can the Generating-Function DP Prune? — trial report

Branch `worktree-genfunc-aggressive-prune` (commits `b31c6d1`, `76ad840`; not pushed at the time of
writing). Harness: `rust/crates/msgf-genfunc/examples/prunelab.rs`.

The question this trial set out to answer was "can we be more aggressive on the DP pruning?" The
answer is **no, and the reason is not the one that was expected**. Exact threshold pruning is worth
keeping (1.31× at realistic thresholds, up to 6.17× at optimistic ones, bit-exact). Every attempt to
prune *harder* than the exact bound lost, and separately it emerged that cell count has stopped
being what the DP costs, which caps the whole family of ideas.

This report keeps the derivations and the raw measurements. `ALGORITHMIDEAS.md` carries the summary.

---

## 1. Setup

| | |
|---|---|
| Corpus | F13, 1,406 spectra after the 200–6000 nominal-mass filter |
| Model | `validation/data/models/HCD_HighRes_Tryp.param` (MS-GF+'s own, not the bundled one) |
| Grid | nominal (`0.999497`) |
| Isotope errors | `-ti 0,1` → 2 graphs per spectrum |
| Cleavage | `credit +2`, `penalty −11`, `prob_cleavage_sites 0.10` |
| Amino acids | 20 standard + oxidized M, all at `prob = 0.05` |
| Threading | single-thread |
| CPU | Intel Xeon E5-2667 v3 @ 3.20 GHz (AVX, no AVX2 path used), 16 MiB L3 |
| Toolchain | rustc 1.94.0, `opt-level = 3`, `lto = "thin"`, `codegen-units = 1` |

Per graph: 1,528 nodes (1,349 reachable), 29,576 edges, 136,351 distribution cells unpruned, mean
score-support width 101.

**Method.** `prepare()` does preprocessing, the scored spectrum, `tables()` and `build_reverse_graph`
once per spectrum, all outside the timed region. The timer covers only `recompute_node_scores` +
the DP call(s) + `merge_group` + one `spectral_probability`. Each (mode, threshold) is run three
times and the fastest kept. The reference tail for each spectrum and each threshold comes from a
`compute_into` pass over the same graphs, so error is measured against this build's own unpruned DP,
not against a stored golden.

**Thresholds.** Two families:

- *MS-GF+'s own RawScore* — the maximum `MSGFScore` per `ScanNum` from the frozen
  `validation/golden/iprg2013_F13.tsv`. This is the threshold a search actually knows once it has
  scored its candidates. 1,253 of the 1,406 spectra have one.
- *DeNovoScore − k* — sweeps match quality without depending on the golden. All 1,406 spectra.

**F13 is the worst case and should not be read as typical.** It identifies essentially nothing (its
own top hits are ~50% decoy), so its DeNovoScore − RawScore gap has median 69 — the far right of the
`− k` table. A corpus with real identifications sits nearer `− 20`.

---

## 2. What was implemented

All of it lives in `rust/crates/msgf-genfunc/`:

| Item | Where |
|---|---|
| `Prune` policy struct, `compute_tail_into`, `compute_tail_with` | `src/lib.rs` |
| `max_remaining` (backward integer sweep), `max_achievable` (forward) | `src/lib.rs` |
| Backward tilted sweeps, Newton saddlepoint solve | `src/tilt.rs` |
| Measurement harness | `examples/prunelab.rs` |
| `tail_prune_is_bit_exact`, `capped_prune_is_one_sided_and_certified`, `avx_matches_scalar_bitwise` | `src/lib.rs` tests |

```rust
pub struct Prune {
    pub threshold: i32,             // lowest score whose tail the caller will read
    pub cap: Option<i32>,           // lossy: score units retained above the threshold
    pub tilt: Option<(f64, f64)>,   // lossy: (θ, absolute error budget)
    pub skip_dead_nodes: bool,      // exact: forward max_achievable sweep
}
```

`Prune::exact(t)` is the only policy that survived. The other three fields are retained so this
report's tables can be regenerated; they are documented in-place as rejected.

---

## 3. The exact bound, and why it is already optimal *for exactness*

Let `P(m,s)` be the DP cell: the total probability-weight of source→`m` paths scoring `s`. Let
`Q_m(r)` be the weight of `m`→sink suffixes scoring at least `r`. The tail at `cut` decomposes as

```
tail(cut) = Σ_s P(m,s) · Q_m(cut − s)      (for any single node m on every path)
```

so cell `(m,s)` matters iff `Q_m(cut − s) > 0`, i.e. iff some suffix scores at least `cut − s`. With
`max_rem[m]` the best suffix score (one descending integer sweep of the CSR — edges are stored by
destination and `prev < i` always, so every node is final when read):

```
drop (m,s)  ⟺  s + max_rem[m] < cut  ⟺  s < cut − max_rem[m]
```

That is a **per-node score floor**, and it is *exactly* the set of cells that cannot contribute. No
exact rule can drop more. The retained cells are reached by the same multiply-adds in the same order
as in the unpruned DP, so they are bit-identical `f64`, not approximately equal.

Two details that make it work end to end:

- **DeNovoScore.** `max_rem[0]` is the source's best full-path score, so the same sweep that supplies
  the bound supplies the exact DeNovoScore; `GenFunc::max_score()` reads `max_rem[0] + credit`
  rather than the distribution. This is what later allowed the (lossy) top cap to be evaluated at
  all — a cap removes the maximum-score cells, so a distribution-derived DeNovoScore would break.
- **`cut` clamping.** `cut = min(threshold − max(credit, penalty), max_rem[0])`. The cleavage
  weighting shifts the merged distribution by `credit`/`penalty`, so the internal threshold is
  lowered by `credit` to keep the credited branch intact; clamping to the DeNovoScore keeps the
  maximum-score cell alive when a caller asks for a threshold no peptide can reach.
- **`valid_from`.** `GenFunc` records the lowest score it answers for; querying below it is
  debug-asserted. `merge_group` takes the `max` over members.

### Measured

Cells are per graph; "vs full" is the reduction against `compute_into`.

| Threshold | DP time | speedup | cells/graph | vs full |
|---|--:|--:|--:|--:|
| full (any) | 2,843 ms | 1.00× | 146,152 / 136,351 † | 1.00× |
| MS-GF+ RawScore | 2,165 ms | **1.31×** | 69,162 | 2.11× |
| DeNovoScore − 40 | 1,668 ms | 1.76× | 22,137 | 6.2× |
| DeNovoScore − 20 | 1,102 ms | 2.68× | 5,973 | 22.8× |
| DeNovoScore − 10 | 691 ms | 4.31× | 1,414 | 96.4× |
| DeNovoScore − 5 | 479 ms | 6.17× | 350 | 390× |

† 146,152 over the 1,253-spectrum golden subset, 136,351 over all 1,406. The RawScore row's `full`
baseline is 2,843 ms over 2,506 graphs; the `− k` rows' is ~2,950 ms over 2,812 graphs.

Run-to-run spread across five separate invocations was ~4% (e.g. the RawScore `exact` row measured
2,113 / 2,153 / 2,162 / 2,167 / 2,200 ms). Treat one significant figure of the speedup as real.

**Reconciling with the earlier 1.42×.** `worktree-genfunc-algo-speedups` measured the same prune at
1.4× on F13's own RawScores, and `measurement-traps.md` quotes that figure. Both are right: that
branch forked from `main` at `1b14c5f`, *before* the redundant-bounds fix (`c8400e3`, PR #9's first
optimization), so its unpruned baseline was slower and the prune looked better against it. This
trial runs on a `main` that already has the faster unpruned convolution, so the same prune shows
1.31×. The cell counts are identical across both (69,162), which is the tell that the algorithm did
not change — only what it is being compared against.

**Validation.** `tail_prune_is_bit_exact` compares `f64::to_bits` of the pruned and unpruned tails
over 40 pseudo-random signed-score DAGs × both cleavage settings × 33 thresholds each. The golden
suite is unchanged: `golden_specprob`, `golden_rescore` (30 PSMs), and the `#[ignore]`d
`golden_search` (1,161 PSMs, exact RawScore and DeNovoScore) all pass.

---

## 4. The finding that caps everything: cells are no longer the cost

At DeNovoScore − 5 the DP retains **390× fewer cells** and runs **6.2× faster**. Those numbers do not
reconcile unless there is a large width-independent term. There is:

```
DP, unpruned, per graph:   1,052 µs        (2,958 ms / 2,812 graphs)
DP, exact prune at −5:       170 µs        (  479 ms / 2,812 graphs)
```

170 µs per graph survives the near-total removal of convolution work. That is **16% of the unpruned
DP**, and it is entirely per-edge:

1. `max_remaining` — one descending pass over all 29,576 edges.
2. The DP's score-range pass — for every node not already excluded, a gather of `NodeDist` per
   incoming edge to compute `cur_min`/`cur_max`.

Both are O(edges) and independent of the score-support width. Two passes × 29,576 edges = 59,152 edge
visits in 170 µs → **2.9 ns per edge visit**, ≈9 cycles at 3.2 GHz.

Consequences:

- **No cell-pruning idea can beat ~6× on this graph shape.** The asymptote is the per-edge floor.
- At realistic thresholds the floor is already about a third of the pruned DP (170 of ~505 µs at the
  RawScore threshold), which is why 2.11× fewer cells buys only 1.31× time.
- The convolution itself, by subtraction, runs at ~0.3 ns per multiply-add unpruned (~1 cycle) and
  ~0.46 ns pruned — narrower slices vectorize worse, so cell reduction has *sub*-linear payoff even
  before the floor.
- This reframes the width-101 observation in the older profiling notes. Width is the opportunity in
  *arithmetic*, but arithmetic is no longer the whole cost.

The actionable follow-on is therefore not more pruning but cheaper edge visits — see §8.

---

## 5. Rejected: top cap on the score axis

**Idea.** Retain only `cap` score units above the threshold: `ceiling = cut + cap`, applied at every
node. The bound looks not merely safe but *tight*: for `s ≥ cut` the suffix requirement `r = cut − s`
is non-positive, so `Q_m(r) = Q_m(0) = B_0(m) ≤ ~1` — every completing suffix works. A discarded
cell can therefore add at most its own probability to the tail, and since `P(m,s)` decays fast above
the threshold, the discarded mass should be geometrically small.

**Why it fails.** *A path's partial score is not monotone.* This is the whole story. Node scores and
edge scores are both signed, and the C-terminal cleavage penalty is −11, so a path routinely climbs
tens of points above its final score partway along and comes back down. Cells at `cut + 30` in the
middle of the graph are not rare high-scoring outliers; they are ordinary paths mid-excursion, and
they carry real tail mass. Clipping them deletes it.

The bound was never wrong — the mass really is bounded by `P(m,s)`, and `err_bound` reported it
faithfully. The *premise* was wrong: `P(m,s)` above `cut` is not small at interior nodes.

### Measured, at the MS-GF+ RawScore threshold

| Top cap | DP time | speedup | cells/graph | vs full | mean \|log10\| | worst \|log10\| | certified rel. err (mean) |
|---|--:|--:|--:|--:|--:|--:|--:|
| none (exact) | 2,165 ms | 1.31× | 69,162 | 2.11× | 0 | 0 | 0 |
| +60 | 2,011 ms | 1.44× | 54,886 | 2.66× | 0.028 | 5.14 | 821 |
| +40 | 1,892 ms | 1.47× | 44,411 | 3.29× | 0.214 | 15.95 | 1,064 |
| +30 | 1,824 ms | 1.55× | 37,812 | 3.87× | 0.503 | 27.24 | 1,157 |
| +20 | 1,743 ms | 1.66× | 30,340 | 4.82× | 1.061 | 39.80 | 1,226 |
| +12 | 1,619 ms | 1.79× | 23,902 | 6.11× | 1.901 | 40.01 | 1,251 |

The project's tolerance is `|log10(rust/java)| ≤ 0.05`. Even `+60` — which buys 1.44× — has a worst
case of 5 log10. Nothing in this table is usable.

At optimistic thresholds the cap is harmless simply because it never binds: at DeNovoScore − 5 the
retained rows are already ~0.25 cells wide, so `+40` changes nothing at all (350 cells, 0 error) and
`+30` touches 4 cells per graph.

### The one genuinely good property

The failure is **self-reported**. The cap is one-sided by construction — probability is only ever
removed — so `p` is a lower bound, and `err_bound` accumulates a true upper bound on what was
discarded. The certified relative error at `+40` averages ~1,000, i.e. the run itself says "this
answer may be off by three orders of magnitude" without any reference DP to compare against. A lossy
scheme that cannot do this should not be merged.

---

## 6. Rejected: Chernoff (exponentially tilted) low-end trim

This was the most promising idea on the pre-existing list, and the one most worth writing down in
full, because the reason it fails is structural rather than incidental.

**Idea.** The exact floor asks *can* a cell reach the threshold — a worst case over the single best
remaining path, and a wildly improbable one. Ask instead how much a cell can *contribute*. Define the
backward tilted sum

```
B_θ(m) = Σ_{m→sink suffixes} weight · e^{θ·score}
       = Σ_{edges m→i} p_e · e^{θ·w_e} · B_θ(i),     B_θ(sink) = 1
```

with `w_e = nodeScore(i) + edgeScore(e)` (sink edges carry `edgeScore = 0`, as in the DP). Markov's
inequality on `e^{θ·score}` gives, for any `θ ≥ 0`,

```
Q_m(r) ≤ e^{−θ·r} · B_θ(m)
```

so cell `(m,s)` contributes at most `P(m,s) · e^{−θ(cut−s)} · B_θ(m)`. Discard from each node's low
end while the running sum of that bound stays inside a budget.

**Implementation notes** (`src/tilt.rs`):

- One descending sweep computes `B_θ` over the same CSR as `max_remaining`. `e^{θ·w}` is a lookup
  table over the graph's integer `w` range, so the sweep is multiply-add only — no `exp` per edge.
- Carrying `dB/dθ` and `d²B/dθ²` through the same recursion gives `K = ln B_θ(0)` and its first two
  derivatives (`B_θ(0)` is the whole graph's MGF). `solve_theta` Newton-solves `K'(θ) = cut` for the
  saddlepoint, tolerance 0.25 score units (the tail is stationary in θ at the saddlepoint, so a
  tighter solve buys nothing), warm-started from the previous graph's θ, clamped to `[0, 5]`.
- The cleavage weighting multiplies the MGF by a scalar, so its cumulants simply add.
- The budget is sized as `ε × tail_est`, `tail_est = e^{K − θ·cut} / (θ·√(2π·K''))`, the leading
  saddlepoint term. It is used **only** to size a budget, never reported.
- Per-node allowance is `remaining_budget / (n − i)`, so an early node with a fat low tail cannot
  spend it all.

**Error accounting is rigorous, not heuristic.** Write `P̃` for the computed (already-trimmed)
values; `P̃ ≤ P` elementwise since only non-negative mass is removed. By a hybrid argument over drop
events, the total deficit is `Σ_events P̃(m,s) · Q_m(cut−s)`, so accumulating the bound *using the
already-pruned values* is valid. Hence `tail ∈ [p, p + err_bound]`, one-sided.

**Trim after convolving, not before.** The row must exist for its discarded mass to be measured. This
is not the waste it looks like: the DP's work is `Σ_nodes retained_width × out-degree`, so a trimmed
row makes each of that node's ~21 successors cheaper and propagates onward through `cur_min` for
free. The loss is the ~1/21 spent at the node itself.

### Measured

At the MS-GF+ RawScore threshold, `exact` = 69,162 cells:

| ε | cells/graph | vs exact | DP time | vs full |
|---|--:|--:|--:|--:|
| 1e-3 | 67,341 | −2.6% | 3,089 ms | 0.92× |
| 1e-2 | 66,872 | −3.3% | 3,104 ms | 0.92× |
| 1e-1 | 66,112 | −4.4% | 3,072 ms | 0.93× |

**A 100× increase in error budget bought 1.8% more cells.** That flatness is the result.

**Why.** The tilt re-centres the distribution on the threshold, so the tilted measure's mass occupies
*the same width the distribution already has*. Dropping a relative ε of it means going out ~4σ, and
the row is only a few σ wide to begin with. Put differently: the exact floor `cut − max_rem[m]` is
already sitting near where the probabilistic bound would put it, so there is almost nothing between
"provably cannot contribute" and "contributes negligibly". The Chernoff bound's own looseness (a
factor of roughly `θ·√(2π·K'')` ≈ 10–15 at the saddlepoint, and worse away from it) eats the rest.

**And it is not free.** Isolating the sweep cost at DeNovoScore − 5: `cap+30` alone is 490 ms,
`cap+30 +tilt` is 2,006 ms, over 2,812 graphs at 3.84 sweeps/graph →
**~140 µs per derivative-carrying sweep, ~13% of the unpruned DP**. The Newton solve wanted 2.5–3.9
sweeps per graph depending on threshold (warm start helps but does not eliminate them). Net effect
across every threshold measured: 0.92×–1.03×, i.e. a loss.

**Budget mis-sizing.** In ~5 of 1,253 spectra at ε = 1e-3 the saddlepoint `tail_est` was far enough
off that the certified error exceeded ε (mean certified error over all spectra is dominated by these
few: 1,197). This is a real caveat for anyone reviving the idea — but it is secondary, since the
method loses on speed before accuracy is even reached.

---

## 7. Rejected: skipping dead nodes with a forward sweep

**Idea.** This is the only attempt that targets the §4 floor directly. Add a forward integer sweep
`max_achievable[m]` (best source→`m` score). Then `best[m] = ach[m] + max_rem[m]` is the best
full-path score *through* `m`, and any node with `best[m] < cut` lies on no path clearing the
threshold — so its ~21 edges never need to be visited at all.

The equivalence is exact: node `m` is empty under the existing rule iff `cut − max_rem[m] ≥ ach[m]+1`
iff `best[m] < cut`. And the optimal prefix path to `m` has `best[p] ≥ best[m]` at every node `p`
along it, so a retained node's `cur_max` really does reach `ach[m]`.

**Measured — 7–15% slower at every threshold:**

| Threshold | exact | exact + skip |
|---|--:|--:|
| MS-GF+ RawScore | 2,162 ms | 2,303 ms |
| DeNovoScore − 5 | 479 ms | 553 ms |
| DeNovoScore − 20 | 1,102 ms | 1,169 ms |
| DeNovoScore − 40 | 1,668 ms | 1,806 ms |

**Why.** The range pass on a dead node is already cheap. Its predecessors are themselves empty, so
each edge costs a `NodeDist` load and a predictable branch — it never reaches the arithmetic. The
forward sweep, by contrast, must do real work (`ach[prev] + node_score + edge_score`, compare) on
every *live* edge in the graph. Even at DeNovoScore − 5, where nearly every node is dead, the sweep
costs more than the passes it saves.

The lesson generalizes: the floor is not "visiting edges of dead nodes", it is "visiting edges at
all". Any fix has to make the *live* edge visit cheaper or eliminate a whole pass.

---

## 8. Neutral changes (kept, but honestly zero)

- **Kernel inlining.** The DP reached its `axpy` through an `unsafe fn` pointer chosen per call and
  invoked once per edge. Replaced with a `#[target_feature(enable = "avx")]` monomorphization of the
  whole DP body, dispatched once per `compute`, so the kernel inlines into the edge loop. Measured
  difference: **none** (2,913 vs 2,958 ms unpruned, 482 vs 505 ms pruned — inside run-to-run
  spread). LLVM was already devirtualizing the call, since the pointer's provenance is visible.
  Kept because it is strictly less indirection, and `avx_matches_scalar_bitwise` now pins the AVX and
  portable bodies bit-identical over 25 random graphs.
- **No-clip fast path under pruning.** A raised floor forces the clipping form of the convolution,
  but only at nodes where the floor actually moved. Tracking that per node and dispatching to the
  unclipped form otherwise measured neutral (2,162 / 2,167 ms vs 2,113–2,200 before). Kept on
  principle: it makes the pruned path degenerate exactly to the unpruned one when the threshold does
  not bite.

Both are recorded because "we tried it and it did nothing" is worth as much as a win, and because
the first one is the obvious first guess for anyone looking at the §4 floor.

---

## 9. Where the remaining opportunity is

Given §4, ideas that narrow the score axis are exhausted. What is left attacks the per-edge visit:

1. **Fuse the DP's two passes over each node's edges.** The score-range pass and the convolution pass
   each gather a `NodeDist` per edge and each re-derive `score_diff = node_score + edge_score`.
   Caching ≤32 resolved descriptors `(src_ptr, len, min_score, score_diff, prob)` on the stack during
   the first pass would remove the second gather entirely. Untested; this is the most direct attack
   on the 2.9 ns/edge figure.
2. **`edge_prob` as a per-amino-acid lookup.** It is one `f64` per edge holding one of ~21 distinct
   values — ~236 KB per graph — read by the DP twice *and* by `max_remaining` *and* by any tilted
   sweep. Replacing it with a small index removes a memory stream from all of them.
3. **Wire `compute_tail_into` into search.** Nothing calls it on this branch. The reordering that
   gives search a threshold (score candidates first, then build the generating function) exists on
   `worktree-genfunc-algo-speedups` and was not ported. Until it is, none of §3's speedup is
   realizable in the product.
4. **Rescore is still unaddressed.** `msgf rescore` caches one generating function per
   `(scan, charge)` across several PSMs, so it needs the per-key minimum RawScore as its threshold.

Item 3 is the one that turns this trial into a shipped improvement; 1 and 2 are the ones that would
raise the ceiling it runs into.

---

## 10. Reproduction

Requires the gitignored `validation/data/` and the MS-GF+-derived F13 golden.

```bash
git checkout worktree-genfunc-aggressive-prune
cd rust

# every table in §3, §5, §6, §7 (default thresholds: RawScore, DeNovo −5/−10/−20/−40)
cargo run -p msgf-genfunc --example prunelab --release

# a subset, faster
cargo run -p msgf-genfunc --example prunelab --release -- 5 20

# bit-exactness, one-sidedness of the lossy policies, AVX/scalar parity
cargo test -p msgf-genfunc --lib

# fidelity gates
cargo test --workspace --release
cargo test -p msgf-search --release --test golden_search -- --ignored
```

The harness prints, per mode: DP wall-time, speedup against the unpruned DP, retained cells per
graph, mean and worst `|log10|` deviation from the unpruned tail, the run's own accumulated certified
relative error, and how many spectra exceeded 1e-3 of it. `DeNovoScore mismatches: 0` on every row is
the check that pruning never disturbed the maximum-score path.
