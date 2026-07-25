# PLAN3 — Spectral p-value acceleration

Execution plan for making the **spectral p-value / SpecEValue stage 5–10× faster** without giving up
the fidelity contract. Strategy context is `PLAN.md` §4 (why the generating function is the whole
point); measured evidence is `PERFORMANCE.md` and `research-trials/`.

**Status: design doc, not started.** Several of its workstreams already have *measured prototypes*
on branches (§2) — consolidating those is task B1 and the fastest route to the first milestone.

**Read `research-trials/measurement-traps.md` before designing any benchmark or test in this plan.**
It records the corpus- and codebase-specific ways measurement here has already produced wrong
conclusions; §3, §4 and §12 below are written to obey it.

**Scope constraint (inherited from the brief this plan is written from):** do not optimize candidate
generation or shrink the candidate set. Work begins once peptide-spectrum RawScores are available.
The targets are the spectrum-specific null-tail calculation, decision-oriented preclassification,
and the final p-value or an equivalent calibrated confidence statistic.

---

## 1. Goal and scope

**Primary success criterion:** 5–10× end-to-end acceleration of the per-spectrum significance stage
on a representative search, **preserving identification power at 1% FDR**, with an explicit,
measured bound on the false-rejection rate of every reject-only shortcut.

**In scope** — everything that runs once per spectrum to turn a null model into a tail probability:

- per-spectrum node tables (`ScoredSpectrum::tables`),
- de novo graph construction (`graph::build_reverse_graph`) — see §3.3, it is *required* for 5×,
- the score-distribution DP (`compute_into`, `merge_group`),
- tail lookup (`ScoreDist::spectral_probability`) and how `search`/`rescore` drive it,
- gates, bounds, estimators and the cascade that decides which of them runs,
- the statistic itself, where an alternative is evaluated on FDR power rather than on agreement.

**Out of scope** — peptide index construction, digestion, modification enumeration, candidate
selection, RawScore computation, decoy database generation, q-value arithmetic (`PLAN2.md`).

**Non-goals.** Replacing the exact generating function as the default. Shipping a learned or
sampled p-value that the user did not ask for. Any change that makes the exact path
non-reproducible.

---

## 2. What is already measured (read before proposing anything)

The brief this plan derives from is generic. This repository has already run a substantial part of
it. Do not re-derive these.

| Brief method | Status here | Evidence |
|---|---|---|
| Prune nodes outside all sink paths | **Implemented, measured** | Draft PR #9: DP 19.1% faster, 21% fewer cells |
| Remove redundant DP range work | **Merged** (`c8400e3`) + PR #9 | `research-trials/generating-function-optimization.md` |
| Threshold-aware (one-sided) score pruning | **Implemented, bit-exact, measured** | **1.31×** on F13's own RawScores against a `main`-like baseline, 2.68× at DeNovoScore−20, 6.17× at −5 (`dp-pruning-limits.md`). The older 1.4× is the same code against the pre-`c8400e3` baseline. |
| Saddlepoint / cumulant DP | **Implemented, measured** | `msgf_genfunc::saddle`: **3.2×** vs exact DP, 96.1% within 0.05 log10; tuning + error tables in `research-trials/saddlepoint-tuning-data.md` |
| Per-spectrum node-table cost | **Fixed** | `worktree-spec-tables-perf`: 820 ms → 111 ms; 314 → 371 spectra/s |
| Two-sided ambiguity corridor (`L_v`, `Z_v`) | Not started — §5.2 | the implemented prune is the `U_v` half only |
| Exact early-exit rejection | Not started — §5.2 | needs `Z_v`; the highest-value missing exact piece |
| Shared multi-sink generating function | **Partly inapplicable** — §5.3 | node scores are candidate-dependent; only structure/bounds can be shared |
| Sparse/dense storage, O(1) tails | Not started — §5.4 | `spectral_probability` re-sums the tail per query |
| Certified coarse-score refinement | Not started — §6 | distinct from the *uncertified* coarsening already rejected |
| Conditional-null / tilted-IS sampling | Not started — §7 | |
| Streaming empirical null | Not started — §8 | cheapest possible gate; see §3.4 for why that matters |
| Distilled survival model | Not started — §8 | research track, off by default |
| Max-score / LR / PEP statistics | Not started — §9 | parallel track, judged on FDR power only |

Also already **evaluated and rejected with reasons** — see §14 before spending time on them: top-cap
truncation, Chernoff/tilted low-end trimming, an extra dead-node sweep, FFT convolution along the
score axis, uncertified score-lattice coarsening, contraction of score-inert graph regions, and
mass-axis reachability pruning on the nominal grid.

---

## 3. What the target number has to mean

### 3.1 The metric

**Unit of measurement:** single-thread throughput of the per-spectrum significance stage
(preprocess → scored spectrum → node tables → graph → DP → tail), in spectra/s, from
`cargo run -p msgf-genfunc --example profile --release`. Report median / p90 / p99 per-spectrum
latency alongside the mean, plus cells, edge visits, graph passes, and allocation volume.

**Report cells and time separately and never convert one into the other** (traps §2): pruning has
already shown 2.11× fewer cells buying 1.42×, and 390× fewer cells buying ~6×, because narrow `axpy`
slices vectorize worse. Allocation counters are anomaly detectors, not cost estimates (traps §3).

**Do not** quote end-to-end `msgf search` wall-clock as the headline. On F13 that run is
index-build dominated (~5 s of ~6 s, `PERFORMANCE.md`), so a 10× significance stage moves it by
under 15%. The significance stage is what scales with spectrum count; the index build is paid once.
Both get reported, separately labelled.

### 3.2 Baseline (task A1 — **measured**)

`cargo run -p msgf-genfunc --example profile --release` on `main` @ `52c3f62`, three runs, 1,406 F13
spectra, single thread, nominal grid, `HCD_HighRes_Tryp.param`, `-ti 0,1`:

| Stage | Time (mean of 3) | Share |
|---|--:|--:|
| preprocess + scored spectrum | 39 ms | 0.9% |
| per-spectrum node tables (`spec-tables`) | 844 ms | 18.8% |
| **graph build** | **558 ms** | **12.5%** |
| **generating-function DP** | **3,041 ms** | **67.9%** |
| merge + final tail | ~1 ms | 0.0% |
| **total** | **4,482 ms** | |

**Throughput: 312 / 311 / 319 spectra/s** (2.6% spread) — i.e. **314 spectra/s, statistically
identical to the 314 recorded at `1b14c5f`.** Nothing measurable has landed on `main` since the
trials began: `c8400e3` is PR #9's *first* optimization only, worth ~3% of the DP (the report
attributes 15.6 of PR #9's 19.1% to the sink-ancestor prune, which is still unmerged). Every win in
§2 is on a branch.

### 3.3 The Amdahl budget — why the DP alone cannot get us there

All scenarios below start from the measured 4,482 ms and assume the `spec-tables` branch lands
(844 → ~116 ms, the measured 7.4× on that stage), because it is finished, bit-exact, and cheap.
CPU-hours are single-thread, per 100,000 spectra; the Java column of reference is **0.3750
CPU-hr/100k** (MS-GF+ `TimeGenFunc`, re-measured on this machine 2026-07-25: 18,979 ms mean over
JIT-warm passes 1–4, 74 spectra/s):

| DP | Graph build | Total | Stage speedup | CPU-hr/100k | vs Java |
|---|---|--:|--:|--:|--:|
| — (`main` as-is, measured) | — | 4,482 ms | 1.00× | 0.0886 | 4.2× |
| 1× | 1× | 3,754 ms | 1.19× | 0.0742 | 5.1× |
| 1.31× — exact tail prune, F13 thresholds | 1× | 3,034 ms | 1.48× | 0.0599 | 6.3× |
| 2.68× — exact tail prune, real-ID corpus | 1× | 1,848 ms | 2.43× | 0.0365 | 10.3× |
| 2.68× | 3× | 1,476 ms | 3.04× | 0.0292 | 12.9× |
| 6× | 3× | 848 ms | **5.29×** | 0.0168 | 22.4× |
| 10× | 3× | 645 ms | 6.95× | 0.0127 | 29.4× |
| 20× | 3× | 493 ms | 9.09× | 0.0097 | 38.5× |
| ∞ | 1× | 713 ms | 6.29× | 0.0141 | 26.6× |
| ∞ | 3× | 342 ms | 13.1× | 0.0068 | 55.5× |

Rows below the first are **projections from separately measured components**, not end-to-end
measurements; the "real-ID corpus" rows use the DeNovoScore−20 prune figure as a stand-in until A2
exists. Multi-core does not scale linearly here — the measured 32-core throughput is 4,209
spectra/s, 13.4× the single-core rate, so quote CPU-hours (which are stable) and derive wall-clock
from the measured scaling factor, not from the core count.

Read off the two target lines: **5× needs DP ≈ 5.5× *and* graph build ≈ 3×. 10× needs DP ≈ 25×,**
which no amount of DP arithmetic delivers — it only comes from not running the DP at all on most
spectra.

**Conclusion, and it is the single most important line in this plan: the 5–10× target is not
reachable by p-value algorithms alone.** Graph construction must come down roughly 3× as well
(§5.6). It is in scope — it is per-spectrum null-model work, not candidate generation.

**Second conclusion, from `research-trials/dp-pruning-limits.md`: the DP columns above past ~6×
cannot come from pruning cells.** That trial measured a per-edge cost floor — at DeNovoScore − 5,
390× fewer cells bought only 6.17×, with ~16% of the unpruned DP left in bound and range passes that
are independent of score-support width. So a 10–20× DP must come from one of three other places:
**(a) not running the DP at all** for most spectra (a sound rejection gate, §5.2/§8.1 — whole-graph
elimination is not bounded by the per-edge floor), **(b) making an edge visit cheaper** (§5.4, §5.6),
or **(c) replacing the DP with an estimator** (§7). Cell pruning is P1/P2 table stakes, not the path
to the target.

On the high-precision grid both the DP and the graph build scale with node count (~274×), so the
split should roughly hold; that is untested and is task A6.

### 3.4 Gate cost is the budget, not gate quality

The cascade's cost ratio is approximately `C_gate/C_exact + f_fallback`. Our measured saddlepoint
gate costs **1/3.2 ≈ 31%** of an exact DP. Even at zero fallback that caps the DP contribution at
3.2×. So:

- a gate that is to buy 5–10× must cost **≤10%** of an exact DP;
- the streaming empirical null (§8.1, essentially free) and the exact early-exit corridor (§5.2,
  ~1% of DP for the bound sweep plus a truncated forward pass) are therefore **higher-priority
  gates** than the saddlepoint, despite the saddlepoint being the more impressive result;
- the saddlepoint's role is the **mid-tier estimator** for spectra that survive the cheap gates,
  not the first gate. Making it cheaper (§7.2, fused-θ sweeps) changes that ranking.

### 3.5 The benchmark corpus problem (blocking)

**F13 is the worst case and must not be the only benchmark.** It identifies essentially nothing
(`PLAN2.md` §4: MS-GF+ itself reports q = 1 for 4132/4133 rows), so its DeNovoScore−RawScore gap has
median 69 — the far end of the pruning table. Every gate in this plan gets *cheaper* as matches get
*better*, which is exactly what F13 does not contain. Measuring the cascade only on F13 understates
it, and measuring FDR power on F13 is meaningless.

**And F13 is inverted for gates versus prunes — neither corpus alone can measure the cascade.** The
threshold prune gets *cheaper* as matches improve, so F13 is its worst case (1.31×). A rejection gate
fires exactly when the best PSM is non-significant, so F13 — where essentially nothing identifies —
is its **best** case, and would flatter it just as badly in the opposite direction. A cascade
benchmarked only on F13 reports a near-perfect gate hit rate that no real run will see; benchmarked
only on a high-quality corpus it reports a gate that never fires. Every cascade number must be
reported on both, with the gate hit rate stated.

Task A2 therefore builds a second benchmark from a held-out MassIVE-KB shard (annotated `SEQ=`
peptides plus mass-identical shuffled decoys — the corpus
`validation/eval_trained_model.py library` already knows how to construct), giving real
identifications, real target-decoy labels, and a usable 1% FDR gate. Benchmark set, per the brief:
narrow tryptic; variable-mod; broad/nonspecific if supported; a null-heavy set (F13 serves); a
threshold-enriched set; and exhaustively enumerated small graphs for exactness.

---

## 4. Workstream A — harness and instrumentation (P0, blocks everything)

| Task | Deliverable |
|---|---|
| A1 | **Done** — see §3.2. Baseline on `main` @ `52c3f62`: 314 spectra/s, 4,482 ms, DP 67.9% / node tables 18.8% / graph build 12.5%. Still to add to `profile`: latency percentiles and edge-visit counts. |
| A2 | The real-identification benchmark corpus (§3.5) + a manifest recording provenance and licence for each set. |
| A3 | A per-spectrum diagnostics record (§10.2) emitted by every mode: tier that produced the value, exactness class, bound/interval, sample count, stopping reason, fallback flag, seed. |
| A4 | Exhaustive small-graph oracle: enumerate every path of tiny synthetic graphs, compare against the DP cell-by-cell with `f64::to_bits`. Extends the existing synthetic-graph tests on `worktree-genfunc-algo-speedups`. Graphs must be **peptide-shaped** (traps §6) — realistic mass-step-to-extent ratios, or path probabilities underflow `f64` and correct code fails the test. |
| A5 | `cargo run -p msgf-genfunc --example pvaluebench` — one binary, mode selected at runtime (not a cargo feature), emitting the comparison table of §12 for every mode over every corpus. |
| A6 | The same profile on the high-precision grid, to check that the §3.3 budget holds where the work actually matters. |

Runtime mode selection, not compile-time: every experiment must be runnable against one binary and
one corpus so that ablations are comparable.

---

## 5. Workstream B — exact dynamic programming

Everything here must stay **bit-exact**: same multiply-adds, same order, `f64::to_bits`-identical
against the unpruned DP. Pruning is allowed only where it provably removes cells that contribute
nothing to the requested tail.

### 5.1 B1 — consolidate the measured exact work (first, and highest value/effort)

Land what already exists and is already validated: draft PR #9 (sink-ancestor prune + redundant
range work), the `spec-tables` node-score cache and ion-major sweep, and threshold-aware
`compute_tail_into` with the `search_at_charge` reorder that scores candidates before building the
generating function.

Two known interactions to resolve while merging: PR #9's `max_remaining == i32::MIN` test *is* the
sink-ancestor prune, so the gains are **not additive**; and the pruned path raises per-node score
floors, so predecessor slices no longer fit whole and overlap clipping returns on that path while
the unpruned path keeps the fast one.

Also finish the piece the trial report leaves open: **`msgf rescore` still uses the unpruned
`compute`.** It caches one generating function per `(scan, charge)` across several PSMs, so its
threshold is the *minimum* RawScore over that key's PSMs. Wire it up (`msgf-cli/src/rescore.rs`).

### 5.2 B2 — two-sided ambiguity corridor and exact early exit

Extend the reverse sweep so each node `v` carries `U_v` (max remaining score to any valid sink,
implemented), `L_v` (min remaining score), and `Z_v` (total probability mass of valid suffix
completions). Then, for a forward state with partial score `s` and mass `q`:

```text
s + U_v <  T   -> discard          (success impossible)
s + L_v >= T   -> answer += q*Z_v  (all completions succeed)
otherwise      -> propagate        (ambiguous)
```

Only the corridor `T - U_v <= s < T - L_v` needs downstream DP. Exact, and it complements the
implemented `U_v` prune directly.

This is also what unlocks **exact early-exit rejection** (the highest-value missing piece): with
`Z_v` accumulating guaranteed-success mass, the running total is an exact *lower bound* on
`p(T_best)`. The moment it exceeds `p_cut`, the best match for that spectrum is proven
non-significant and the whole spectrum is abandoned — no approximation, no false rejections
possible. Combined with the monotonic spectrum-level shortcut (start from the best retained PSM;
if it fails, every lower-scoring PSM for that spectrum fails too), this is the cheapest sound gate
available.

Tests: `U/L/Z` verified against exhaustive enumeration on small graphs (A4); equality with the
baseline p-value; and a census of guaranteed-fail / guaranteed-succeed / ambiguous states by graph
depth, so we can see *where* the corridor closes.

### 5.3 B3 — sharing across the isotope-error candidates (scope correction)

The brief proposes one graph with each permitted precursor mass as a sink. **That does not apply
here as written.** In `build_reverse_graph`, `node_score[m] = round(prefix[complement_mass − m] +
suffix[m])` — the node scores are a function of the candidate mass, which is why
`Graph::recompute_node_scores` exists. Two isotope-error candidates share edges and node masses but
have *different node scores everywhere*, so a single multi-sink DP would not reproduce
`merge_group`, and MS-GF+ sums the candidates independently anyway (no double-counting risk, and no
sharing opportunity in the DP itself).

What can genuinely be shared, and is worth measuring:

- the reverse reachability mask and the `max_remaining` / `L_v` / `Z_v` sweeps — same shape, adjacent
  sinks, currently rebuilt per candidate. Measure the sweep cost first; it is small next to
  convolution;
- `Graph::recompute_node_scores` currently zeroes and rewrites the whole array per candidate;
- the per-node `sinks.contains(&i)` test in the DP's hot loop (`compute_into`) — tiny, but it is per
  node, per candidate.

Record the negative result explicitly so the multi-sink idea is not proposed again.

### 5.4 B4 — score-distribution representation

- **Fuse the DP's two per-edge passes** — the direct attack on the 2.9 ns/edge floor, and the single
  most concrete unexplored exact idea (`dp-pruning-limits.md` §9.1). The score-range pass and the
  convolution pass each gather a `NodeDist` per incoming edge and each re-derive
  `score_diff = node_score + edge_score`. Caching ≤32 resolved descriptors
  `(src_ptr, len, min_score, score_diff, prob)` on the stack during the first pass removes the second
  gather entirely. Untested. Note §7 of that report: the *live* edge visit is the cost, so anything
  that adds a pass to remove work (the forward `max_achievable` sweep) measured 7–15% **slower**.
- **O(1) tails.** `ScoreDist::spectral_probability` re-sums the tail on every query. One cumulative
  suffix-sum array after the DP makes every PSM lookup O(1) — irrelevant for one PSM per spectrum,
  material for `rescore` and for top-N search output. Must sum in an order that reproduces the
  current per-query sum bit-for-bit, or be gated behind an explicitly-named API.
- **Sparse/dense hybrid.** Measure early-node occupancy *before* implementing; switch to dense above
  a measured threshold; conversion must not change addition order.
- **`edge_prob` is effectively constant** (one value per amino acid) yet stored per edge — ~236 KB
  per graph written at build time and streamed on every DP edge. Replacing it with a per-AA lookup
  removes a memory stream from both stages.
- Predecessor aggregation, tight active score bounds after each update, and wider contiguous
  vectorized shifts, subject to the no-reassociation rule.
- **`f32` is out.** It is a fidelity change, not a representation change.

### 5.5 B5 — bidirectional / meet-in-the-middle GF (stretch)

Split near the mass midpoint, compute prefix distributions from the source and suffix tails from the
sinks, and join across a separator using cumulative suffix tails so the join is linear in score
width. The obstacle specific to this codebase is context: edge scores depend on resolved node masses
and the cleavage term rides on the source edge, so a separator must carry enough context to
reproduce them. Do this only after B1–B4 are understood; compare several cut masses; combine with
the corridor on both halves.

### 5.6 B6 — graph build (required for the 5× target, §3.3)

At 14.2% and ~682 MB of allocation per run it is now the second stage. Targets: the `edge_prob`
removal above; building edge arrays once per spectrum rather than per candidate (partly done);
avoiding the full-array rewrite in `recompute_node_scores`; and checking whether the two counting
passes can be collapsed given that the per-node edge count is a pure function of `m` and the amino
acid set.

---

## 6. Workstream C — certified coarse-score bounds

Quantize score contributions to width `b`, propagating both the floor-rounded and ceil-rounded
contribution to obtain `p_lower <= p_exact <= p_upper`, refined along `b = 16 → 8 → 4 → 2 → 1`.
Reject when `p_lower > p_cut`; accept as promising when `p_upper` is below the cutoff; refine or go
exact only when the interval straddles a decision boundary.

**This is not the coarsening already rejected.** That one binned scores and reported the result as
if it were exact (≈0.15–0.2 log10 of error per unit of threshold misplacement, outside tolerance at
`b = 2`). The certified version reports an interval that provably contains the exact value; the same
error becomes interval width, which is sound. Whether it is *useful* depends on how fast the
interval closes — measure accumulated per-residue rounding, and test checkpoint quantization as an
alternative to rounding every transition. Kill it early if the interval does not close.

---

## 7. Workstream D — approximate estimators

All of these live behind explicitly named APIs. Nothing in `search`, `rescore`, or the FDR path may
select them silently.

### 7.1 D1 — productionize the saddlepoint estimator

`msgf_genfunc::saddle` already exists and is measured: 3.2× the exact DP, median log10 ratio −0.002,
96.1% within 0.05, 100% within 0.30, with a documented error-vs-tail-depth table. Remaining work is
not the mathematics: a stable public API, the error assertion as a test (per the fidelity rule that
an approximate estimator must assert its measured error), a guard band derived from that table, an
uncertainty estimate per spectrum, and automatic exact fallback outside the validated domain
(accuracy decays as the threshold approaches the maximum achievable score).

### 7.2 D2 — fewer saddlepoint sweeps

~2.5 tilted sweeps per graph is near the Newton floor. Fusing three or four θ values into one AVX
sweep and recovering `K'`/`K''` by finite differences could approach one sweep-equivalent, roughly
doubling the 3.2× — which by §3.4 is what would move the saddlepoint from mid-tier estimator into
gate territory.

### 7.3 D3 — sequential conditional-null sampling (reject-only)

Backward scalar partition function `Z_v`, then sample null paths conditioned on reaching an allowed
sink, and estimate `P(S >= T_best)` as a Bernoulli rate over 16 → 32 → 64 → … paths. Reject when an
anytime-valid lower confidence bound exceeds `p_cut`; otherwise fall through. Deliberately
reject-only: strong or boundary matches always get exact treatment. Record samples, exceedances,
interval, stopping reason, and seed.

Note the synergy: `Z_v` is the same quantity B2 needs. Build it once.

### 7.4 D4 — exponentially tilted importance sampling

Same tail, tilted path measure, `p(T) = M(λ)·E_λ[e^{−λS}·1{S ≥ T}]`. λ chosen so the tilted expected
score sits near `T` — the saddlepoint solver already computes exactly that λ, so D1 and D4 share
their expensive part. Return estimate, standard error, effective sample size, and interval; fall
back to exact on poor ESS or when the interval straddles a threshold. Test tilt mixtures for
multimodal high-score path families.

### 7.5 D5 — adaptive multilevel splitting (stretch, likely skip)

Only if tilted sampling demonstrably fails on multimodal cases. Requires a progress coordinate that
includes remaining-score potential, not raw partial score. Report estimator variance, resampling
bias, and seed stability.

---

## 8. Workstream E — amortized preclassification

### 8.1 E1 — streaming empirical null (cheapest gate; high priority per §3.4)

During the normal database scan, maintain compact score histograms or quantile sketches over
null-like candidates (decoys especially), and use
`p_emp = (1 + #{S_null >= T}) / (1 + N_null)` as a near-free rejection signal. Stratify by
mass/length, charge, modification class, and cleavage class where sample size allows; use it only
within its supported resolution (deep tails fall through to a real method); treat it as a
*database-conditional* statistic, never as the reported SpecEValue. Validate target-decoy
exchangeability, and keep winner selection out of the null calibration to avoid leakage.

Interaction with `PLAN2.md`: decoys are already produced by the same pass, so the sketch costs a
counter update per scored candidate.

### 8.2 E2 — distilled conditional survival model (research track, off by default)

Train a **monotonic-in-score** model on exact-GF outputs to predict `−log10 p` or a conservative
interval, from features that are cheap at scoring time (threshold, precursor mass, charge, graph
size, score moments and bounds, peak count, spectrum entropy, search configuration, and diagnostics
from one scalar DP pass). Split by spectrum/run/dataset; optimize conservative decision error near
the operational cutoffs rather than mean regression loss; produce intervals by conformal
calibration or quantile models; reject only on the conservative side; go exact inside the band or
outside the validated domain.

Two repo-specific cautions. This introduces a trained artifact — `docs/models.md` and `LICENSING.md`
rules apply: it must be trained from a corpus we may ship, or fetched rather than vendored. And it
must never become the default: this project's value proposition is that the number is exact.

---

## 9. Workstream F — alternative confidence statistics (parallel track)

Evaluated on **identification power and calibration at fixed FDR**, never on agreement with the
exact p-value: distribution of the best null score (GEV / generalized Pareto / conditional quantile
models); a spectrum-peptide likelihood ratio; and a PEP / q-value mixture model over score, score
gap, precursor error, matched-ion evidence, peptide properties, and spectrum quality.

**Evaluation rule.** Correlating with the exact p-value proves nothing. An alternative must be
calibrated on held-out runs and must be non-inferior in PSM/peptide yield at 1% FDR on *every*
primary benchmark, including modification-heavy modes. And per `CLAUDE.md`: judge on ground truth
(held-out MassIVE-KB shard with `SEQ=` annotations and mass-identical decoys), not on agreement with
MS-GF+, whose own F13 top hits are 50% decoy.

If this track wins, it does not replace the exact path — it becomes a fourth mode (§10.1).

---

## 10. The cascade

### 10.1 Modes

```text
T_best (best retained PSM for the spectrum)
  │
  ├─ exact early-exit lower-bound DP (§5.2)      ── lower bound > p_cut ─► reject spectrum
  ├─ streaming empirical null / distilled model  ── safely non-significant ─► reject
  ├─ saddlepoint or tilted IS (§7)               ── interval resolves category ─► return
  └─ exact optimized GF (§5): corridor + shared bounds + hybrid storage
```

Four user-selectable modes, exact being the default:

| Mode | Contract |
|---|---|
| `exact` (default) | always the exact spectral p-value; bit-exact to MS-GF+ |
| `reject-only` | shortcuts may only *prove non-significance*; everything else exact |
| `approx` | approximate values with confidence metadata + exact fallback near boundaries |
| `fast` | alternative calibrated statistics permitted, once §9 validates them |

Selected at runtime (`--pvalue-mode`), so one binary serves the benchmark harness.

### 10.2 Result provenance (non-negotiable)

Every value carries which layer produced it and what kind of number it is:

```rust
pub enum Exactness {
    Exact,                                   // full distribution
    ExactAboveThreshold { valid_from: i32 }, // pruned DP; tail below valid_from not computed
    Bounded { lo: f64, hi: f64 },            // certified interval (§6)
    Estimated { se: f64, ess: f64 },         // sampled (§7.3, §7.4)
    Approximated { guard_log10: f64 },       // saddlepoint (§7.1)
    Predicted { .. },                        // model (§8.2)
    RejectedBelowCutoff { bound: f64 },      // gate proved non-significance; no value computed
}
```

The existing `GenFunc::valid_from` (on `worktree-genfunc-algo-speedups`) is the first member of this
enum and the model for the rest: querying below the validity threshold is a programming error, not a
silently wrong answer. The search/rescore TSV gains a provenance column; the diagnostics record
(A3) gets the full struct.

---

## 11. Why shortcuts are safe: the FDR-invariance argument

This is the correctness argument the whole cascade rests on, and it must be *tested*, not asserted.

1. **Rejection only moves rows that are already below the cutoff.** A spectrum rejected by a sound
   gate has `p(T_best) > p_cut`; its PSMs sort below every reported identification.
2. **q-values above the cutoff count only rows above the cutoff.** MS-GF+'s TDA (`PLAN2.md` §1.4)
   accumulates target and decoy counts down a SpecEValue-sorted list, so perturbing the *order* of
   rows below the cutoff cannot change the q-value of any row above it.
3. **Therefore the 1% FDR identification list is invariant** under any gate that is sound at `p_cut`,
   provided `p_cut` is chosen safely below the operating point.
4. **But the reported values for rejected rows do change.** `msgf search` writes top-N matches for
   every spectrum, including non-significant ones. In `reject-only` mode those rows carry
   `RejectedBelowCutoff { bound }` rather than a number. That is a deliberate, documented divergence
   from MS-GF+ output — the default `exact` mode must not have it.

**Gate G-FDR (blocking for any mode above `exact`):** on every benchmark corpus, `exact` and the
cascade must produce (a) an identical identification list at 1% FDR, (b) identical q-values for
every row with q below the reporting threshold, and (c) zero exact-significant PSMs rejected, plus a
stated statistical upper bound on that rate for the sampling-based tiers.

---

## 12. Acceptance criteria

### Correctness

- Exact methods: `f64::to_bits`-identical to the baseline DP on exhaustive small graphs, on synthetic
  graphs across thresholds/seeds/cleavage settings, and on every benchmark spectrum; the F13 search
  output TSV stays byte-identical; `golden_specprob` and the `#[ignore]`d `golden_search`
  (1161/1161 exact RawScore + DeNovoScore) unchanged.
- Bounded methods: contain the exact value at the promised coverage on every tested spectrum.
- Sampled methods: unbiasedness (or quantified finite-sample bias) on tractable cases; reproducible
  under a fixed seed; stable across seeds at the stated sample size.
- Every shortcut decision reconstructible from the stored diagnostics.

### Performance and power

| Path | Bar |
|---|---|
| Exact | ≥ **2×** the A1 baseline on the aggregate benchmark, numerically equivalent |
| Hybrid cascade | **5–10×** aggregate on the per-spectrum stage (§3.1 metric) |
| `reject-only` | zero observed false rejections + an explicit upper bound on the rate |
| `approx` | error reported per **relative** tail-depth decade (`tail / Z(0)`, traps §5 — absolute SpecEValue is not a depth scale); exact fallback retained near unstable regions |
| Alternative statistics | non-inferior PSM/peptide yield at 1% FDR on *every* primary benchmark |

### Required ablations

Baseline; representation improvements only; corridor only; shared-bounds only; full exact stack;
each gate alone; each gate + exact fallback; full cascade. Each measured on each corpus, with
fallback fraction and mean gate cost as a fraction of exact-GF cost reported alongside speed.

---

## 13. Milestones

| # | Milestone | Contents | Gate |
|---|---|---|---|
| **P0** | Harness | A1–A5 | baseline reproducible; A4 oracle green |
| **P1** | Exact consolidation | B1 (PR #9 + spec-tables + threshold prune + rescore threshold) | exact ≥ 2× baseline; goldens bit-exact |
| **P2** | Exact corridor + early exit | B2, B3 | corridor exact on A4; early-exit gate sound; G-FDR |
| **P3** | Constants | B4, B6 | graph build ≈3×; allocation down; still bit-exact |
| **P4** | Cheap gates | E1, cascade skeleton + provenance (§10) | gate cost ≤10% of exact DP; G-FDR |
| **P5** | Estimators | D1, D2, then C; D3/D4 if the budget still needs them | published error tables; fallback triggers documented |
| **P6** | Alternatives | E2, F | non-inferior at 1% FDR, or explicitly shelved with the measurement |

Milestones are ordered by (measured value)/(effort), not by the brief's order: P1 is
already-validated code, P2 is the only *sound* large gate, P4 is the cheapest gate, and the
estimators — despite being the most interesting — come after, because §3.4 says a 31%-of-exact gate
cannot carry the target on its own.

---

## 14. Do not redo these

Recorded with evidence in `research-trials/generating-function-optimization.md` (synthesis) and
`research-trials/dp-pruning-limits.md` (derivations + full instrument output for the aggressive
prunes):

- **Top-cap truncation** — invalid: partial path scores are not monotone (signed node/edge scores +
  cleavage penalty). At cap +30: mean |log10| error 0.50, worst 27.2.
- **Chernoff / tilted low-end trimming** — removed 3% more cells at ε = 1e-3 while the sweep itself
  cost ~13% of the DP. Net 0.92×–1.03×.
- **Extra dead-node sweep** — measured 7–15% *slower*.
- **FFT convolution** — does not apply: each edge is a shift, not a general convolution, and the mass
  axis is not translation-invariant because node scores differ per node.
- **Uncertified score-lattice coarsening** — outside the 0.05 tolerance by `b = 2`. (§6 is a
  different, certified construction.)
- **Contracting score-inert regions** — the amino-acid kernel is already sparser than its contraction
  (~186 entry-node terms per exit node vs 21 taps).
- **Mass-axis reachability pruning on the nominal grid** — everything at/above mass 57 is reachable.
  May be worth revisiting at the low-mass end of the high-precision grid.
- **Shared multi-sink DP across isotope candidates** — inapplicable; node scores are
  candidate-dependent (§5.3).

---

## 15. Deliverables

- [ ] Baseline profile + benchmark manifest (A1, A2).
- [ ] Exhaustive small-graph oracle (A4).
- [ ] One branch per method, one runtime flag per method, one binary (A5).
- [ ] Unified result schema with exactness provenance (§10.2), surfaced in CLI output.
- [ ] Benchmark tables: runtime, latency percentiles, fallback fraction, error by tail decade, FDR
      power, with the full ablation grid (§12).
- [ ] G-FDR invariance test in CI-runnable form (§11).
- [ ] Default thresholds, guard bands, and documented failure domains with automatic fallback.
- [ ] A trial report per method under `research-trials/`, including the ones that fail.
- [ ] `PERFORMANCE.md` and `ALGORITHMIDEAS.md` updated with what landed.

---

## 16. Rules this plan does not get to break

- The **default path stays bit-exact** to MS-GF+ (`CLAUDE.md`): integer scores exact, SpecEValue
  within `|log10(rust/java)| ≤ 0.05`. No reassociated summation, no FMA contraction, no `f32`, no
  probability cutoffs on the exact path.
- An approximate estimator lives behind an **explicitly named API**, documents its measured error
  against the exact DP, and **asserts that error in a test**. Search, rescore, and FDR code never
  select one silently.
- Goldens are regenerated deliberately, never as a side effect
  (`validation/reference/build_all_golden.sh --with-java`); new goldens are gitignored and wired into
  that script or their tests skip forever.
- Nothing here touches the clean-room boundary (`LICENSING.md`) — but §8.2's model artifact would
  introduce a new trained file, which does.
