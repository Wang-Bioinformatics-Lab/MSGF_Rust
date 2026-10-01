# MSGF_Rust speed scout (2026-10-01)

Branch `speed-scout` in `.worktrees/speed-scout`, from `main` `3073604`. Nothing was merged or pushed.
Scope: find, measure and rank speed-ups after the clean-room speed round (`b034177`), and prototype
the cheap ones. Every prototype marked exact passes the rebuilt byte-identity harness (section 6).

## 1. Summary

On a 16-core EPYC 9634, the branch tip with the optional `mimalloc` feature turned on gives these
times. All outputs are byte-identical to the `0fb0738` reference.

- **Full HeLa R01 search** (113,314 spectra, 16 threads): 21.9 s → 16.0 s (−27 %).
- **3,000-spectrum search, 1 thread**: 3.87 s → 3.09 s (−20 %).
- **Rescore of the HeLa subset, 1 thread**: 1.41 s → 1.13 s (−20 %).
- **Full-run rescore of 542k PSMs**: 103 s → 9.2 s (11×). This is mostly the new `--threads`.

The levers, best first:

1. **`rescore --threads`**: 8–10× on the full-run rescore. Committed on its own as `5d19794`. The
   user has already approved it.
2. **Generating function (GF):** the per-row descriptor build is now a vectorised pass, and the
   convolution uses explicit SIMD with 32-cell blocks. Rescore is 16–27 % faster.
3. **mimalloc, as an opt-in build feature.** Multi-threaded search is 8–26 % faster.
4. **Parallel MGF parsing and parallel TSV formatting.** A full-run search on 16 threads is 8–9 %
   faster.
5. **Allocation-free peptide keys** in search: 6 %.

The generating function is not limited by floating-point throughput, so the following gave no gain:

- FMA;
- AVX-512 over AVX2 (2–3 %);
- `target-cpu` builds;
- PGO, which made one path slower.

Making the DP f32 is not worth the risk of underflow.

## 2. Ranked opportunities

Machines:

- **G**: HPCC EPYC 9634 (Genoa, AVX-512), node r44 shared with other jobs, pinned cores.
- **M**: HPCC EPYC 7713 (Milan). The nodes were heavily loaded (load 100–340 on 128 threads).
- **R**: MRB Ryzen 9 9900X, quiet. Measured before the power loss.
- **W**: MRB worker, a KVM guest with no AVX.
- **O**: HPCC Opteron 6376, which has AVX but not AVX2, so it runs the SSE2 path. Loaded.

Each cell is the minimum wall time over 3–5 interleaved runs.

| # | Opportunity | Exact? | Measured effect | Effort | Risk | User decision? | Status |
|--:|---|:-:|---|:-:|:-:|:-:|---|
| 1 | **`msgf rescore --threads <N>`** (default: all cores), parallel over (scan, charge) groups | yes | Full-run rescore (542k PSMs): G 103 → 10.6 s (16 t); M 138 → 15.8 s (16 t); R 60.2 → 6.3 s (22 t). HeLa subset: R 0.83 → 0.10 s (22 t); W 2.91 → 0.28 s (44 t). | S | low | already approved | **`5d19794`**, own commit |
| 2 | **GF: vectorised descriptor pass + explicit SSE2/AVX2/AVX-512 convolution with 32-cell blocks** | yes | Rescore, 1 thread: G −16…−19 % (r01 1.41 → 1.15, r12 2.42 → 1.95); M −12…−17 %; R −22…−27 % (r01 0.83 → 0.62, r12 1.43 → 1.05). Full-run rescore, 1 thread, G: 103 → 85 s. Search, 1 thread, G: −6 %. | M | low on AVX2 and AVX-512; **O (SSE2 path) +2…+11 % slower** on rescore | no | **`ee679ec`** |
| 3 | **mimalloc global allocator** (opt-in cargo feature) | yes | Search 1 t: G −10 %, M −12 %. Search 16 t on the subset: G 0.64 → 0.49 s, M 0.96 → 0.71 s. Full search 16 t: G −8 %, M −10 %, and −8 % more on top of #4. Rescore −1…−2 %. RSS +35 MB. | S | low | **yes**: a new dependency that builds bundled C code (MIT) | **`12060d3`**, off by default |
| 4 | **Parallel MGF parsing + parallel TSV formatting** | yes | Full search 16 t: G 19.3 → 17.5 s. Full rescore 16 t: G 10.6 → 9.8 s. Subset search 16 t: 0.63 → 0.58 s. Neutral on 1 thread. **Full-rescore RSS 0.60 → 1.04 GB** (the whole MGF is held while it is parsed). | S | low | no (the RSS cost is worth a look) | **`bafd57c`** |
| 5 | **Search: allocation-free peptide keys, cached delta text, Fx hasher, reused map** | yes | Search 1 t: G 3.65 → 3.43 s (−6 %). Full search 22 t: R 11.38 → 10.68 s (−6 %). | S | low | no | **`f24c1e1`** |
| 6 | Per-candidate RawScore in search: share prefix sums and vertex lookups across candidates of the same window | yes | Not built. Measured share (HPCC perf, full run): `match_and_terminal_parts` 7–10 %, `edge` 2–3 %, `search_spectrum` 6 %. Estimate: 5–8 % of search. | M–L | medium | no | estimate |
| 7 | Less allocation in search (`proteins` Vec per hit, `Psm` strings) without mimalloc | yes | Not built. malloc/free is ~5 % of the full-run profile; mimalloc recovers most of it (#3). | M | low | no | estimate |
| 8 | GF: share work across PSMs, sinks and spectra; tighter pruning | yes | **Nothing left.** The GF is built once per (scan, charge) and the edge table is shared by every sink, both already. Vertex weights depend on the sink (V_p(p − m)), so the DP rows cannot be shared across sinks. Null distributions never repeat across spectra. The cut is already exact at the lowest reported RawScore. A higher cut changes reported rows, so it would only be possible behind an output-restricting option. | — | — | yes, for any output restriction | analysed |
| 9 | Index build | yes | 1 thread: 1.8–2.7 s of a 150 s full run (1–2 %). 22 threads: 0.39 s of 10.6 s (4 %). It is already rayon-parallel. Its 23 % share of the 3,000-spectrum benchmark comes from the small input. The `par_sort_unstable` order of equal masses decides grouping order, so it must not change. | M | medium | no | not worth it |
| 10 | AVX-512 kernel (vs the AVX2 kernel) | yes | R: 2–3 % (r01 0.62 vs 0.64 s, r12 1.05 vs 1.07 s). | — | low | no | inside #2 |
| 11 | FMA in the convolution | **no** | No gain: G r01 1.15 → 1.15 s; R ≤ 3 %. 14 of 542,461 rows differ on the full run at ti −1,2; ti 0,0 and ti 0,1 on the subset change md5. | S | — | (yes) | **rejected**, branch `speed-scout-fma` `a2cb57f` |
| 12 | f32 DP (8/16 lanes) | **no** | Not built. The kernel is not FLOP-bound (#10, #11), so the estimate is ≤ 15–20 % of GF time. The real full run has 4 PSMs with SpecEValue < 1e-38, the minimum being 1.5e-41, below f32's normal range, so per-row power-of-two scaling would be needed. | M–L | high | (yes) | not recommended |
| 13 | `-C target-cpu=x86-64-v3 / v4 / znver5` | yes | G and M: ±1 % on every workload once the kernels are explicit. Before #2, R gained 5 %. | S | low (portability) | (yes) | not recommended |
| 14 | PGO | yes | It can only be built through musl: zig's linker rejects `-u __llvm_profile_runtime`. M: neutral (±3 %). **G: rescore 1.7× slower**, because the training ran on the AVX2 dev VM and the AVX-512 path was laid out as cold code. The musl allocator makes 16-thread search 7× slower (131 vs 19 s). | M | high | (yes) | rejected |
| 15 | AVX-only tier for pre-AVX2 CPUs | yes | O: 0–6 % slower than the SSE2 path, because Piledriver splits 256-bit operations. | S | — | no | rejected (patch: `/data/researchdata-botnet/dia_gate/msgf_scout/patches/avx1_tier_rejected.patch`) |
| 16 | SSE2-path regression from #2 (generic KVM CPUs, pre-2013 x86) | yes | O: r01 5.23 → 5.79 s, r12 9.73 → 10.18 s, rcomp +10 %, s1 −4 %. W on PAD16 earlier: r12 6.40 → 5.40. Proposed fix: on the SSE2 tier only, keep the per-row branchy descriptor build and a 16-cell block. | S–M | low | no | open |
| 17 | aarch64 / NEON | yes | Not measured; there is no ARM box. The branch compiles for `aarch64-unknown-linux-gnu`. The portable kernel there is scalar Rust that the compiler may vectorise. An explicit `float64x2` kernel mirroring the SSE2 one is about 60 lines and exact. It needs a Graviton or Apple box to time. | S | low | no | estimate |

## 3. Profiles

### Phase cycles

These come from an instrumented build on R (rdtsc phase timers, feature `prof`, not committed).
The instrumented build is the PAD-16 version with the vectorised descriptor pass.

| Phase | rescore r01, before #2 | r01 with the new descriptor pass | search, 3,000 spectra, 1 t | full search, 1 t |
|---|--:|--:|--:|--:|
| Edge table | 5.9 % | 7.0 % | 3.1 % | 2.0 % |
| Integer passes (vertex, forward, backward, layout) | 9.1 % | 10.6 % | 4.8 % | 3.2 % |
| **Descriptor build** | **32.5 %** | merged into the next row | 19.1 % (desc + kernel) | 23.6 % (desc + kernel) |
| **Convolution kernel** | **41.2 %** | 69.6 % (desc + kernel) | — | — |
| Mixture + merge | 0.2 % | 0.2 % | 0.2 % | 0.2 % |
| Prepare + RawScore (caller) | 11.0 % | 12.6 % | 5.8 % (prepare) | 2.6 % |
| Search: window + key + RawScore | — | — | 32.3 % | 34.4 % |
| Search: sort/truncate hits | — | — | 6.2 % | 4.2 % |
| Search: PSM build | — | — | 1.0 % | 0.6 % |

Counters for r01 rescore:

- 2,960 GF builds over 5,920 sinks;
- 3.59 M non-empty rows;
- 55.4 M descriptors, about 15 edges per row.

Before #2, the descriptor build cost about 16 cycles per edge. The 16-cell block was compiled to
**scalar** `vmulsd`/`vaddsd` with stack spills, because LLVM did not SLP-vectorise the
`[f64; 16]` accumulator. That is why the earlier AVX-512 and `target-cpu` experiments did nothing,
and why explicit intrinsics gained.

Wall time by stage of the full-run search on R with 22 threads (instrumented):

| Stage | Time |
|---|--:|
| Index build | 0.39 s |
| MGF read + engine setup, serial | 0.93 s |
| Search | 12.2 s |
| TSV write, serial | 0.79 s |

The serial stages are about 20 % of an uninstrumented 22-thread run. That motivated #4.

### HPCC `perf` (EPYC, debug-symbol `b034177`)

| Workload | GF (`build_avx2`) | Other top symbols |
|---|--:|---|
| Rescore subset | 82 % | `site_table` 5.7 %, memmove 1.9 %, `prepare_with_cache` 1.8 %, MGF parse ~2 % |
| Rescore full run | 90 % | `site_table` 2.9 % |
| Search subset, `-c 2`, 1 t | 24 % | index-build closure 13.4 %, `peptide_string` 11.9 % (13.7 % inclusive), `rayon::sort::recurse` 8.9 %, RawScore 4.7 %, malloc/free ~4 % |
| Search full run, 4 t | 50 % | `peptide_string` 10.6 % (14.6 % inclusive, of which `format!` is 4 %), RawScore 7.2 %, `search_spectrum` 6.5 %, `edge` 2.6 %, SipHash 1.7 %, malloc/free ~5 % |
| Search full run, 1 t | 46 % | `peptide_string` 12.2 %, RawScore 9.8 %, `search_spectrum` 6.8 % |

## 4. Cumulative timings

### G: EPYC 9634, pinned (jobs 29331482, 29331750, 29332845, 29332004)

| Workload | `b034177` (main) | +GF #2 | tip `f24c1e1` | tip + mimalloc | + I/O #4 | + I/O + mimalloc |
|---|--:|--:|--:|--:|--:|--:|
| rescore ti 0,1 (1 t) | 1.41 | 1.15 | 1.15 | 1.13 | 1.17 | — |
| rescore ti 0,0 (1 t) | 0.88 | 0.74 | 0.74 | 0.73 | — | — |
| rescore ti −1,2 (1 t) | 2.42 | 1.96 | 1.95 | 1.93 | — | — |
| rescore comp + oxM (1 t) | 1.45 | 1.18 | 1.18 | 1.17 | — | — |
| search, 3k spectra, 1 t | 3.87 | 3.65 | 3.43 | 3.09 | 3.45 | — |
| search, 3k spectra, 16 t | 0.66 | — | 0.64 | 0.49 | 0.58 | 0.46 |
| search full R01, 16 t | 21.90 | — | 19.00 | 17.55 | 17.47 | **16.01** |
| rescore full (542k), 1 t | 103.2 | — | 85.1 | — | — | — |
| rescore full, 16 t | — | — | 10.57 | — | 9.81 | **9.22** |

### M: EPYC 7713 (loaded nodes; read the ratios, not the absolute times)

| Workload | main | tip | tip + mimalloc |
|---|--:|--:|--:|
| rescore ti 0,1 | 1.89 | 1.59 | 1.53 |
| rescore ti −1,2 | 3.29 | 2.72 | 2.65 |
| search, 1 t | 5.38 | 4.84 | 4.27 |
| search, 16 t | 0.98 | 0.96 | 0.71 |
| search full, 16 t | 32.5 | 27.9 | 25.1 |
| rescore full | 137.8 (1 t) | 115.4 (1 t) / 32.4 (4 t) / 15.8 (16 t) | — |

### R: Ryzen 9900X (quiet, one pinned core)

| Workload | main (`b034177`) | rescore threads only (`5d19794`) | + GF #2 | + keys #5 |
|---|--:|--:|--:|--:|
| rescore ti 0,1 | 0.83 | 0.83 | 0.62 | 0.65\* |
| rescore ti 0,0 | (0.54†) | 0.51 | 0.40 | — |
| rescore ti −1,2 | 1.45 | 1.43 | 1.05 | 1.07\* |
| rescore comp + oxM | (0.90†) | 0.85 | 0.63 | 0.64\* |
| search, 1 t | (2.34†) | 2.26 | 2.10 | 1.99 |
| search full, 22 t | 13.92‡ | 12.08 | 11.38 | 10.68 |

\* Measured after the reboot, with the host under other load.
† From the speed-round report (same machine, same build), not re-measured here; for single-thread runs `5d19794` runs the same code as main.
‡ A single, non-interleaved run.

Rescore thread scaling on R, HeLa subset ti 0,1:

| Threads | 1 | 2 | 4 | 8 | 12 | 22 (default) |
|---|--:|--:|--:|--:|--:|--:|
| Wall time | 0.84 s | 0.44 s | 0.24 s | 0.14 s | 0.11 s | 0.10 s |

RSS grows from 21 MB to 73 MB across that range.

## 5. Prototypes on `speed-scout`

| Commit | What | Byte-identity |
|---|---|---|
| `5d19794` | `rescore --threads` (+ unit test 1 vs 2/3/8 threads, + integration test 1/2/5/default threads) | Original harness 62/62 at default/1/3/7 threads; rebuilt harness 51/51 at default/1/7 threads (run on the tip, which contains it) |
| `ee679ec` | GF: vectorised descriptor pass, explicit SIMD kernels, 32-cell blocks, AVX-512 dispatch | Original harness 62/62; md5 identical on the AVX-512 (R, G), AVX2 and SSE2 (W, O) paths |
| `f24c1e1` | Search keys: no per-candidate allocation, cached delta text, Fx hasher | Original harness 62/62; rebuilt harness 51/51 ×3 |
| `3c74dee` | `MSGF_GF_MAX_ISA=avx2\|sse2` diagnostic, replacing `MSGF_NO_AVX512` | Rebuilt harness 51/51 with the cap unset and at sse2 |
| `bafd57c` | Parallel MGF parse (`read_mgf_file`) + parallel TSV formatting | Rebuilt harness 51/51 (cap unset and sse2); error messages identical; new parser-equivalence unit test |
| `12060d3` | Opt-in `mimalloc` feature (off by default) | Output md5 identical on every timed workload (G, M) |

These were also built and measured but not committed to `speed-scout`:

- **FMA**: branch `speed-scout-fma` (`a2cb57f`, REJECTED).
- **AVX-only tier**: patch only, in `/data/researchdata-botnet/dia_gate/msgf_scout/patches/`.
- **Phase-timer instrumentation**: not kept. Its patch was in the scratchpad that the reboot wiped; the measurements are in section 3.
- **`target-cpu` and PGO builds**: build flags only, nothing to commit.

`cargo test --workspace --release` passes 131 tests, with 0 failed and 1 ignored (the F13 human
search). The clean-room test-vector tests now skip, because the vector package was lost (see
section 6). `cargo fmt --check` is clean.

## 6. Validation, and what was lost

A reboot of the dev VM wiped the session scratchpad. That took with it:

- the clean-room test-vector package;
- the oracle binary;
- the original 62-comparison harness;
- the previous `MSGF_SPEED_REPORT.md`.

The harness was rebuilt:

- **Reference binary:** `0fb0738`, rebuilt from `git archive` with compiler output discarded and run
  as a black box only, as in the speed round. md5 `a664fcd5…`.
- **Harness location:** `/data/researchdata-botnet/dia_gate/msgf_scout/harness/compare.sh`. The
  reference outputs are cached in `harness/ref/`, and the oracle is in `oracle/`, also on HPCC
  under `…/msgf_scout/oracle/`.
- **What it covers (51 comparisons):**
  - the HeLa subset, 5 rescore configurations;
  - the full R01 run, 2 rescores of 542k PSMs;
  - 8 search configurations, including the full R01 search;
  - the MS-GF+ F13 models, 4 searches and 4 rescores;
  - decoy, and fdr on 3 inputs.
- **No longer covered:** the synthetic-spectra edge-case set and the two retrained alt-model rescores
  (both lost). The unit test in `5d19794` uses its own synthetic spectra.

Logs are in `harness/cmp_tip.log`, `cmp_isa_*.log` and `cmp_io_*.log`.

The synthetic MGF can be regenerated only by re-running `make_synthetic.py`, which was lost too.
To restore the full original harness, regenerate the test-vector package from its source.

## 7. Housekeeping and cleanup

These have to be cleaned up later. The MRB boxes lost power during this work.

- **mrb-pve5-lxc1:** `/tmp/scout/` (binaries, the 380 MB MGF, FASTA and outputs, about 1 GB).
- **mrb-botnet-worker-1:** `/tmp/scout/` (binaries and inputs).
- **NFS `/data/researchdata-botnet/dia_gate/msgf_scout/`:**
  - `bin/`, about 30 binaries, can go;
  - keep `in/`, `oracle/` and `harness/` for future byte-identity checks.
- **HPCC `/bigdata/mxwanglab/mingxunw/dia_gate/msgf_scout/`:**
  - `bin/`, the job logs and the sbatch scripts can go;
  - keep `in/` and `oracle/`.
  - The jobs removed their node-local `/tmp/msgf_scout_*` directories.

## 8. Decisions for the user

1. **Merge the exact prototypes:** `5d19794`, then `ee679ec`, `f24c1e1`, `3c74dee`, `bafd57c`. The
   first can be merged on its own.
   - For `bafd57c`, decide whether to accept the extra transient memory in `rescore` (about 0.4 GB
     on a 380 MB MGF), or to have rescore stream instead.
2. **Turn mimalloc on by default** (`12060d3`) for release builds. It adds a C dependency.
3. **SSE2-path regression** (#16): accept a ≤ 11 % slowdown on CPUs without AVX2, or keep the old
   descriptor build on that tier.
4. **Non-exact modes** (FMA, f32): measured or estimated as not worth offering.

## Decisions (user, 2026-10-01)

- Every exact commit was merged into `main`, including the extra ~0.4 GB of rescore memory from parallel MGF parsing.
- **`mimalloc` is now the default** for `msgf-cli`. To build without it, use `--no-default-features`. The library crates (`msgf-scorer`, `msgf-genfunc`, …) do not depend on it.
- The slowdown of up to 11 % on CPUs without AVX2 is accepted, since those CPUs are too old to matter.
- There is no inexact "fast" mode.
