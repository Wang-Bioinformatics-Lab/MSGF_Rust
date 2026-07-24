# rust/ — MSGF_Rust workspace

Cargo workspace for the Rust reimplementation of MS-GF+ significance scoring. Every crate is
validated against the frozen golden corpus in `../validation/` (see `../PLAN.md`).

```bash
cargo test --workspace     # unit tests + golden validation against ../validation/golden/
cargo clippy --workspace --all-targets
cargo fmt --all
```

## Crates

| Crate | Status | Validated by |
|---|---|---|
| `msgf-chem` | ✅ implemented | `validation/golden/chemistry/` — atomic/residue/peptide masses, b/y ions, tolerance, mass scaling |
| `msgf-io` | ✅ MGF reader | `validation/golden/spectra/` — byte-for-byte peak-list hashes over F13 (1,406) + test.mgf (5,760) |
| `msgf-scorer` | ✅ full RawScore (node + edge) | `.param` decoded; node scores match on 95,306 values; preprocessing (incl. deconvolution) exact; **full per-peptide RawScore (node + edge, `DBScanScorer`) matches MS-GF+ 30/30 across charge 2/3 + mods**. Generating function (SpecEValue) next |
| `msgf-genfunc` | ✅ **SpecEValue p-value (bit-exact)** | Generating-function DP (ScoreDist + GeneratingFunctionGroup) over the de novo graph. **DeNovoScore + SpecEValue match MS-GF+ 30/30; score distributions agree to ~2e-8** (`golden_specprob.rs`) |
| `msgf-cli` | ✅ `msgf rescore` (RawScore/DeNovoScore/SpecEValue) | `tests/golden_rescore.rs` — the `msgf` binary reproduces MS-GF+ **30/30** on F13 (RawScore + DeNovoScore exact, SpecEValue to f64 noise) |
| `msgf-search` | ⬜ later | Sage-inspired search engine |

## Benchmarks

`cargo bench -p msgf-scorer` (criterion; needs `validation/data/`). Current single-threaded
numbers on the F13 high-res set (1,406 spectra, HCD/QExactive model), measuring the per-spectrum
work a search does once — preprocess + build scored spectrum + compute `prefixScore`/`suffixScore`:

| Benchmark | Time |
|---|---|
| `read_param_model` (one-time model load) | ~0.78 ms |
| `preprocess_one` (per spectrum) | ~6.3 µs |
| `score_one_spectrum` (preprocess + node scores) | ~0.81 ms |
| `preprocess_and_score_all` (1,406 spectra) | **~1.16 s → ~1,210 spectra/s** |

Reference: MS-GF+ (Java, JIT-warm) does the identical work at **~770 spectra/s** on this machine,
so the Rust port is ~1.5× faster **single-threaded** — before any SIMD or multi-core (the
per-spectrum work is embarrassingly parallel). Caching the per-segment partition lookup (instead of
re-running it per nominal mass) was a +160% throughput win. The generating function — the real
high-res hot path — is where the larger gains will come.

## Notes

- Our code is MIT. The reference MS-GF+ models/data are UC-licensed and are **not** vendored
  (kept in the gitignored `../validation/data/`). Golden JSON is derived numeric facts.
- Tests that read `../validation/data/` (gitignored) **skip** gracefully when the data is absent,
  so a fresh clone still passes; run `../validation/fetch_reference_data.sh` to enable them.
- `msgf-chem::scaling` already encodes the nominal (`0.999497`) and high-precision (`274.335215`)
  mass grids — the discretization the generating-function DP will run on.
