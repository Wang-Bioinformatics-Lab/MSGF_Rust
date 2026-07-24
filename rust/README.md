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
| `msgf-scorer` | ✅ full RawScore (node + edge) | `.param` decoded; node scores match on 95,306 values; preprocessing (incl. deconvolution) exact; **full per-peptide RawScore (node + edge, `DBScanScorer`) matches MS-GF+ 30/30 across charge 2/3 + mods** |
| `msgf-genfunc` | ✅ **SpecEValue p-value (bit-exact)** | Generating-function DP (ScoreDist + GeneratingFunctionGroup) over the de novo graph. **DeNovoScore + SpecEValue match MS-GF+ 30/30; score distributions agree to ~2e-8** (`golden_specprob.rs`) |
| `msgf-db` | ✅ FASTA, decoys, digestion | **Target-decoy FASTA byte-identical to MS-GF+'s `.revCat.fasta`** for both references (`golden_decoy_fasta.rs`, PLAN2 TD-1) |
| `msgf-fdr` | ✅ PSM + peptide q-values | **`QValue`/`PepQValue` reproduce MS-GF+ exactly for all 1,610 F13 PSMs** (`golden_fdr.rs`, PLAN2 TD-2 Gate 1) and the `TargetDecoyAnalysis` map bit-for-bit over 14 synthetic cases (`golden_fdrmap.rs`, Gate 2). Both goldens are MS-GF+-derived, so both are regenerated locally and both tests skip without them |
| `msgf-search` | ✅ database search | End-to-end vs MS-GF+ on F13 (`golden_search.rs`): **0 scans where our best candidate scores lower**; on the 1,161 scans with the same top peptide, **RawScore and DeNovoScore exact 1,161/1,161** and SpecEValue within tolerance |
| `msgf` | ✅ umbrella facade | One dependency re-exporting the whole pipeline (`msgf::{chem, io, scorer, genfunc, db, fdr, search}`); `default-features = false` drops the search engine |
| `msgf-cli` | ✅ `search` / `rescore` / `decoy` / `fdr` | `msgf-cli/tests/golden_rescore.rs` — the binary reproduces MS-GF+ **30/30** on F13 |

## Using it as a library

```toml
[dependencies]
msgf = { git = "https://github.com/mwang87/MSGF_Rust" }
# scoring only, without the search engine and its rayon dependency:
msgf = { git = "https://github.com/mwang87/MSGF_Rust", default-features = false }
```

```rust
use msgf::prelude::*;
let model = msgf::scorer::read_param_file("HCD_HighRes_Tryp.param")?;
let db = msgf::db::fasta::ProteinDb::read("human.revCat.fasta", "XXX_")?;
let index = msgf::search::PeptideIndex::build(&db, &DigestParams::default(), &Default::default());
```

The individual `msgf-*` crates can also be depended on directly; the facade adds no wrappers.

## Command line

```bash
msgf search  -s run.mgf -p HCD_HighRes_Tryp.param -d human.revCat.fasta -o psms.tsv
msgf rescore -s run.mgf -p HCD_HighRes_Tryp.param -i psms.tsv
msgf decoy   -d human.fasta -o human.revCat.fasta
msgf fdr     -i psms.tsv -o psms.q.tsv
```

`msgf <command> --help` lists every flag. Two defaults worth knowing: **missed cleavages are
unlimited** (MS-GF+'s own default — pass `-c 2` for a conventional, far smaller index), and
q-values need decoys (`--tda`, or search a concatenated `*.revCat.fasta`).

### Heavy tests

The F13-vs-MS-GF+ search comparison builds a ~48M-candidate index (~2 GB) and is `#[ignore]`d:

```bash
cargo test -p msgf-search --release -- --ignored --nocapture
```

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
