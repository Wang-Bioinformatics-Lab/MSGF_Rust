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
| `msgf-scorer` | ✅ `.param` reader + scoring math | `validation/golden/models/` — all 4 models decoded & sentinel-aligned; `node_score`/`missing_ion_score` match MS-GF+ on 17,644 values (worst Δ 4e-7). Spectrum preprocessing + RawScore next |
| `msgf-genfunc` | ⬜ next | DeNovoScore + SpecEValue from `iprg2013_F13.golden.json` (the generating-function core) |
| `msgf-cli` | ⬜ later | end-to-end differential test vs the F13 golden |
| `msgf-search` | ⬜ later | Sage-inspired search engine |

## Notes

- Our code is MIT. The reference MS-GF+ models/data are UC-licensed and are **not** vendored
  (kept in the gitignored `../validation/data/`). Golden JSON is derived numeric facts.
- Tests that read `../validation/data/` (gitignored) **skip** gracefully when the data is absent,
  so a fresh clone still passes; run `../validation/fetch_reference_data.sh` to enable them.
- `msgf-chem::scaling` already encodes the nominal (`0.999497`) and high-precision (`274.335215`)
  mass grids — the discretization the generating-function DP will run on.
