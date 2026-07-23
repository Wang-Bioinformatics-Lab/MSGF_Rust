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
| `msgf-scorer` | ✅ `.param` reader + scoring + scored-spectrum | models decoded & sentinel-aligned; `node_score` matches MS-GF+ on 17,644 values; **`prefixScore`/`suffixScore` match on 95,306 values across 30 real high-res spectra (worst Δ 9e-7)**. Preprocessing (deconvolution) + edge scores + RawScore summation next |
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
