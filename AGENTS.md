# Repository Guidelines

## Project Structure & Module Organization

The Rust workspace lives in `rust/`. Its crates separate chemistry, spectrum I/O, scoring, generating functions, database handling, FDR, search, training, the `msgf` facade, and the `msgf-cli` binary under `rust/crates/`. Keep library code in each crate's `src/`, integration tests in `tests/`, and Criterion benchmarks in `benches/`.

`validation/` contains Python, shell, and Java tooling for reference-data retrieval and golden generation. Technical specifications are in `docs/`; consult `docs/param-format.md` before changing model serialization and `docs/training.md` before changing training counts. See `LICENSING.md` before adding data or models.

## Build, Test, and Development Commands

Run Rust commands from `rust/`:

- `cargo build --release -p msgf-cli` builds `target/release/msgf`.
- `cargo run -p msgf-cli -- --help` runs the CLI locally.
- `cargo test --workspace` runs unit and golden integration tests.
- `cargo test -p msgf-genfunc --test golden_specprob` runs one integration test target.
- `cargo clippy --workspace --all-targets` checks common correctness and style issues.
- `cargo fmt --all --check` verifies formatting; use `cargo fmt --all` to apply it.
- `cargo bench -p msgf-genfunc --bench genfunc` runs the generating-function benchmark.

From `validation/`, run `python3 regression/run_regression.py` to re-derive and verify available fixtures.

## Coding Style & Naming Conventions

Use Rust 2021 conventions and rustfmt defaults (four-space indentation). Name modules, functions, and test files in `snake_case`; types and traits in `UpperCamelCase`; constants in `SCREAMING_SNAKE_CASE`. Keep crate responsibilities narrow and public APIs documented. In fidelity-sensitive scoring code, matching Java arithmetic and operation order takes precedence over stylistic refactoring.

## Testing Guidelines

Add focused `#[test]` unit tests beside implementation code and cross-crate or golden checks under `crates/<crate>/tests/`, using descriptive `snake_case` names. No percentage coverage target is defined; every behavior change should include regression coverage. Golden tests must skip cleanly when unvendored reference data is absent. Regenerate goldens only deliberately through `validation/reference/`; never commit `validation/data/`, MS-GF+ models/JARs, or UC-derived outputs.

## Commit & Pull Request Guidelines

History follows scoped Conventional Commit subjects such as `fix(search): report decoy status` and `perf(genfunc): reduce edge work`. Keep commits focused and use `feat`, `fix`, `perf`, `docs`, or `chore` with a relevant scope.

Pull requests should explain the behavior change, affected crates, validation commands run, and any numerical or performance impact. Link related issues and call out fixture, licensing, or clean-room implications. Include CLI output examples when user-facing behavior changes.
