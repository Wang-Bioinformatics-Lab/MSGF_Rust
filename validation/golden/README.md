# golden/ — frozen reference outputs

Two kinds of golden live here, and only one of them is committed.

## Committed (ours, no UC input)

| Family | What it is | Generator |
|---|---|---|
| `chemistry/` | atomic/residue/peptide masses, b/y ions — computed from published physical constants | `../reference/make_chemistry_golden.py` |
| `models/param_inventory.golden.json` | size, SHA-256 and identity header of each `.param` file — a **manifest**, no trained numbers | `../reference/make_param_inventory.py` |

These carry no MS-GF+ content, so they stay in git and their tests assert (rather than skip).

## Not committed (derived from MS-GF+ or its UC-licensed data)

Everything else: `iprg2013_F13.*`, `worked_example.golden.json`, `rawscore/`, `spectra/`,
`models/*.model.golden.json`, `models/node_scores.golden.json`, `fdr/**`.

These are the project's numeric oracle and they remain essential — they are simply **generated on
your machine instead of shipped**, exactly like `../data/`. The reason is licensing: they are
outputs of MS-GF+ (Copyright UC Regents, non-commercial) or embed fragments of its UC-licensed test
data (e.g. the spectra goldens store peak values from `F13.mgf`, the model goldens sample trained
table rows). MSGF_Rust ships MIT, so it distributes none of it. See `../../LICENSING.md`.

`UC_DERIVED.sha256` records what those files hashed to when they were removed from git
(MS-GF+ v2024.03.26) — provenance evidence, not a gate: a different MS-GF+ release may legitimately
produce different bytes.

## Regenerating

```bash
cd validation
./fetch_reference_data.sh --all                    # inputs + iPRG FASTA + the MS-GF+ jar (needs Java 11+)
cd reference && ./build_all_golden.sh --with-java  # every golden, then the regression suite
```

**`--with-java` is the flag that matters.** Without it only the no-Java families are built, and the
tests backed by everything else keep skipping. The MS-GF+-derived families are a dependency chain
(F13 search → scored spectrum → RawScore / SpecEValue / HighRes), so run the script rather than the
`make_*.sh` pieces individually; it sequences them correctly.

Without this step the golden-backed tests **skip** — `cargo test` on a clean checkout still passes,
it just validates less. `cargo test -- --nocapture | grep skip:` shows what was skipped, and a run
that reports zero skips is the one that actually exercised the fidelity contract.

The no-Java families (`chemistry/`, `spectra/`, `param_inventory`, `worked_example`) can be rebuilt
from `../reference/make_*.py` + `parse_msgf_tsv.py` alone; the rest need the jar. `fdr/` is the
cheapest of those — `../reference/make_fdr_golden.sh` needs only the jar and a JVM, no spectra,
models or database.
