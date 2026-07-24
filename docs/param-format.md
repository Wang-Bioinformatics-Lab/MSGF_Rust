# The MS-GF+ `.param` scoring-model binary format

Reverse-engineered from the reference reader/writer and validated by round-trip: this repo's
`msgf_scorer::read_param` decodes it and `msgf_scorer::write_param` re-encodes all four high-res UC
models **byte-for-byte identical** (`tests/roundtrip_write.rs`). This document is the normative spec
for anyone emitting a `.param` (e.g. a trainer) — the format is an interface; see `docs/models.md`
§2 for the licensing boundary between the *format* (documented here, unencumbered) and UC's trained
*bytes* (test-only, gitignored).

## Conventions

- **Endianness:** big-endian, everywhere (Java `DataOutputStream`).
- **Primitives:** `i32` = 4 bytes; `f32` = 4 bytes IEEE-754; `bool` = 1 byte (`0`/`1`); `u8` = 1 byte.
- **String** (`jstring`): 1 byte unsigned length `L` (count of UTF-16 code units), then `L × 2` bytes
  of UTF-16BE. (Java `writeByte(len)` + `writeChars`.) Model identity strings are ASCII, so `L` =
  character count. Max length 255.
- **Optional string:** a `jstring` with `L = 0` (a single `0x00` byte) means **absent** (`None`).
- Everything is a flat positional stream — there are no tags, offsets, or lengths beyond the counts
  called out below. A misread desyncs the whole tail; the trailing sentinel is the only integrity
  check.

## Layout (in stream order)

| # | Section | Fields |
|---|---|---|
| 1 | **Header / identity** | `version:i32`; `activation:jstring`; `instrument:jstring`; `enzyme:jstring?`; `protocol:jstring?`; `mme_is_ppm:bool`; `mme_value:f32`; `apply_deconvolution:bool`; `deconvolution_error_tolerance:f32` |
| 2 | **Charge histogram** | `n:i32`, then `n × (charge:i32, count:i32)` |
| 3 | **Partitions** | `n:i32`; `num_segments:i32`; then `n × (charge:i32, parent_mass:f32, seg:i32)` |
| 4 | **Precursor offsets** | `n:i32`, then `n × (charge:i32, reduced_charge:i32, offset:f32, tol_is_ppm:bool, tol_value:f32, frequency:f32)` |
| 5 | **Fragment offsets** | **one block per partition** (count = §3 `n`, no separate length). Each block: `size:i32`, then `size × (is_prefix:bool, charge:i32, offset:f32, frequency:f32)` |
| 6 | **Rank distributions** | `max_rank:i32`. Let `C = max_rank + 1`. For **each partition whose §5 block is non-empty**, in partition order: for each fragment ion in that block (§5 order), `C × f32`; then one trailing `C × f32` **noise** row. Empty-block partitions contribute nothing. |
| 7 | **Error / isotope dists** | `error_scaling_factor:i32` (`E`). If `E > 0`, let `W = 2·E + 1`; then for **every** partition (§3 `n` of them): `signal: W × f32`, `noise: W × f32`, `ion_existence: 4 × f32`. If `E == 0`, this section is just the one int. |
| 8 | **Sentinel** | `i32` = `0x7FFFFFFF` (`Integer.MAX_VALUE`). Reader asserts this; a mismatch means the parse desynced upstream. |

### Ordering constraints (must hold for a reader to decode)

- **Partitions are canonically sorted** by `(charge, seg, parent_mass)` ascending, unique. MS-GF+
  writes them from a `TreeSet`, so files are already in this order; the reader re-sorts defensively.
  §5/§6/§7 are all **parallel to this sorted partition list** — emit them in the same order.
- **§6 skips empty partitions** but **§7 does not** — §6 iterates only partitions with ≥1 fragment
  ion, §7 iterates all partitions. Getting this wrong desyncs the sentinel.
- Within a §6 partition, ion rows are in the **same order as the §5 block**, with **noise last**.

## Two read-side transforms (not stored in the bytes)

A faithful re-encode must account for these, or byte-identity/round-trip breaks:

1. **Fragment ion name is derived, not stored.** The reader synthesizes
   `name = "{P|S}_{charge}_{round(offset)}"` where `P`=prefix, `S`=suffix and `round` is Java
   `Math.round` = `floor(x + 0.5)`. The writer must **not** emit the name (it isn't in the format);
   the reader will reproduce it on the next read.
2. **Zero ion-existence is floored.** The reader replaces any `ion_existence` entry equal to `0.0`
   with `0.001`. So `write(read(f))` is byte-identical to `f` **iff** `f` has no exact-zero
   ion-existence entries (all four shipped high-res models satisfy this); otherwise the re-encode is
   still *struct*-identical (`read(write(m)) == m`) but a few bytes differ at those slots.

## Where the trained numbers live

For a trainer, "producing a model" = filling these and serialising §1–§8:

- **§5 `frequency` / `offset`** — which ion types are worth scoring per partition, and how often each
  is observed. *(Fragment offset frequencies.)*
- **§6 ion rows + noise row** — the load-bearing scores. The per-node score is
  `ln( ionFreq[rank] / (noiseFreq[rank] · min(ionCharge, num_segments)) )`
  (`ScoringModel::score_from_table`); the `max_rank` column is the "ion absent" bin. *(Rank
  distributions.)*
- **§7 signal/noise/ion_existence** — high-res mass-error term used by edge scoring.
- **§4** — precursor offset frequencies (used in preprocessing/filtering).
- **§1–§3** identity + partition scheme + tolerances the above are conditioned on.

## Worked size sanity

`HCD_QExactive_Tryp.param` = 741,431 bytes, 140 partitions, `max_rank = 150` (⇒ `C = 151`),
`error_scaling_factor = 100` (⇒ `W = 201`). §6 alone is ~125,632 `f32` (`rank_dist_total_floats`);
§7 is `140 × (201 + 201 + 4) × 4` bytes ≈ 227 KB. These match the committed model-inventory golden
(`validation/golden/models/*.model.golden.json`).
