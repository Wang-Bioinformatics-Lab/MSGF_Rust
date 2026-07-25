# Plans

Design and execution plans. `PLAN.md` is the authoritative strategy document; the numbered plans are
self-contained execution plans for one workstream each and carry their own status line.

| Doc | Workstream | Status |
|---|---|---|
| [PLAN.md](PLAN.md) | Authoritative design doc: phases, decisions D1–D5, algorithm derivation | living |
| [PLAN1.md](PLAN1.md) | Model ownership — own the `.param`, train our own model | in progress |
| [PLAN2.md](PLAN2.md) | Target-decoy + FDR: decoy FASTA, q-values, `msgf-search` wiring | TD-1/2/3 implemented; §4 oracle problem open |
| [PLAN3.md](PLAN3.md) | Spectral p-value acceleration: 5–10× on the significance stage | design doc, not started |
| [PLAN4.md](PLAN4.md) | Desktop UI: an embedded, zero-dependency web server on `127.0.0.1` | design doc, not started |
| [PLAN5.md](PLAN5.md) | Nextflow scale-out: scatter spectra, gather FDR | design doc, not started |
| [PLAN6.md](PLAN6.md) | timsTOF (Bruker `.d`) DDA support: direct reader, fragment tolerance, TOF model | design doc, not started |

Related, and not plans:

- `../PERFORMANCE.md` — measured Rust-vs-Java timings for what already landed.
- `../ALGORITHMIDEAS.md` — index of algorithm/performance research; `../research-trials/` holds the
  detailed reports. PLAN3 is the plan; those are its evidence.
- `../docs/` — normative specs (`param-format.md`), model provenance (`models.md`), trainer
  semantics (`training.md`).
