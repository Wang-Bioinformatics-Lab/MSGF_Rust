# PLAN4 — Desktop UI over a loopback web server

Execution plan for a **desktop-feeling graphical front-end** to the `msgf` binary: a web server bound
to `127.0.0.1` that serves an embedded UI and drives the existing search / decoy / FDR paths.

**Status: design doc, not started.**

The governing constraint is that this must not change what MSGF_Rust *is*. The workspace has exactly
**one** runtime dependency (`rayon`) across ten crates; `release.yml` already ships a single static
`msgf` binary for five targets. The UI is required to preserve both properties: **zero new
dependencies, no JavaScript build step, still one file to distribute.**

| Decision | Choice | Why |
|---|---|---|
| HTTP stack | hand-rolled on `std::net::TcpListener` | zero new deps; loopback needs no TLS, no HTTP/2, no async |
| Frontend | vanilla HTML/CSS/JS, `include_str!`-embedded | `cargo build` remains the entire build; `release.yml` untouched |
| Scope (v1) | run a search + browse results | covers `search`, `decoy`, `fdr`; no spectrum viewer, no project management |
| Packaging | `msgf ui` subcommand, not a second binary | one artifact to distribute, one version number |

---

## 1. Goal and scope

**Primary success criterion:** a user who has never seen a terminal can double-click the binary, pick
an MGF and a FASTA, press Run, watch progress, and read a sorted PSM table with an FDR summary —
and the TSV they export is **byte-identical** to what `msgf search` would have written for the same
parameters.

That last clause is the fidelity contract restated for the UI: the UI is a *front-end*, not a second
implementation. It must not have its own scoring, its own FDR, or its own output formatting.

**In scope**

- `msgf ui` subcommand: bind loopback, open the browser, serve the app.
- Server-side file browsing so multi-GB MGF/FASTA files are never uploaded through the browser.
- A search form covering the flags in `msgf search --help` that matter to a human (§6.1).
- One running job at a time, with live progress, a phase indicator, and cancellation.
- A paged, sortable, filterable results table; FDR summary; TSV export.
- Decoy-FASTA construction (`msgf decoy`) and q-value annotation (`msgf fdr`) as secondary screens.

**Out of scope for v1** — each is a plausible PLAN4.1, none is needed for the criterion above.

Spectrum/PSM viewer with annotated b/y peaks. Multi-run project management, run history, or
model/parameter comparison. A `rescore` screen. A `msgf-train` screen. Anything multi-user, remote,
or authenticated beyond the loopback token. mzML input (`msgf-io` is MGF-only; that is PLAN-level
work, not UI work). Live-updating results while the search is still running — results appear when
the job completes.

**Non-goals.** Serving on a non-loopback interface. Bundling a browser engine (Tauri/Electron/webview
— all of them break "pure Rust distribution, one file"). Reimplementing any scoring or FDR logic in
JavaScript.

---

## 2. Architecture

A new crate `msgf-ui` holds the server as a **library**; `msgf-cli` gains a `ui` subcommand that
calls it. The UI is a default-on feature so `cargo build -p msgf-cli --no-default-features` still
yields the lean CLI.

```
                    browser (127.0.0.1:PORT)
                            │  HTTP/1.1 + SSE
   ┌────────────────────────┴───────────────────────────┐
   │ msgf-ui                                            │
   │  http/    TcpListener, request parse, response,    │
   │           SSE writer, fixed thread pool            │
   │  guard/   Host + Origin + token checks             │
   │  assets/  index.html, app.js, style.css (embedded) │
   │  api/     JSON endpoints, hand-rolled encoder      │
   │  jobs/    background thread, progress, cancel      │
   │  fs/      rooted directory listing                 │
   └────────────────────────┬───────────────────────────┘
                            │  plain function calls
        msgf-db · msgf-search · msgf-fdr · msgf-scorer · msgf-io
                    (unchanged, except §5's progress hook)
```

New crate `rust/crates/msgf-ui`, dependencies: the `msgf-*` crates it drives, plus `rayon` (already
in the tree) for the search itself. **Nothing else.**

---

## 3. The HTTP layer

What has to be written, and — equally important — what deliberately is not.

**Supported:** HTTP/1.1 `GET` and `POST`; request line + headers; `Content-Length` request bodies;
`Connection: keep-alive`; responses with `Content-Length`; `text/event-stream` for progress;
`Content-Type` by extension for the three static assets.

**Not supported, and rejected explicitly with a clear status code:** TLS, HTTP/2, chunked *request*
encoding (411), `multipart/form-data` (uploads are designed out, §4), range requests, compression,
cookies, redirects, virtual hosts.

A hand-rolled server is only safe if it is bounded. These are requirements, not niceties:

| Bound | Value | Failure mode it prevents |
|---|---|---|
| Max request line + headers | 16 KiB | unbounded header read → OOM |
| Max request body | 1 MiB (JSON only) | same; uploads are not a route |
| Read/write timeout | 30 s | a half-open socket parking a pool thread forever |
| Worker threads | fixed pool, 8 | connection-per-thread exhaustion |
| Concurrent SSE streams | 4 | pool starvation by long-lived streams |

The SSE stream is the one long-lived response: it writes `data: {...}\n\n` per progress tick and a
`: keepalive` comment every 15 s so an idle proxy or a sleeping laptop does not silently drop it.
Because SSE holds a pool thread for the life of the job, the pool must be larger than the SSE cap —
hence 8 vs 4.

---

## 4. Security

A server on `127.0.0.1` is still a network service, and any web page the user visits can attempt to
talk to it. Three real attacks apply; all three get an explicit defense.

1. **DNS rebinding.** A remote page resolves its own hostname to `127.0.0.1` and then scripts our
   API from a foreign origin. Defense: reject any request whose `Host` header is not exactly
   `127.0.0.1:<port>` or `localhost:<port>`. This is the single most important check in the server.
2. **Cross-origin drive-by.** Any page can issue a simple `POST` without a preflight. Defense:
   require `Origin` to be absent or match our own on every state-changing route, and require the
   session token (below). Send no permissive CORS headers — ever.
3. **Any other local process reading the user's disk.** The filesystem-listing route is the sharpest
   edge in this design. Defenses: a **session token** generated per launch (see below), embedded in
   the URL that gets opened, kept in `sessionStorage`, and required on every `/api/*` route; plus
   root-scoping and traversal rejection (§5.3).

**Token generation without a dependency.** There is no `rand` in the tree and there will not be.
Use the OS: read 32 bytes from `/dev/urandom` on Unix and `BCryptGenRandom` (or
`RtlGenRandom`) on Windows, hex-encode. If the source cannot be read, **fail to start** rather than
falling back to a timestamp — a predictable token is worse than no UI.

Also: bind `127.0.0.1` explicitly (never `0.0.0.0`), default to an ephemeral port (`:0`) and print
the real one, send `Content-Security-Policy: default-src 'self'` with no inline scripts, and document
plainly in `--help` and the README that this is a **single-user local tool** and the port must not be
forwarded or exposed.

---

## 5. Changes required outside `msgf-ui`

Three upstream changes. All are additive; none may alter any number the CLI produces.

### 5.1 A progress hook in `msgf-search`

`SearchEngine::run` is today a blocking `par_iter` with no visibility
(`msgf-search/src/search.rs:238`). Add, without
changing the existing signature:

```rust
pub fn run_with(&self, spectra: &[Spectrum], obs: &(dyn Observer + Sync)) -> Vec<Psm>;
pub trait Observer { fn spectrum_done(&self, n: usize); fn cancelled(&self) -> bool { false } }
```

`run` becomes `run_with(spectra, &())`, with a no-op impl for `()`. The closure increments an
`AtomicUsize` and checks an `AtomicBool`. This is a pure-plumbing change — the golden search tests
must stay green untouched, and that is the gate.

`PLAN5` wants the same hook for per-chunk logging, so build it once here.

### 5.2 Honest phases, not a fake percentage

A search is: read FASTA → build peptide index → search spectra → assign q-values → write. Only the
third has a natural per-item denominator. Index build is a large, silent, non-decomposable chunk
(`msgf-search/src/index.rs:46`, a `par_iter` over proteins — it *could* report per-protein progress,
and should).
The UI reports a **phase name plus a bar only where a denominator exists**, and a spinner elsewhere.
Do not invent a percentage for the index phase.

### 5.3 A rooted, traversal-proof path resolver

`GET /api/fs?path=…` lists a directory. Requirements: canonicalize the requested path, verify it is
inside one of the session's allowed roots (the launch CWD plus any root the user explicitly adds),
reject symlinks that escape, never follow `..` textually — canonicalize first, compare after. Return
name, size, is_dir, modified; filter files by extension for the picker (`.mgf`, `.fasta`/`.fa`,
`.param`, `.tsv`). Hidden files off by default.

---

## 6. The API

All routes require the token; all return JSON with a hand-rolled encoder (§7).

| Method | Route | Purpose |
|---|---|---|
| GET | `/` `/app.js` `/style.css` | embedded assets |
| GET | `/api/hello` | version, bundled-model identity, core count |
| GET | `/api/fs?path=` | directory listing (§5.3) |
| POST | `/api/search` | start a job → `{job_id}` |
| GET | `/api/jobs/{id}` | status snapshot (poll fallback for SSE) |
| GET | `/api/jobs/{id}/events` | SSE progress stream |
| POST | `/api/jobs/{id}/cancel` | request cancellation |
| GET | `/api/jobs/{id}/psms?offset&limit&sort&dir&q&max_q` | paged results |
| GET | `/api/jobs/{id}/tsv` | export — streams `report::write_tsv` |
| POST | `/api/decoy` | `msgf decoy` equivalent |
| POST | `/api/fdr` | `msgf fdr` equivalent |

### 6.1 What the search form exposes

Mirror `msgf search --help`, grouped as the help text already groups it: **required** (spectra,
FASTA), **search space** (enzyme, ntt, missed cleavages, min/max length, precursor tolerance,
isotope-error range, matches per spectrum, charge range), **modifications** (fixed/variable mod
rows with a `C+57.021464` free-text form plus a short menu of common mods, max mods per peptide),
**target-decoy** (`--tda` toggle, decoy prefix), **other** (model: bundled or a `.param` path;
threads).

Two form behaviors are load-bearing rather than cosmetic:

- **Warn where the CLI warns.** No decoys in the database → the CLI prints that q-values are not an
  FDR estimate. The UI must say the same thing, in the results header, not swallow it.
- **Surface the model identity.** `crate::model::load` distinguishes bundled from `--param`. Show
  which one ran, because the bundled model is deliberately *not* bit-exact to MS-GF+ (CLAUDE.md) and
  a user comparing against Java output needs to see that at a glance.

### 6.2 Results paging

A search over a human FASTA can produce >100k rows; sending them all as JSON would stall the tab.
Sorting and filtering are **server-side** over the in-memory `Vec<Psm>`, returning pages of ~200.
Filters: substring on peptide/protein, `max_q`, target-only toggle.

The summary line comes from `report::summary` — the same string the CLI prints to stderr — so the UI
cannot drift from the CLI's account of the same run.

---

## 7. JSON without serde

`serde_json` is a dev-dependency only, and adding it as a runtime dep would break the posture. Write
a ~120-line encoder in `msgf-ui::api::json`. Two details are where hand-rolled JSON usually breaks,
and both occur in real data here:

- **String escaping.** Protein accessions come from arbitrary FASTA headers and routinely contain
  `"`, `\`, and non-ASCII. Escape `"` `\` and every control byte as `\u00XX`; emit valid UTF-8
  through. Fuzz this function against a corpus of real FASTA headers.
- **Non-finite floats.** `q_value` is `f32::NAN` until `assign_q_values` runs (`report.rs` renders it
  as `NA`), and `NaN`/`Infinity` are **not valid JSON**. Emit `null` and render it in the UI as
  `NA` — matching the CLI's spelling.

---

## 8. Frontend

`crates/msgf-ui/assets/`: `index.html`, `app.js` (~600 lines, no framework), `style.css`. Embedded
with `include_str!`.

Three screens behind a tab strip: **Search** (form → progress → results), **Decoy** (FASTA in →
revCat FASTA out), **FDR** (PSM TSV in → annotated TSV out). Light/dark via
`prefers-color-scheme`. No inline scripts, so the CSP in §4 holds.

**Dev-loop affordance:** when `MSGF_UI_ASSET_DIR` is set, read the assets from disk instead of the
embedded copies. Otherwise every CSS tweak costs a `cargo build`, and the UI will be developed
badly. Embedded remains the default and the shipped path.

---

## 9. Testing

Integration tests bind `127.0.0.1:0` and speak HTTP over a real `std::net::TcpStream` — no new test
dependencies. Follow CLAUDE.md's skip pattern for anything that needs `validation/data/`.

| Test | Asserts |
|---|---|
| `routing` | every route's status, content type, 404 for unknown |
| `host_header_rejected` | `Host: evil.com` → 403 (the §4.1 defense) |
| `token_required` | `/api/*` without the token → 401 |
| `fs_traversal_rejected` | `../../etc/passwd` and an escaping symlink → 403 |
| `body_limit` | 2 MiB body → 413, connection not wedged |
| `sse_framing` | well-formed `data:` frames, keepalive comments |
| `json_escaping` | quotes/backslashes/control chars/UTF-8; `NaN` → `null` |
| `cancel_stops_a_job` | job ends, and is reported **cancelled**, not complete |
| **`ui_tsv_matches_cli_tsv`** | **the export is byte-identical to `msgf search` for the same params** |

The last one is the plan's real gate; the others are hygiene. A cancelled job must never present
partial results as a finished search — the `par_iter` yields whatever it had, and labelling that
"done" would be a quietly wrong FDR.

---

## 10. Packaging

`release.yml` needs **no change** — still `cargo build --release -p msgf-cli`, still five targets,
still one file. Add `ui` to the subcommand table in `main.rs`'s `USAGE`.

Browser launch: `xdg-open` (Linux), `open` (macOS), `cmd /c start` (Windows); `--no-open` to
suppress. Always print the URL, so a headless or SSH session still works via a port-forward the user
sets up deliberately.

`--port <N>` for a fixed port; default ephemeral. `--open-timeout`/retry is not needed — print and
move on.

---

## 11. Milestones

| # | Deliverable | Gate |
|---|---|---|
| U1 | `msgf-ui` crate; server, router, embedded assets, `msgf ui` launches and serves a page | `routing` green; zero new deps in `Cargo.lock` |
| U2 | Security layer: token, Host/Origin checks, limits, timeouts | all §9 security tests green |
| U3 | `/api/fs` + file-picker UI | traversal tests green; can pick a real MGF/FASTA |
| U4 | Search job: `run_with` hook, background thread, SSE progress, cancel | golden search tests still green **untouched**; `cancel_stops_a_job` |
| U5 | Results: paged table, filters, summary, TSV export | **`ui_tsv_matches_cli_tsv`** |
| U6 | Decoy + FDR screens; README + `--help` | end-to-end run on F13 from a clean checkout |

U4's gate is the one to watch: if adding the observer changes any golden output, the hook is wrong.

---

## 12. Risks

**The hand-rolled server is the main risk.** It is ~400 lines that must be correct about header
parsing, `Content-Length`, keep-alive framing, and timeouts. Mitigation: the bounds table in §3 is
mandatory, the surface is deliberately tiny, and every limit has a test. If the server turns out to
need chunked requests, WebSockets, or uploads, that is the signal to revisit the `tiny-http` option
rather than to keep growing this.

**Memory.** Results live in RAM for the life of the job; a large search holds `Vec<Psm>` plus the
peptide index at once. Show the row count and let the user export and start a new run; do not add a
database.

**Scope drift toward a workbench.** The spectrum viewer is the obvious next ask and is genuinely
useful — but it needs peak-annotation plumbing the CLI has never needed. Keep it out of v1 so U5's
gate stays reachable.

---

## 13. Open questions

1. **Concurrent jobs.** v1 runs one at a time (the search already saturates every core). If users
   want to queue runs, that is a job-queue design, not a tweak.
2. **`msgf ui` vs a double-clickable artifact.** A `.app`/`.exe` shim that just runs `msgf ui` would
   complete the "desktop" story. Cheap, but it does touch `release.yml` — deferred until U6 lands.
3. **Where the file-browser roots start.** CWD is the obvious default; whether to add `$HOME` by
   default is a real security-vs-convenience call worth making explicitly.
