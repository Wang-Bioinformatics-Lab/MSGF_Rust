# Spectral-probability scoring for DDA PSMs: functional specification

Version 1, 2026-09-30. Target module: `rust/src/dda/specprob/` in DIA_Proteomics_Rust (MIT/Apache).

This document says **what** the component computes, precisely enough that an independent
implementation reproduces the reference scorer's numbers. It says nothing about **how** the
reference program is organised, and you should not try to mirror it. Choose your own data
structures, loop order and names. Where the order of floating-point operations can change an
integer result, this document states the precision of each intermediate quantity and the
tolerance that is accepted.

The ground truth is the black-box oracle in `oracle_bin/`, run with the model in `model/`. The
test vectors in `test_vectors/` were produced by that oracle. See §11 for the acceptance criteria.

---

## 0. Scope and vocabulary

The component takes one centroided MS/MS spectrum, its precursor m/z and a charge, a trained
scoring model, and one or more candidate peptides. It returns three things:

| output | type | meaning |
|---|---|---|
| **RawScore** | integer | log-likelihood-ratio match score of one candidate against the spectrum, including the enzymatic-terminus terms |
| **DeNovoScore** | integer | the largest score *any* amino-acid string of the right nominal mass could obtain against this spectrum |
| **SpecEValue** | real in [0, 1] | spectral probability: the total probability, under an i.i.d. residue model, of the strings whose score is at least RawScore |

This is the MS-GF / MS-GF+ spectral-probability approach, from Kim, Gupta & Pevzner (JPR 2008),
Kim et al. (MCP 2010) and Kim & Pevzner (Nat Commun 2014):

- a spectrum becomes a vector of integer vertex scores over nominal masses;
- a high-resolution edge term scores the mass error between consecutive fragment peaks;
- the score distribution over all residue strings is computed by dynamic programming over an
  amino-acid graph whose edges are weighted by residue probabilities.

This document fixes the concrete choices of the reference implementation where the papers leave
them open. §12 lists which parts come from the papers and which were fixed from observed behaviour.

Out of scope:

- training models (the model is an input);
- database digestion and candidate generation;
- FDR;
- N-terminal-cleaving enzymes, discussed in §8.7.

**Notation.**

- `f32(·)` is rounding to IEEE-754 binary32 and `f64(·)` to binary64.
- ⊕ ⊖ ⊗ ⊘ are binary32 operations. Each rounds its result to f32, and they are evaluated left to right as written.
- Plain + − × / are exact or binary64, as stated in context.
- `ln` is the natural logarithm, evaluated in binary64.
- `S` is the number of mass segments in the model, `Rmax` the model's maximum rank, and `E` the model's error-scaling factor (§2).

---

## 1. Constants and rounding conventions

### 1.1 Physical constants (Da)

| name | value | use |
|---|---|---|
| m_H | 1.0078250319 | atomic masses (for residues and water) |
| m_C | 12.0 | |
| m_N | 14.0030740052 | |
| m_O | 15.9949146221 | |
| m_S | 31.9720707300 | |
| water W | 2·m_H + m_O (binary64) | peptide = Σ residues + W |
| proton p⁺ | 1.0072764669 | precursor charge removal (§3.1) |
| charge carrier c⁺ | 1.00727649 | precursor-peak suppression and isotope-cluster reduction (§3.3, §3.5); used as f32 |
| first isotope spacing Δ₁ | 13.00335483 − 12 = 1.00335483 | isotope-cluster reduction; used as f32 |
| second isotope spacing Δ₂ | 14.003241 − 13.00335483 = 0.99988617 | isotope-cluster reduction; used as f32 (this odd value is what the reference uses) |
| nominal scaler κ | 0.999497 (as f32) | real mass → nominal-mass grid |

**Residue masses** are computed in binary64 as elemental sums `c·m_C + h·m_H + n·m_N + o·m_O + s·m_S`
from the formulas below. The formula is [C, H, N, O, S] of the residue, which is the free amino acid minus H₂O.

```
G 2 3 1 1 0   A 3 5 1 1 0   S 3 5 1 2 0   P 5 7 1 1 0   V 5 9 1 1 0
T 4 7 1 2 0   C 3 5 1 1 1   L 6 11 1 1 0  I 6 11 1 1 0  N 4 6 2 2 0
D 4 5 1 3 0   Q 5 8 2 2 0   K 6 12 2 1 0  E 5 7 1 3 0   M 5 9 1 1 1
H 6 7 3 1 0   F 9 9 1 1 0   R 6 12 4 1 0  Y 9 9 1 2 0   W 11 10 2 1 0
```

Only these 20 letters are residues. Any other letter (B, J, O, U, X, Z, …) makes a candidate
unscorable.

### 1.2 Rounding functions (all three are used, deliberately)

- **R(x)**, *score rounding*: `floor(x ⊕ 0.5)`, computed on an f32 `x`. Halves go toward +∞, so
  R(−2.5) = −2. Applied to every score that turns from real into integer.
- **Rn(x)**, *nominal rounding*: round to nearest, **halves away from zero**, on an f32 value. Used
  only for mass → nominal conversion.
- **T(x)**, truncation toward zero, f32 → integer. Used for segment indices.

Real → integer conversions saturate. NaN converts to **0**, +∞ to i32::MAX and −∞ to i32::MIN.
The NaN rule is observable (§6.2, §10).

### 1.3 Nominal grid

- A real mass `m` maps to the nominal integer `ν(m) = Rn(f32(m) ⊗ κ)`.
- A nominal integer `k` maps back to the representative real mass `μ(k) = f32(k) ⊘ κ`.

Nominal masses of peptides are always **sums of per-residue nominal masses**, never ν of the summed
mass.

---

## 2. The scoring model (`.param` file) as an interface

The binary layout is specified in `model/param-format_from_MSGF_Rust.md`. The layout is big-endian
throughout. It holds a header, a charge histogram, partitions, precursor offsets, fragment ion types
per partition, rank tables, error tables and a sentinel. `tools/param_reader.py` is an independent
reader, written from that document, that you may use as a cross-check. The layout is not repeated
here. This section says only what a decoded model must provide and how its tables are indexed.

A decoded model **M** provides the following.

**Header scalars**

- **fragment tolerance τ**: a value and a unit, Da or ppm. The Da-equivalent at a theoretical m/z x is
  - `tol(x) = f32(value)` for Da, or
  - `tol(x) = f32(x·value·1e-6)` in binary64 for ppm.
- **ρ-scale**: the raw tolerance *number* `f32(value)`, whatever its unit (used in §4.4).
- **reduce isotope clusters** flag, plus its tolerance `τ_d` (f32).
- **S** (segments per precursor), **Rmax** and **E**.

**Partitions** `π_0 … π_{P−1}`. Each is a triple (charge `z_π`, lower parent-mass boundary `B_π`,
segment `s_π`).

- Order them by the key `(z_π, s_π, B_π)`, ascending, and remove duplicates.
- All per-partition tables below are parallel to this sorted list.

**Ion types of partition π.** An ordered list. Each entry is `(prefix?: bool, charge c ≥ 1, offset o: f32, frequency φ: f32)`.

- The theoretical m/z of a type for a real prefix or suffix residue mass `x` is `x ⊘ f32(c) ⊕ o`.
- The inverse mapping, from peak m/z `y` back to residue mass, is `(y ⊖ o) ⊗ f32(c)`.
- The type's **name** is `P` or `S`, then `_c_`, then `floor(o + 0.5)` as an integer (e.g. `S_1_19`).
  The name is not stored in the file. It is the join key for the rank tables. Two types of one
  partition never share a name in valid models. If they do, both use the row of the first.

**Rank tables.** These exist only for partitions whose ion-type list is non-empty.

- There is one row of `Rmax+1` f32 values per ion type, plus one **noise** row.
- For a row `row`, bin `b ∈ [0, Rmax−1]` stands for peak rank `b+1`. Bin `Rmax` stands for "no peak".

**Error tables.** These exist only if `E > 0`, and there is one per partition:

- `signal` and `noise`, each with `2E+1` f32 values;
- `existence`, with 4 f32 values. **Replace every existence value equal to 0.0 by 0.001 on read.**

**Precursor offsets.** A list of `(charge z_o, reduced charge r_o, offset δ_o, …, frequency)`.
Only `z_o`, `r_o` and `δ_o` are used. The per-entry tolerance fields and the frequency are
**ignored**.

**Bundled model.** `model/MSGFRust_HCD_HighRes_Tryp_v1.param`, SHA-256 `5f6ab76f…2d02`, has:

- τ = 0.5 Da, reduce = true, τ_d = 0.02, S = 2, Rmax = 150, E = 100;
- 236 partitions over charges 2–6;
- every (charge, segment) group has a partition with B = 0.0;
- precursor offsets for charges 2–6.

It is the only model the acceptance tests use.

---

## 3. Spectrum preparation

**Input:** peaks `(mz_i, I_i)` as read from the source, each converted to f32; the precursor m/z
`q` (binary64); and a charge `z ≥ 1`.

### 3.1 Masses

- **Parent (neutral precursor) mass**, in f32:
  `M = (f32(q) ⊗ f32(z)) ⊖ (f32(z) ⊗ f32(p⁺))`.
- **Candidate nominal mass**: `N₀ = Rn((M ⊖ f32(W)) ⊗ κ)`.
- If `N₀ < 50` or `N₀ > 10000`, the spectrum is not scored at this charge: there are no outputs.
  All later sections assume `50 ≤ N₀ ≤ 10000`.

### 3.2 Peak order

The working order of the peaks is m/z ascending. If the input is already non-decreasing in m/z,
keep the input order unchanged, including the order among equal m/z values. Otherwise, sort
stably by (m/z, then intensity) ascending. **No peak is ever removed** anywhere in §3, including
zero-intensity peaks, duplicate peaks and m/z 0.

### 3.3 Precursor-peak suppression

1. Choose the precursor-offset charge key. The trained charges are the distinct values `z_o` in
   the table. Take the largest trained charge that is ≤ z. If there is none, take the smallest
   trained charge. If the table is empty, skip this step.
2. For every table entry with that key, in table order:
   - Let `c = z − r_o`. Skip the entry if `c = 0`.
   - Compute the centre `x = ((M ⊕ (f32(c) ⊗ f32(c⁺))) ⊘ f32(c)) ⊕ δ_o`.
   - Compute `t = tol(x)`, using the model's fragment tolerance.
   - Set the intensity of every peak with `x ⊖ t ≤ mz ≤ x ⊕ t` (inclusive at both ends) to **0**.

### 3.4 Intensity ranks

Give each peak a rank from 1 to n, a permutation over all n peaks, zeros included:

- higher intensity gets a smaller rank;
- for equal intensity, the **higher m/z** gets the smaller rank;
- for full ties, the peak earlier in the working order gets the smaller rank.

Ranks are computed on the intensities from §3.3, **before** §3.5. They then stay attached to their
peaks unchanged.

### 3.5 Isotope-cluster charge reduction ("deconvolution")

This step runs only if the model's reduce flag is set. It rewrites m/z values in place. It never
changes intensities or ranks, and never removes peaks. Let `a[0..n)` be the m/z values in working
order. Some entries get rewritten during the pass, and **every comparison reads the current,
possibly already rewritten, values**. Keep a per-peak "consumed" mark, initially unset.

```
for i = 0 .. n−1 in order, skipping consumed i:
    base ← a[i]
    for each trial charge u = 2, 3, … while u < z and u < 4:
        step₁ ← f32(Δ₁) ⊘ f32(u)
        scan j = i+1, i+2, … < n (consumed peaks are scanned too):
            d ← (a[j] ⊖ base) ⊖ step₁
            if −τ_d < d < τ_d:                               # strict
                mark j consumed
                a[i] ← f32(u) ⊗ a[i] ⊖ f32(u−1) ⊗ f32(c⁺)    # i.e. (u·a[i]) − ((u−1)·c⁺), each product in f32
                base₂ ← a[j]                                  # value before j is rewritten
                step₂ ← f32(Δ₂) ⊘ f32(u)
                scan k = j+1, … < n:
                    d₂ ← (a[k] ⊖ base₂) ⊖ step₂
                    if −τ_d < d₂ < τ_d: mark k consumed; rewrite a[k] the same way; stop scanning k
                    else if d₂ > τ_d: stop scanning k
                rewrite a[j] the same way
                stop scanning j, and stop trying further charges for this i
            else if d > τ_d: stop scanning j (try the next u)
```

The pseudocode fixes the semantics, not the structure of your code.

- A rewrite maps a u-charged m/z `y` to its singly-charged equivalent: `u·y − (u−1)·c⁺`.
- Because the scans read rewritten values, a peak that was moved to a large m/z earlier ends a later
  scan early. That early stop is part of the specified behaviour.
- Consumed peaks stay in the list, with their rewritten m/z.

Afterwards, re-sort all peaks stably by (m/z, then intensity) ascending. This is the **final peak
list** `Π`. If the flag is not set, `Π` is the working order from §3.4.

**Invariants of Π:**

- it has the same length and the same multiset of (intensity, rank) as the input;
- it is sorted by m/z;
- every intensity inside a suppression window of §3.3 is 0.

The test vectors carry Π for about 75 spectra (§11).

---

## 4. Per-spectrum scoring context

### 4.1 Partition for each segment

For each segment `s ∈ [0, S)`, find a partition `π(s)`. Lexicographic "floor" on `(charge, segment, boundary)`:

1. `F(zq, s) =` the last partition in sorted order whose key `(z_π, s_π, B_π)` is ≤ `(zq, s, M)`,
   comparing lexicographically, with B compared as f32.
2. Compute `F(z, s)`.
   - If there is none: let `z₀` be the charge of the first partition and use `F(z₀, s)`.
   - If there is one and its charge is `z`: use it.
   - If there is one with a different charge `z'`: use `F(z', s)`.

With the bundled model, the effect is this:

- The charge is the largest trained charge ≤ z. If z is below every trained charge, the smallest
  trained charge is used: z = 1 → 2, and z ≥ 7 → 6.
- The partition is the one of that charge and segment with the largest boundary `B ≤ M`.

If `π(s)` has no ion types, segment s contributes nothing.

### 4.2 Peak match rule

Used everywhere a theoretical m/z `x` is looked up in Π:

- The window is `[x ⊖ t, x ⊕ t]` with `t = tol(x)`, inclusive at both ends.
- The **matched peak** is the peak in the window with the greatest intensity. If several share it,
  take the one that comes **last** in Π, i.e. the highest m/z.
- If the window is empty, there is no match.

A zero-intensity peak is a valid match. Suppressed precursor peaks can therefore be matched, and
they carry large ranks.

### 4.3 Rank log-likelihoods

For partition π, ion type t of charge c, and bin b ∈ [0, Rmax]:

```
λ(π,t,b) = f32( ln( f64( row_t[b] ⊘ ( noise[b] ⊗ f32(min(c, S)) ) ) ) )
```

The bin of a lookup:

- matched peak with rank r ≤ Rmax → b = r − 1;
- matched peak with r > Rmax → b = Rmax − 1, so every deep rank shares the last observed-rank bin;
- no match → b = Rmax (the "missing ion" score).

This is the paper's RankScore(ion, rank), with RankScore(ion, ∞) as the missing-ion case.

### 4.4 Peak density ρ and the edge-term tables

Peak density, all in f32:

```
ρ = f32(max(n, 1)) ⊘ max( (M ⊖ f32(W)) ⊘ (ρ-scale ⊗ 2), 1 )
```

- n counts **all** peaks of Π, including zero-intensity ones.
- ρ-scale is the raw tolerance number from the header. For the bundled model it is 0.5, so the
  divisor is the peptide mass in Da.

The **edge partition** is `π(S−1)`, the partition of the last segment.

**Existence scores.** Index `i ∈ {0, 1, 2, 3}` means: bit 0 set = "the current node has an anchor
peak", bit 1 set = "the previous node has an anchor peak". The null probabilities, in f32:

- `q₀ = (1 ⊖ ρ) ⊗ (1 ⊖ ρ)`
- `q₁ = q₂ = ρ ⊗ (1 ⊖ ρ)`
- `q₃ = ρ ⊗ ρ`

Then `ε_i = f32( ln( f64(existence_i) / f64(q_i) ) )`, using the edge partition's existence values.

When ρ > 1, `q₁` and `q₂` are negative, so ε₁ and ε₂ are **NaN**. This happens on very dense
spectra, n > peptide mass in Da. See §6.2 for how NaN becomes an integer.

**Error scores.** For `j ∈ [0, 2E]`: `η_j = f32( ln( f64(signal_j) / f64(noise_j) ) )`, from the
edge partition's error tables.

**Models with E = 0** have no error tables, and edge scoring is undefined in the reference. Such
models are out of scope. An implementation may reject them.

### 4.5 Anchor ion type

Aggregate over the partitions `π(0) … π(S−1)`, one per segment:

- Sum the frequencies φ of ion types in f32, grouping types that have the same (name, exact offset
  bits).
- The **anchor type** A is the group with the largest sum. Ties do not occur with the bundled
  model; if you need a tie rule, break ties by first appearance.

With the bundled model, A is always a singly-charged suffix type:

- either the y ion `S_1_19` (offset ≈ 19.0178),
- or its ¹³C isotope `S_1_20` (offset ≈ 20.0212).

Which one it is depends on the charge and the parent-mass partition. A sweep over charges 2–6 and
M = 300–12000 Da gives roughly half each. Do **not** hard-code y.

**Anchor mass** of nominal node k:

- `α(0) = 0.0`;
- otherwise, let `x = μ(k) ⊘ f32(c_A) ⊕ o_A`. If a peak matches at x (§4.2),
  `α(k) = (mz_peak ⊖ o_A) ⊗ f32(c_A)`. If no peak matches, `α(k) = −1`.
- If there is no anchor type at all, `α(k) = −1` for every k > 0.

A node is **anchored** iff `α(k) ≥ 0`. Node 0 is therefore always anchored.

---

## 5. Vertex (node) scores

For nominal mass `k ≥ 1` and polarity `side ∈ {prefix, suffix}`:

```
Site(k, side) = f32-sum, starting from 0.0, over segments s = 0 … S−1 ascending,
                over the ion types t of π(s) in list order with prefix?(t) matching side,
                of λ(π(s), t, bin(t, k))  — but only for (s, t) where  seg(x) = s
    where  x = μ(k) ⊘ f32(c_t) ⊕ o_t
           seg(x) = min( T( (x ⊘ M) ⊗ f32(S) ), S − 1 )
           bin(t, k) = the bin of the §4.2 lookup at x (§4.3)
```

So each ion type contributes only in the segment where its own theoretical m/z falls, using that
segment's partition. The "segment" is where the fragment's m/z sits relative to the parent mass M.

Define `Pre(k) = Site(k, prefix)` and `Suf(k) = Site(k, suffix)` for k ≥ 1, with `Pre(0) = Suf(0) = 0`.

- Both depend only on the spectrum, the charge and the model, not on any candidate.
- Compute them once per (spectrum, charge) for `k ∈ [0, N₀ − ι_lo]`, where `ι_lo` is the lower
  isotope bound (§8.1).
- Candidates of other masses (§7) may need values outside that range. Those are computed the same
  way.

**Combined vertex score** of a cleavage that splits a peptide of nominal mass P at prefix mass a:

```
V_P(a) = R( Pre(a) ⊕ Suf(P − a) )
```

Precision note: a sum of up to ~12 f32 terms per site, rounded once. Changing the summation order
moves the f32 result by an ulp or so. That flips R(·) only when the value is within ~1e-6 of a
half-integer, which in practice is ≪ 0.1 % of vertices.

---

## 6. Edge scores (high-resolution term)

### 6.1 Definition

An edge joins a "previous" node k′ to a "current" node k and carries a theoretical residue mass
`m_a` (f32):

```
i  = [α(k) ≥ 0] + 2·[α(k′) ≥ 0]
if i = 3:  δ = (α(k) ⊖ α(k′)) ⊖ m_a
           j = clamp( R(δ ⊗ f32(E)), −E, E ) + E
           value = ε₃ ⊕ η_j
else:      value = ε_i
EdgeScore(k, k′, m_a) = R(value)      (an integer)
```

With E = 100, the error bins are 0.01 Da wide, as in the Nat Commun 2014 spectral DAG. Only
edges between two anchored nodes consult the mass error. Every other edge gets one of three
per-spectrum constants.

### 6.2 NaN rule

If `value` is NaN, `EdgeScore = 0`. This happens for i ∈ {1, 2} when ρ > 1 (§4.4). The test
vectors include such a spectrum (synthetic scan 9014).

---

## 7. RawScore of a candidate

A candidate is a residue string `r₁ … r_L` (L ≥ 1). Each position carries a residue letter and a
total modification delta (binary64, 0 if unmodified). Fixed and variable modifications are both
just deltas. There is also a **terminal context** (§7.3).

### 7.1 Masses

- `m_i = mres(r_i) + delta_i`, where mres is the §1.1 residue mass (binary64).
- `n_i = ν(m_i)`, meaning `Rn(f32(m_i) ⊗ κ)`.
- `N_j = n₁ + … + n_j` (integers), with `N_0 = 0` and `P = N_L`.
- `A_j = m₁ + … + m_j`, accumulated left to right in binary64, with `A_0 = 0`.

A **modified position** is one with `delta_i ≠ 0`.

### 7.2 Match score (vertices + edges)

If L = 1, the match score is 0. Otherwise it is:

```
Match = Σ_{j=1}^{L−1} V_P(N_j)                                         (vertices)
      + Σ edges, in the direction given by the anchor type:
          anchor is a suffix type:   for j = 1 … L−1:
                 EdgeScore( k = P − N_j,  k′ = P − N_{j+1},  m_a = f32(A_{j+1} − A_j) )
          anchor is a prefix type:   for j = 1 … L−1:
                 EdgeScore( k = N_j,      k′ = N_{j−1},      m_a = f32(A_j − A_{j−1}) )
      + 0 × (number of modified positions)        (modification penalty is zero)
```

- `m_a` is the binary64 difference of the cumulative sums, rounded to f32. It is **not** f32 of
  `m_i` directly. The two can differ in the last bit.
- For the suffix direction, the edges cover residues r₂ … r_L. The edge of r_L starts at node 0,
  which is always anchored.
- For the prefix direction, the edges cover r₁ … r_{L−1}.
- In both cases, the one residue whose edge would touch the full-mass node P is left out. This
  matches the graph in §8, where edges into the sink score 0.
- Pre/Suf are evaluated at this candidate's own P. P may differ from N₀ because of isotope error,
  modifications or a wrong precursor.

### 7.3 Terminal (enzymatic) terms

The terminal terms depend on the enzyme's cleavage residue set `𝒦`, which is {K, R} for trypsin,
a credit `κ⁺ = +2` and a penalty `κ⁻ = −11`.

- **C-terminal term** `T_C`: κ⁺ if `r_L ∈ 𝒦`, else κ⁻. A peptide ending at the protein C-terminus
  with a non-𝒦 last residue still gets **κ⁻**. The protein end does not substitute.
- **N-terminal term** `T_N`: κ⁺ if any of the following holds, else κ⁻:
  - the residue before r₁ is in 𝒦;
  - r₁ is the first residue of its protein;
  - r₁ is the second residue of its protein and the first residue is M (excised initiator Met).

**RawScore = Match + T_N + T_C.**

The two oracle modes express the N-terminal context differently:

- *Search mode* (§11) uses the real protein context and applies the initiator-Met rule.
- *Rescore mode* takes the context from the peptide string `X.PEPTIDE.Y`:
  - `X ∈ 𝒦` or `X = '-'` gives credit, and any other X gives the penalty;
  - a string **without** flanks is treated as N-terminal credit;
  - rescore mode has no initiator-Met rule, so `M.SEGAYQR.L` gets κ⁻.

Your API should take the N-terminal credit decision as a boolean from the caller. The DDA engine
knows the protein context.

**Cleavage disabled** (unspecific enzyme): `T_N = T_C = 0`, and §8 changes as noted there.

---

## 8. Score distribution and SpecEValue

### 8.1 Isotope sinks

The isotope-error range is `[ι_lo, ι_hi]`, with default `[0, 1]`. The candidate nominal masses
("sinks") are:

```
p ∈ { N₀ − ι_hi, …, N₀ − ι_lo }, keeping only p > 0
```

An isotope error of +1 means the precursor was picked one ¹³C peak high, so the true nominal
mass is one lower. The reference builds sinks from the isotope range only. It does **not** add
sinks for the precursor-tolerance window.

### 8.2 Residue alphabet

The alphabet is an ordered list of entries `(letter, real mass m_a, nominal n_a = ν(m_a), probability w_a)`.

- The letter is the unmodified residue, used for the cleavage test.
- `m_a` includes any modification delta.
- Letters may repeat, e.g. an oxidised-M entry beside M. Nominal masses may coincide, e.g. I/L and
  K/Q. Each entry is a separate edge.
- The mixture weight of §8.5 is `π = Σ_{letters ℓ ∈ 𝒦} w(ℓ)`, where w(ℓ) is the background
  probability of the residue letter ℓ. That is `w_K + w_R` for trypsin. Each letter counts once,
  even if the alphabet also has modified entries for it.

The oracle alphabets are listed in each test-vector README.

- **Uniform (rescore default):** the 20 residues, each `w = 0.05`, unmodified masses (C = 103.00919),
  in the order `G A S P V T C L I N D Q K E M H F R Y W`. Then w_𝒦 = 0.1.
- **Rescore with composition and oxidation:** the same 20 letters, with probabilities from a supplied
  table (0.05 for any letter the table lacks), plus one entry `M + 15.994915` carrying w_M.
- **Search mode:**
  - the 20 residues at their *fixed-modified* masses (e.g. C + 57.021464);
  - then, for each position-unrestricted variable modification and each residue it targets, one
    entry at the modified mass carrying the base residue's probability;
  - probabilities are the target-database composition (count of each of the 20 letters ÷ total
    count of the 20 letters; decoy proteins not counted; missing letters get 0);
  - position-restricted variable mods (N-term etc.) are **not** added to the alphabet.

### 8.3 The amino-acid graph for one sink p

Vertices are the integers `0 … p`. Vertex 0 is the **source** and p is the **sink**.

**Vertex weight:**

- `v(0) = v(p) = 0`;
- for `0 < m < p`: `v(m) = V_p(p − m) = R( Pre(p − m) ⊕ Suf(m) )`.

The graph index m is a **suffix** (C-terminal) nominal mass. Paths run from the C-terminus
towards the N-terminus.

**Edges:** for every vertex `m ∈ [1, p]` and every alphabet entry a, in alphabet order, with
`m′ = m − n_a ≥ 0`, there is an edge `m′ → m`. Its weight:

- `m < p`: `e = EdgeScore(k = m, k′ = m′, m_a) + [m′ = 0]·T(a)`, where
  - T(a) = κ⁺ if letter(a) ∈ 𝒦, else κ⁻;
  - T(a) is 0 if cleavage is disabled.
  - The source edge carries the peptide's C-terminal cleavage term, and the anchor lookups are
    α(m) and α(m′) of §4.5.
- `m = p`: `e = 0`. Edges into the sink carry no edge score and no cleavage term.

The probability of an edge is `w_a`. A path from 0 to p spells a residue string of nominal mass p,
read from the C-terminus. Its **score** is the sum over its non-source vertices of `v(m) + e(in-edge of m)`.

For the true peptide with P = p and a suffix-type anchor, this path score equals `Match + T_C` of
§7, up to the last-bit difference between `f32(m_a)` and `f32(A_{j+1} − A_j)`. That identity is a
useful self-test.

### 8.4 Dynamic program

Let `D_m` be a distribution: a map from integer score to probability.

```
D_0 = {0 ↦ 1}
for m = 1 … p (increasing):
    D_m(σ) = Σ_{edges m′→m, D_{m′} non-empty}  w_a · D_{m′}(σ − v(m) − e)
```

- A vertex with no reachable predecessor has an empty D_m.
- Use binary64 for probabilities.
- The support of D_m is the contiguous integer range from the minimum to the maximum shifted
  predecessor support. **Structural** support means cells with probability exactly 0 still count
  toward the range. This matters for DeNovoScore only if some `w_a = 0`.
- If `D_p` is empty, the sink is unreachable and contributes nothing.

### 8.5 Neighbour-cleavage mixture

The residue *before* the peptide is unknown under the null. Mix, with π = w_𝒦:

```
G_p(σ) = π · D_p(σ − κ⁺) + (1 − π) · D_p(σ − κ⁻)
```

If cleavage is disabled, `G_p = D_p`.

### 8.6 Merge, DeNovoScore, SpecEValue

- `H = Σ_p G_p` over the reachable sinks. This is a sum, not an average, so the total mass can
  exceed 1.
- If no sink is reachable, the (spectrum, charge) has **no result** at all. That happens when
  N₀ − ι_hi … N₀ − ι_lo are not reachable with the alphabet.
- **DeNovoScore** = the largest σ in the structural support of H. Equivalently, it is the maximum
  over reachable sinks of (the best path score to p) + κ⁺, or + 0 if cleavage is disabled.
- **SpecEValue(RawScore)** = `min(1, Σ_{σ ≥ RawScore} H(σ))`.
  - If RawScore is above the support, the sum is empty and SpecEValue = 0. The reference prints
    `-0.000000e0`, a negative zero; sign is irrelevant.
  - This happens whenever RawScore > DeNovoScore, e.g. a candidate with a residue mass that is
    absent from the alphabet.
  - If RawScore is below the support, the whole mass is summed and capped at 1.
- **EValue** (engine level, optional) = SpecEValue × (number of candidate peptides searched).
  It does not affect ranking or FDR.

DeNovoScore and H depend only on (spectrum, charge, model, isotope range, alphabet, cleavage
rule). Build them **once per (spectrum, charge)** and read every candidate's SpecEValue off the
same H.

### 8.7 Exact tail pruning (optional optimisation, recommended)

A caller that knows the lowest RawScore θ it will query can skip work below it. For each vertex m,
let `best_rem(m)` be the best score attainable from m to the sink. That is one backward integer
pass over the graph, with sink in-edges scoring 0.

- A cell (m, σ) can be dropped when `σ + best_rem(m) < θ − κ⁺`.
- Cells at or above the threshold are then computed by exactly the same operations as without
  pruning, so the tail values are **identical**.
- Clamp the cut at the best full-path score, so DeNovoScore stays exact.
- The reference uses this with θ = the minimum RawScore of the candidates sharing a
  (spectrum, charge).

### 8.8 N-terminal enzymes

The graph is built in the C-terminal direction, which suits trypsin-like enzymes. For enzymes
that cleave N-terminal to their residues (Lys-N, Asp-N), the reference **disables** cleavage
scoring entirely (§7.3, §8.3, §8.5). Match this behaviour. Proper support is future work.

---

## 9. API boundary (proposal)

The following interface gives a DDA engine everything it needs. Names are suggestions; types are
Rust-flavoured.

```text
ScoreModel::from_bytes(&[u8]) -> Result<ScoreModel, ModelError>      // §2; reject bad sentinel / E = 0
ScoreModel::bundled() -> ScoreModel                                   // the CC0-trained default

struct Ms2 { peaks: &[(f64 /*mz*/, f64 /*intensity*/)], precursor_mz: f64, charge: u8 }
fn prepare(model, &Ms2) -> Option<PreparedSpectrum>                   // §3–§5; None if N₀ ∉ [50, 10000]
    PreparedSpectrum::peaks() -> &[(f32 mz, f32 intensity, u32 rank)] // Π, for tests/debug
    PreparedSpectrum::n0() -> i32                                     // N₀

struct Residue { letter: u8, delta: f64 }                             // delta includes fixed + variable mods
struct Candidate<'a> { residues: &'a [Residue], n_term_credit: bool }
fn match_and_terminal_score(&PreparedSpectrum, &Candidate, &Cleavage) -> i32         // §7

struct AlphabetEntry { letter: u8, mass: f64, prob: f64 }
struct Cleavage { sites: Vec<u8>, credit: i32 /*+2*/, penalty: i32 /*-11*/, enabled: bool }
struct NullModel { alphabet: Vec<AlphabetEntry>, cleavage: Cleavage, isotope: (i32, i32) }
fn score_distribution(&PreparedSpectrum, &NullModel, min_query: Option<i32>) -> Option<NullTail>  // §8
    NullTail::best_possible() -> i32  // DeNovoScore
    NullTail::tail_mass(raw: i32) -> f64    // SpecEValue; assert raw >= min_query when pruned
```

Guidance:

- `PreparedSpectrum` should cache `Pre`, `Suf` and `α` over `[0, max sink]`. They are
  candidate-independent.
- The per-spectrum edge constants for i ∈ {0, 1, 2} can be cached too.
- `match_and_terminal_score` must work for any P, including P > max sink.
- Thread safety: everything is immutable after construction. Keep one DP scratch arena per worker.
- Performance reference: the black box rescores about 4,000 PSMs/s single-threaded on this data
  (17.7 k PSMs over about 3 k spectrum–charge groups in about 4 s). With exact pruning, a DP is
  about 1 ms per (spectrum, charge) on the 0.5 Da / ~1 Da nominal grid.

---

## 10. Edge cases (observed; all are in the test vectors)

| case | behaviour |
|---|---|
| empty peak list | scored. ρ uses max(n, 1). Every lookup misses, so every site takes its missing-ion bin. RawScore and DeNovoScore are negative. SpecEValue is still defined (synthetic 9009). |
| single peak | as above, with one possible match (9008) |
| N₀ < 50 or > 10000 | no output for that (spectrum, charge) (9010) |
| charge 1 | the charge-2 partitions are used (§4.1). The reduce step does nothing, because u < z is never true (9006) |
| charge above the trained charges (7, 8) | the charge-6 partitions are used (9007). Reduce tries u = 2, 3 |
| unsorted input peaks | sorted per §3.2 (9003) |
| duplicate m/z, zero intensity, m/z 0 | kept. Zeros rank last, with higher m/z first among equal intensities (9011) |
| equal intensities everywhere | ties broken by m/z descending (9002) |
| ρ > 1 (very dense spectrum) | ε₁ and ε₂ are NaN, so those edges score 0 (9014) |
| residue mass not in the alphabet (e.g. C+57 with the uniform alphabet, M+16 without the oxidation entry) | RawScore can exceed DeNovoScore, and then SpecEValue = 0 (9004, 9005, 9012 z = 3) |
| isobaric swap (I/L) | identical scores |
| non-standard letter, or a modification written before the first residue (`+42.011PEPTIDE`) | unscorable in the oracle's peptide parser. Your API takes residues directly; express N-terminal mods as a delta on r₁ |
| precursor one ¹³C high | the true peptide sits in sink N₀ − 1, so the SpecEValue still reflects it (9013) |
| candidate P ≠ any sink | RawScore is computed at its own P. SpecEValue is read off H as usual |
| modification deltas printed with 3 decimals | parsed literally. A 0.0005 Da change can flip an error bin, and a few PSMs move by 1–9 points between "search" (full-precision delta) and "rescore" (3-decimal delta) |

---

## 11. Conformance

Acceptance criteria against the oracle, all with the bundled model:

1. **Model decode.** All header scalars and table sizes as in §2. Sentinel verified.
2. **Preprocessing.** For every spectrum in `test_vectors/intermediate/*.jsonl`, the m/z (f32),
   intensity and rank arrays equal the recorded Π exactly (same length and order). Allow ≤ 1 f32
   ulp on rewritten m/z values.
3. **Per-spectrum context.**
   - The segment partitions equal the recorded ones.
   - ρ is within 1 ulp.
   - The three edge constants (neither, current-only, previous-only anchored) are equal.
   - The both-anchored edge integers over the recorded error sweep are equal.
   - `Pre`, `Suf` and `α` agree within 1e-4 absolute, and exactly in ≥ 99.9 % of entries.
4. **Integer scores.** RawScore and DeNovoScore are equal for ≥ 99.5 % of the PSMs in
   `test_vectors/rescore_*` and `search_*`. Every mismatch is off by 1, and each is attributable to
   an f32 rounding tie; list them.
5. **SpecEValue.**
   - Spearman ρ ≥ 0.99 on log₁₀ SpecEValue per test set.
   - Where the integer scores agree, |log₁₀ ours − log₁₀ oracle| ≤ 0.01. Both being 0 counts as a
     match.
   - The full merged distributions in `intermediate/` agree cell-wise to relative 1e-9 where
     > 1e-300.
6. **Identifications.** Take `test_vectors/search_hela_r01` and the top-1 PSM per spectrum, ordered
   by SpecEValue. Target-decoy q ≤ 1 % should give the same PSM set up to ±0.5 %.

**How to measure:** `test_vectors/README.md` describes every file. `oracle_bin/USAGE.md` explains
how to run the oracle on new inputs.

---

## 12. Provenance of each rule

**From the published papers** (MS-GF 2008, MS-GFDB 2010, MS-GF+ 2014 main text and methods):

- the spectral-vector and rank-score construction;
- ion types as (charge, offset, prefix/suffix);
- the missing-ion score;
- the 0.9995-type nominal rescaling (the reference uses 0.999497);
- the log-likelihood vertex and edge scores of the spectral DAG;
- the 0.01 Da error bins;
- the amino-acid graph and generating-function DP;
- spectral probability as the tail mass;
- EValue as SpecEValue × database size;
- trained parameters conditioned on charge, parent mass and m/z segment;
- the enzyme-dependent scoring of peptide termini.

**Fixed from observed reference behaviour** (black-box outputs and MSGF_Rust documentation). The
papers do not determine these:

- all f32/f64 precisions and the three rounding functions (§1.2);
- the NaN → 0 rule;
- the precursor-suppression rule: which table entries, the centre formula, the use of the global
  tolerance, the inclusive window;
- that ranks are computed after suppression and before isotope reduction, and the tie order;
- the exact isotope-reduction procedure, including live rewritten values, the second-isotope
  constant Δ₂ = 0.99988617, the charge limits u < z and u < 4, the strict tolerance and no peak
  removal;
- the partition floor rule and the charge fallback;
- segment assignment by fragment m/z relative to M, with truncation and a clamp to S−1;
- the peak-match rule (max intensity, last-in-order tie-break);
- deep ranks sharing bin Rmax − 1;
- the `min(charge, S)` factor in the rank log-likelihood;
- existence values of 0 floored to 0.001;
- the ρ formula using the raw tolerance number;
- edge-partition = last segment's partition;
- the anchor-type choice by summed frequency;
- α(0) = 0 counting as anchored;
- the edge direction and the omitted sink-adjacent edge;
- the credit/penalty values +2 / −11;
- T_C gets no protein-end credit, while T_N gets protein-start and initiator-Met credit;
- sink in-edges scoring 0;
- the neighbour-cleavage mixture with π = w_K + w_R;
- isotope sinks only, with no precursor-tolerance sinks;
- summing (not averaging) the per-sink distributions;
- the structural DeNovoScore;
- SpecEValue capped at 1;
- the per-sink uniform/composition alphabets;
- the unsupported N-terminal enzymes.
