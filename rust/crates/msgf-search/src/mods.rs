//! Fixed and variable modifications, and enumeration of a peptide's modified forms.
//!
//! Specs follow the documented `Mods.txt` interchange format:
//!
//! ```text
//! NumMods=3
//! C2H3N1O1,C,fix,any,Carbamidomethyl
//! O1,M,opt,any,Oxidation
//! ```
//!
//! i.e. `Mass-or-CompositionStr, Residues, ModType, Position, Name`, where `Residues` is a set of
//! upper-case residues or `*` (any), `ModType` is `fix`/`opt`, and `Position` is one of
//! `any`, `N-term`, `C-term`, `Prot-N-term`, `Prot-C-term`. `#` starts a comment.
//!
//! Fixed mods are folded into the residue mass everywhere they apply. Variable mods generate the
//! combinatorial set of placements, capped at [`ModSet::max_var_mods`] per peptide (`NumMods`).

use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// Where a modification may be attached.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ModPosition {
    /// Anywhere in the peptide.
    Any,
    /// The peptide's N-terminal residue.
    PeptideNTerm,
    /// The peptide's C-terminal residue.
    PeptideCTerm,
    /// The N-terminal residue of the protein (only when the peptide starts the protein).
    ProteinNTerm,
    /// The C-terminal residue of the protein (only when the peptide ends the protein).
    ProteinCTerm,
}

impl ModPosition {
    fn parse(s: &str) -> Result<ModPosition, String> {
        // "-" is optional and matching is case-insensitive, per the format's own documentation.
        let k: String = s
            .chars()
            .filter(|c| *c != '-' && !c.is_whitespace())
            .flat_map(char::to_lowercase)
            .collect();
        Ok(match k.as_str() {
            "any" => ModPosition::Any,
            "nterm" => ModPosition::PeptideNTerm,
            "cterm" => ModPosition::PeptideCTerm,
            "protnterm" => ModPosition::ProteinNTerm,
            "protcterm" => ModPosition::ProteinCTerm,
            _ => return Err(format!("unknown modification position `{s}`")),
        })
    }

    /// Whether this position is a terminus (as opposed to anywhere in the peptide).
    #[inline]
    pub fn is_terminal(&self) -> bool {
        !matches!(self, ModPosition::Any)
    }
}

/// One modification specification.
#[derive(Debug, Clone, PartialEq)]
pub struct ModSpec {
    /// Monoisotopic delta in Da.
    pub mass: f64,
    /// Residues this mod applies to; **empty means `*`** (any residue).
    pub residues: Vec<u8>,
    pub is_fixed: bool,
    pub position: ModPosition,
    pub name: String,
}

impl ModSpec {
    /// Whether this mod can attach to `residue` at peptide index `idx` of a peptide of length
    /// `len`, given whether the peptide starts/ends its protein.
    #[inline]
    pub fn applies(&self, residue: u8, idx: usize, len: usize, prot_n: bool, prot_c: bool) -> bool {
        if !self.residues.is_empty() && !self.residues.contains(&residue) {
            return false;
        }
        match self.position {
            ModPosition::Any => true,
            ModPosition::PeptideNTerm => idx == 0,
            ModPosition::PeptideCTerm => idx + 1 == len,
            ModPosition::ProteinNTerm => idx == 0 && prot_n,
            ModPosition::ProteinCTerm => idx + 1 == len && prot_c,
        }
    }

    /// Parse one `Mass-or-CompositionStr, Residues, ModType, Position, Name` line.
    pub fn parse(line: &str) -> Result<ModSpec, String> {
        let f: Vec<&str> = line.split(',').map(str::trim).collect();
        if f.len() < 5 {
            return Err(format!(
                "modification `{line}`: expected Mass-or-Composition,Residues,ModType,Position,Name"
            ));
        }
        let mass = parse_mass_or_composition(f[0])
            .ok_or_else(|| format!("modification `{line}`: bad mass/composition `{}`", f[0]))?;
        let residues = if f[1] == "*" {
            Vec::new()
        } else {
            if !f[1].bytes().all(|b| b.is_ascii_uppercase()) {
                return Err(format!(
                    "modification `{line}`: Residues must be upper-case residues or `*`"
                ));
            }
            f[1].as_bytes().to_vec()
        };
        let is_fixed = match f[2].to_ascii_lowercase().as_str() {
            "fix" | "fixed" => true,
            "opt" | "optional" | "variable" => false,
            "custom" => {
                return Err(format!(
                    "modification `{line}`: custom amino acids are not supported yet"
                ))
            }
            other => return Err(format!("modification `{line}`: unknown ModType `{other}`")),
        };
        let position =
            ModPosition::parse(f[3]).map_err(|e| format!("modification `{line}`: {e}"))?;
        Ok(ModSpec {
            mass,
            residues,
            is_fixed,
            position,
            name: f[4..].join(",").trim().to_string(),
        })
    }

    /// Parse a convenience shorthand `<residues><+|-><mass>` (e.g. `C+57.021464`, `M+15.994915`),
    /// or a full comma-separated spec. `is_fixed` sets the type for the shorthand form.
    pub fn parse_short(spec: &str, is_fixed: bool) -> Result<ModSpec, String> {
        if spec.contains(',') {
            return ModSpec::parse(spec);
        }
        let i = spec.find(['+', '-']).ok_or_else(|| {
            format!("modification `{spec}`: expected <residues>+<mass>, e.g. C+57.021464")
        })?;
        let (res, delta) = spec.split_at(i);
        let mass: f64 = delta
            .parse()
            .map_err(|_| format!("modification `{spec}`: bad mass `{delta}`"))?;
        let residues = if res == "*" || res.is_empty() {
            Vec::new()
        } else {
            if !res.bytes().all(|b| b.is_ascii_uppercase()) {
                return Err(format!(
                    "modification `{spec}`: residues must be upper-case"
                ));
            }
            res.as_bytes().to_vec()
        };
        Ok(ModSpec {
            mass,
            residues,
            is_fixed,
            position: ModPosition::Any,
            name: format!("{res}{delta}"),
        })
    }
}

/// Parse either an explicit monoisotopic mass (`15.994915`) or an elemental composition string
/// (`C2H3N1O1`, `H-1N-1O1`). Element counts may be negative; a bare element means one atom.
fn parse_mass_or_composition(s: &str) -> Option<f64> {
    if let Ok(m) = s.parse::<f64>() {
        return Some(m);
    }
    let b = s.as_bytes();
    let mut i = 0usize;
    let mut mass = 0.0f64;
    while i < b.len() {
        if !b[i].is_ascii_uppercase() {
            return None;
        }
        // Two-letter elements first (Br, Cl, Fe, Se), then single-letter.
        let (elem_len, unit) = match (b[i], b.get(i + 1)) {
            (b'B', Some(b'r')) => (2, 78.918_337_6),
            (b'C', Some(b'l')) => (2, 34.968_852_68),
            (b'F', Some(b'e')) => (2, 55.934_936_3),
            (b'S', Some(b'e')) => (2, 79.916_521_8),
            (b'C', _) => (1, msgf_chem::mass::C),
            (b'H', _) => (1, msgf_chem::mass::H),
            (b'N', _) => (1, msgf_chem::mass::N),
            (b'O', _) => (1, msgf_chem::mass::O),
            (b'S', _) => (1, msgf_chem::mass::S),
            (b'P', _) => (1, msgf_chem::mass::P),
            _ => return None,
        };
        i += elem_len;
        let start = i;
        if i < b.len() && b[i] == b'-' {
            i += 1;
        }
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
        }
        let count: i32 = if i == start {
            1 // a bare element symbol means one atom
        } else {
            s[start..i].parse().ok()?
        };
        mass += count as f64 * unit;
    }
    Some(mass)
}

/// The full modification configuration for a search.
#[derive(Debug, Clone)]
pub struct ModSet {
    pub mods: Vec<ModSpec>,
    /// `NumMods` — the maximum number of **variable** modifications on one peptide.
    pub max_var_mods: usize,
}

impl Default for ModSet {
    fn default() -> Self {
        Self {
            mods: Vec::new(),
            max_var_mods: 2,
        }
    }
}

impl ModSet {
    /// Read a `Mods.txt`-format configuration file.
    pub fn read_file(path: &Path) -> io::Result<ModSet> {
        let mut set = ModSet::default();
        for (n, line) in BufReader::new(File::open(path)?).lines().enumerate() {
            let line = line?;
            let body = line.split('#').next().unwrap_or("").trim();
            if body.is_empty() {
                continue;
            }
            let bad = |e: String| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("{}:{}: {e}", path.display(), n + 1),
                )
            };
            if let Some(v) = body.strip_prefix("NumMods=") {
                set.max_var_mods = v
                    .trim()
                    .parse()
                    .map_err(|_| bad(format!("bad NumMods `{}`", v.trim())))?;
                continue;
            }
            set.mods.push(ModSpec::parse(body).map_err(bad)?);
        }
        Ok(set)
    }

    pub fn fixed(&self) -> impl Iterator<Item = &ModSpec> {
        self.mods.iter().filter(|m| m.is_fixed)
    }

    pub fn variable(&self) -> impl Iterator<Item = (usize, &ModSpec)> {
        self.mods.iter().enumerate().filter(|(_, m)| !m.is_fixed)
    }

    /// Total fixed-mod delta on `residue` at peptide index `idx`, summed over every applicable
    /// fixed mod (MS-GF+ allows a residue-level and a terminal fixed mod to stack).
    #[inline]
    pub fn fixed_delta(
        &self,
        residue: u8,
        idx: usize,
        len: usize,
        prot_n: bool,
        prot_c: bool,
    ) -> f64 {
        self.fixed()
            .filter(|m| m.applies(residue, idx, len, prot_n, prot_c))
            .map(|m| m.mass)
            .sum()
    }

    /// The fixed-mod delta applied to `residue` **everywhere it occurs** — i.e. only the
    /// position-independent (`any`) fixed mods. This is what folds into the de novo graph's
    /// residue masses, where per-position terminal mods have no representation.
    #[inline]
    pub fn fixed_residue_delta(&self, residue: u8) -> f64 {
        self.fixed()
            .filter(|m| m.position == ModPosition::Any)
            .filter(|m| m.residues.is_empty() || m.residues.contains(&residue))
            .map(|m| m.mass)
            .sum()
    }
}

/// One placed variable modification: which residue of the peptide carries which mod.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct PlacedMod {
    /// Index of the modified residue within the peptide.
    pub pos: u8,
    /// Index into [`ModSet::mods`].
    pub mod_idx: u8,
}

/// Maximum number of variable modifications carried on one candidate. `NumMods` above this is
/// clamped (and reported) rather than silently allowing an unbounded candidate explosion.
pub const MAX_PLACED_MODS: usize = 4;

/// The variable-mod placements of one candidate peptide, stored inline (no per-candidate heap
/// allocation — a proteome search holds tens of millions of these).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ModPlacement {
    pub n: u8,
    pub slots: [PlacedMod; MAX_PLACED_MODS],
}

impl ModPlacement {
    #[inline]
    pub fn placed(&self) -> &[PlacedMod] {
        &self.slots[..self.n as usize]
    }

    /// The mod delta on peptide residue `pos`, if any.
    #[inline]
    pub fn delta_at(&self, pos: usize, set: &ModSet) -> f64 {
        self.placed()
            .iter()
            .filter(|p| p.pos as usize == pos)
            .map(|p| set.mods[p.mod_idx as usize].mass)
            .sum()
    }
}

/// Enumerate every variable-mod placement for one peptide, calling `emit(placement, delta_sum)`.
/// The unmodified form is always emitted first.
///
/// At most one variable mod may occupy a given residue. Placements are capped at
/// `min(set.max_var_mods, MAX_PLACED_MODS)` mods and at `max_variants` total forms; the return
/// value is the number of forms **not** emitted because of the `max_variants` cap, so a caller can
/// report truncation instead of silently under-searching.
pub fn enumerate_placements(
    seq: &[u8],
    set: &ModSet,
    prot_n: bool,
    prot_c: bool,
    max_variants: usize,
    mut emit: impl FnMut(ModPlacement, f64),
) -> usize {
    let len = seq.len();
    // Sites: every (residue index, variable mod) pair that could be placed.
    let mut sites: Vec<PlacedMod> = Vec::new();
    for (mod_idx, spec) in set.variable() {
        if mod_idx > u8::MAX as usize {
            break;
        }
        for (i, &r) in seq.iter().enumerate() {
            if i > u8::MAX as usize {
                break;
            }
            if spec.applies(r, i, len, prot_n, prot_c) {
                sites.push(PlacedMod {
                    pos: i as u8,
                    mod_idx: mod_idx as u8,
                });
            }
        }
    }

    emit(ModPlacement::default(), 0.0);
    let max_k = set.max_var_mods.min(MAX_PLACED_MODS);
    if sites.is_empty() || max_k == 0 {
        return 0;
    }

    // Depth-first over combinations of distinct sites (ascending index), rejecting any pair that
    // would put two variable mods on the same residue. The walk carries its context in a struct so
    // the recursion stays a two-argument step.
    struct Walk<'a, F> {
        sites: &'a [PlacedMod],
        set: &'a ModSet,
        max_k: usize,
        max_variants: usize,
        emitted: usize,
        skipped: usize,
        emit: F,
    }
    impl<F: FnMut(ModPlacement, f64)> Walk<'_, F> {
        fn step(&mut self, from: usize, depth: usize, cur: &mut ModPlacement, delta: f64) {
            if depth == self.max_k {
                return;
            }
            for i in from..self.sites.len() {
                let site = self.sites[i];
                if cur.placed().iter().any(|p| p.pos == site.pos) {
                    continue; // one variable mod per residue
                }
                if self.emitted >= self.max_variants {
                    self.skipped += 1;
                    continue;
                }
                cur.slots[depth] = site;
                cur.n = depth as u8 + 1;
                let d = delta + self.set.mods[site.mod_idx as usize].mass;
                (self.emit)(*cur, d);
                self.emitted += 1;
                self.step(i + 1, depth + 1, cur, d);
                cur.n = depth as u8;
            }
        }
    }

    let mut walk = Walk {
        sites: &sites,
        set,
        max_k,
        max_variants,
        emitted: 1, // the unmodified form, already emitted
        skipped: 0,
        emit: &mut emit,
    };
    walk.step(0, 0, &mut ModPlacement::default(), 0.0);
    walk.skipped
}

#[cfg(test)]
mod tests {
    use super::*;

    fn oxidation() -> ModSpec {
        ModSpec::parse("O1,M,opt,any,Oxidation").unwrap()
    }

    #[test]
    fn parses_composition_and_mass() {
        assert!((oxidation().mass - 15.994_914_622).abs() < 1e-6);
        let carbamidomethyl = ModSpec::parse("C2H3N1O1,C,fix,any,Carbamidomethyl").unwrap();
        assert!((carbamidomethyl.mass - 57.021_463_7).abs() < 1e-6);
        assert!(carbamidomethyl.is_fixed);
        // Negative counts (deamidation) and bare element symbols.
        let deamid = ModSpec::parse("H-1N-1O1,NQ,opt,any,Deamidated").unwrap();
        assert!((deamid.mass - 0.984_016).abs() < 1e-5);
        // Explicit mass instead of a formula.
        let explicit = ModSpec::parse("15.994915,M,opt,any,Oxidation").unwrap();
        assert!((explicit.mass - 15.994915).abs() < 1e-12);
    }

    #[test]
    fn parses_positions_leniently() {
        for s in ["Prot-N-term", "ProtNTerm", "prot-n-Term"] {
            assert_eq!(ModPosition::parse(s).unwrap(), ModPosition::ProteinNTerm);
        }
        assert_eq!(ModPosition::parse("any").unwrap(), ModPosition::Any);
        assert!(ModPosition::parse("middle").is_err());
    }

    #[test]
    fn shorthand_specs() {
        let m = ModSpec::parse_short("C+57.021464", true).unwrap();
        assert_eq!(m.residues, b"C");
        assert!(m.is_fixed);
        assert!((m.mass - 57.021464).abs() < 1e-12);
        let n = ModSpec::parse_short("M-15.5", false).unwrap();
        assert!((n.mass + 15.5).abs() < 1e-12);
        assert!(!n.is_fixed);
    }

    #[test]
    fn position_restricted_application() {
        let nterm = ModSpec::parse("C2H2O,*,opt,Prot-N-term,Acetyl").unwrap();
        assert!(nterm.applies(b'A', 0, 5, true, false));
        assert!(!nterm.applies(b'A', 0, 5, false, false)); // not at the protein N-term
        assert!(!nterm.applies(b'A', 1, 5, true, false)); // not the first residue
        let cterm = ModSpec::parse("1.0,K,opt,C-term,Test").unwrap();
        assert!(cterm.applies(b'K', 4, 5, false, false));
        assert!(!cterm.applies(b'R', 4, 5, false, false)); // wrong residue
    }

    #[test]
    fn enumerates_variable_placements() {
        let set = ModSet {
            mods: vec![oxidation()],
            max_var_mods: 2,
        };
        let mut forms: Vec<(usize, i64)> = Vec::new();
        let skipped = enumerate_placements(b"MAMK", &set, false, false, 10_000, |p, d| {
            forms.push((p.n as usize, (d * 1000.0).round() as i64))
        });
        assert_eq!(skipped, 0);
        // unmodified, M0, M2, M0+M2
        assert_eq!(forms.len(), 4);
        assert_eq!(forms.iter().filter(|(n, _)| *n == 0).count(), 1);
        assert_eq!(forms.iter().filter(|(n, _)| *n == 1).count(), 2);
        assert_eq!(forms.iter().filter(|(n, _)| *n == 2).count(), 1);
    }

    #[test]
    fn num_mods_limits_placements() {
        let set = ModSet {
            mods: vec![oxidation()],
            max_var_mods: 1,
        };
        let mut n = 0;
        enumerate_placements(b"MMM", &set, false, false, 10_000, |_, _| n += 1);
        assert_eq!(n, 4); // unmodified + one per M
    }

    #[test]
    fn variant_cap_is_reported() {
        let set = ModSet {
            mods: vec![oxidation()],
            max_var_mods: 3,
        };
        let mut n = 0;
        let skipped = enumerate_placements(b"MMMMM", &set, false, false, 3, |_, _| n += 1);
        assert_eq!(n, 3);
        assert!(skipped > 0);
    }

    #[test]
    fn one_variable_mod_per_residue() {
        // Two different variable mods on M must never land on the same residue.
        let set = ModSet {
            mods: vec![
                oxidation(),
                ModSpec::parse("32.0,M,opt,any,Dioxidation").unwrap(),
            ],
            max_var_mods: 2,
        };
        let mut seen: Vec<ModPlacement> = Vec::new();
        enumerate_placements(b"MK", &set, false, false, 10_000, |p, _| seen.push(p));
        for p in &seen {
            let mut pos: Vec<u8> = p.placed().iter().map(|x| x.pos).collect();
            pos.sort_unstable();
            let n = pos.len();
            pos.dedup();
            assert_eq!(pos.len(), n, "duplicate residue in {p:?}");
        }
        // unmodified + oxidation on M + dioxidation on M
        assert_eq!(seen.len(), 3);
    }
}
