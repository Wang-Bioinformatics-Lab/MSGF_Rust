//! The candidate peptide index: digest every protein, apply modifications, and sort the resulting
//! candidates by neutral mass so a precursor window is a binary-search slice.
//!
//! Candidates are addressed by `(protein, absolute offset, length)` into [`ProteinDb::seq`] plus an
//! inline [`ModPlacement`], so nothing per-candidate is heap-allocated. Building is parallel over
//! proteins; the sort is the only global step.

use crate::mods::{enumerate_placements, ModPlacement, ModSet};
use msgf_chem::mass;
use msgf_db::enzyme::{digest, DigestParams};
use msgf_db::fasta::ProteinDb;
use rayon::prelude::*;

/// Ceiling on the number of modified forms generated from one digested peptide. A peptide with
/// many modifiable residues would otherwise blow up combinatorially; hitting this is counted in
/// [`PeptideIndex::truncated_peptides`] and reported rather than silently under-searching.
pub const MAX_VARIANTS_PER_PEPTIDE: usize = 4096;

/// One searchable candidate: a peptide sequence plus a specific modification placement.
#[derive(Debug, Clone, Copy)]
pub struct Candidate {
    /// Neutral monoisotopic peptide mass, including water and every fixed and variable mod.
    pub mass: f64,
    /// Index into [`ProteinDb::proteins`].
    pub protein: u32,
    /// Absolute offset of the first residue in [`ProteinDb::seq`].
    pub start: u32,
    pub len: u16,
    /// Number of enzymatic termini (2 = fully enzymatic).
    pub n_termini: u8,
    pub placement: ModPlacement,
}

/// The mass-sorted candidate database.
pub struct PeptideIndex {
    /// Every candidate, ascending by `mass`.
    pub candidates: Vec<Candidate>,
    /// Digested peptides whose modified forms hit [`MAX_VARIANTS_PER_PEPTIDE`].
    pub truncated_peptides: usize,
    /// Distinct digested peptides (before modification placement).
    pub n_peptides: usize,
}

impl PeptideIndex {
    /// Digest `db` and enumerate every modified candidate, sorted by neutral mass.
    pub fn build(db: &ProteinDb, digest_params: &DigestParams, mods: &ModSet) -> PeptideIndex {
        let per_protein: Vec<(Vec<Candidate>, usize, usize)> = (0..db.proteins.len())
            .into_par_iter()
            .map(|pi| {
                let protein = &db.proteins[pi];
                let seq = &db.seq[protein.start..protein.start + protein.len];
                let mut out: Vec<Candidate> = Vec::new();
                let (mut truncated, mut n_pep) = (0usize, 0usize);

                digest(seq, digest_params, |d| {
                    n_pep += 1;
                    let pep = &seq[d.start..d.start + d.len];
                    let prot_n = d.start == 0;
                    let prot_c = d.start + d.len == protein.len;

                    // Base mass: residues + water + every applicable fixed mod.
                    let mut base = mass::WATER;
                    for (i, &r) in pep.iter().enumerate() {
                        base += msgf_chem::residue_mass(r)
                            .expect("digest rejects non-standard residues")
                            + mods.fixed_delta(r, i, d.len, prot_n, prot_c);
                    }

                    let skipped = enumerate_placements(
                        pep,
                        mods,
                        prot_n,
                        prot_c,
                        MAX_VARIANTS_PER_PEPTIDE,
                        |placement, delta| {
                            out.push(Candidate {
                                mass: base + delta,
                                protein: pi as u32,
                                start: (protein.start + d.start) as u32,
                                len: d.len as u16,
                                n_termini: d.n_termini,
                                placement,
                            });
                        },
                    );
                    if skipped > 0 {
                        truncated += 1;
                    }
                });
                (out, truncated, n_pep)
            })
            .collect();

        let total: usize = per_protein.iter().map(|(c, _, _)| c.len()).sum();
        let mut candidates = Vec::with_capacity(total);
        let (mut truncated_peptides, mut n_peptides) = (0usize, 0usize);
        for (c, t, n) in per_protein {
            candidates.extend(c);
            truncated_peptides += t;
            n_peptides += n;
        }
        // Ascending mass: a precursor window is then a contiguous slice. Masses are finite, so the
        // partial comparison never sees NaN.
        candidates.par_sort_unstable_by(|a, b| a.mass.partial_cmp(&b.mass).expect("finite masses"));
        PeptideIndex {
            candidates,
            truncated_peptides,
            n_peptides,
        }
    }

    /// The candidates whose neutral mass falls in `[lo, hi]`.
    #[inline]
    pub fn window(&self, lo: f64, hi: f64) -> &[Candidate] {
        let a = self.candidates.partition_point(|c| c.mass < lo);
        let b = self.candidates.partition_point(|c| c.mass <= hi);
        &self.candidates[a..b]
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.candidates.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.candidates.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mods::ModSpec;
    use msgf_db::enzyme::Enzyme;
    use msgf_db::fasta::{DecoyStrategy, ProteinDb, DEFAULT_DECOY_PREFIX};
    use std::io::Write;

    fn db_from(body: &str, name: &str) -> ProteinDb {
        let dir = std::env::temp_dir().join("msgf-search-tests");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join(name);
        std::fs::File::create(&p)
            .unwrap()
            .write_all(body.as_bytes())
            .unwrap();
        ProteinDb::read(&p, DEFAULT_DECOY_PREFIX).unwrap()
    }

    fn params() -> DigestParams {
        DigestParams {
            enzyme: Enzyme::builtin(1).unwrap(),
            max_missed_cleavages: 0,
            min_len: 3,
            max_len: 30,
            min_termini: 2,
            cleave_initiator_met: false,
        }
    }

    #[test]
    fn candidate_masses_are_sorted_and_correct() {
        let db = db_from(">P1\nSAMPLERPEPTIDEK\n", "i1.fasta");
        let idx = PeptideIndex::build(&db, &params(), &ModSet::default());
        assert_eq!(idx.len(), 2); // SAMPLER, PEPTIDEK
        assert!(idx.candidates.windows(2).all(|w| w[0].mass <= w[1].mass));
        let want = msgf_chem::peptide_neutral_mass("SAMPLER").unwrap();
        assert!(idx.candidates.iter().any(|c| (c.mass - want).abs() < 1e-9));
    }

    #[test]
    fn fixed_mod_shifts_every_occurrence() {
        let db = db_from(">P1\nCCAGGK\n", "i2.fasta");
        let mods = ModSet {
            mods: vec![ModSpec::parse("C2H3N1O1,C,fix,any,Carbamidomethyl").unwrap()],
            max_var_mods: 0,
        };
        let idx = PeptideIndex::build(&db, &params(), &mods);
        assert_eq!(idx.len(), 1);
        let bare = msgf_chem::peptide_neutral_mass("CCAGGK").unwrap();
        assert!((idx.candidates[0].mass - bare - 2.0 * 57.021_463_7).abs() < 1e-5);
    }

    #[test]
    fn variable_mod_creates_extra_candidates() {
        let db = db_from(">P1\nMAMGGK\n", "i3.fasta");
        let mods = ModSet {
            mods: vec![ModSpec::parse("O1,M,opt,any,Oxidation").unwrap()],
            max_var_mods: 2,
        };
        let idx = PeptideIndex::build(&db, &params(), &mods);
        // unmodified, +1 on either M, +2 on both
        assert_eq!(idx.len(), 4);
        let bare = msgf_chem::peptide_neutral_mass("MAMGGK").unwrap();
        let heaviest = idx.candidates.last().unwrap().mass;
        assert!((heaviest - bare - 2.0 * 15.994_914_622).abs() < 1e-5);
    }

    #[test]
    fn window_selects_the_mass_slice() {
        let db = db_from(">P1\nSAMPLERPEPTIDEKAAAGGGKMMMK\n", "i4.fasta");
        let idx = PeptideIndex::build(&db, &params(), &ModSet::default());
        let target = msgf_chem::peptide_neutral_mass("SAMPLER").unwrap();
        let w = idx.window(target - 0.01, target + 0.01);
        assert_eq!(w.len(), 1);
        assert!((w[0].mass - target).abs() < 1e-9);
        assert!(idx.window(0.0, 1.0).is_empty());
    }

    #[test]
    fn decoys_double_the_index() {
        let mut db = db_from(">P1\nSAMPLERPEPTIDEK\n", "i5.fasta");
        db.add_decoys(DecoyStrategy::Reverse, DEFAULT_DECOY_PREFIX);
        let idx = PeptideIndex::build(&db, &params(), &ModSet::default());
        // The reversed protein KEDITPEPRELPMAS digests to two peptides as well.
        assert_eq!(idx.len(), 4);
        assert_eq!(
            idx.candidates
                .iter()
                .filter(|c| db.proteins[c.protein as usize].is_decoy)
                .count(),
            2
        );
    }
}
