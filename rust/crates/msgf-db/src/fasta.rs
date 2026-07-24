//! FASTA reading, the in-memory protein database, and target-decoy generation.
//!
//! Sequences are concatenated into one `Vec<u8>` and each protein records an `(start, len)` slice
//! into it, so a digested peptide is addressed by `(protein index, offset, length)` with no
//! per-peptide string allocation. Residues are upper-cased on read; non-standard characters
//! (`X`, `B`, `Z`, `U`, `O`, `J`, `*`) are kept in the buffer but peptides spanning them are
//! rejected at digestion time (their masses are undefined — [`msgf_chem::residue_mass`] is `None`).

use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// MS-GF+'s default decoy accession prefix (`-tda 1` writes `XXX_<accession>`). Already
/// normalised — see [`crate::decoy::normalize_prefix`] for the `XXX`/`XXX_` equivalence.
pub const DEFAULT_DECOY_PREFIX: &str = "XXX_";

/// How to build the decoy half of a target-decoy database.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DecoyStrategy {
    /// Reverse each target protein sequence (MS-GF+ `-tda 1`).
    Reverse,
    /// Do not generate decoys — either the FASTA already contains them (e.g. a `.revCat.fasta`)
    /// or the caller does not want an FDR estimate.
    None,
}

/// One protein entry: its accession, description, and the slice of [`ProteinDb::seq`] it owns.
#[derive(Debug, Clone)]
pub struct Protein {
    /// Accession — the first whitespace-delimited token of the FASTA header.
    pub name: String,
    /// The remainder of the header line.
    pub desc: String,
    pub start: usize,
    pub len: usize,
    pub is_decoy: bool,
}

/// The searched protein database: all sequences in one buffer plus the per-protein index.
#[derive(Debug, Clone, Default)]
pub struct ProteinDb {
    /// Every protein sequence concatenated, upper-cased.
    pub seq: Vec<u8>,
    pub proteins: Vec<Protein>,
}

impl ProteinDb {
    /// Parse a FASTA file. Blank lines and lines before the first `>` header are ignored;
    /// whitespace inside sequences is stripped.
    ///
    /// A protein whose accession already starts with `decoy_prefix` is flagged as a decoy on read,
    /// which is how a pre-built concatenated database (`*.revCat.fasta`) is searched without
    /// generating a second set of decoys.
    pub fn read<P: AsRef<Path>>(path: P, decoy_prefix: &str) -> io::Result<ProteinDb> {
        let file = File::open(path)?;
        let mut db = ProteinDb::default();
        let mut cur: Option<Protein> = None;

        for line in BufReader::new(file).lines() {
            let line = line?;
            let line = line.trim_end();
            if let Some(header) = line.strip_prefix('>') {
                if let Some(mut p) = cur.take() {
                    p.len = db.seq.len() - p.start;
                    db.proteins.push(p);
                }
                let header = header.trim_start();
                let (name, desc) = match header.find(char::is_whitespace) {
                    Some(i) => (&header[..i], header[i..].trim_start()),
                    None => (header, ""),
                };
                cur = Some(Protein {
                    name: name.to_string(),
                    desc: desc.to_string(),
                    start: db.seq.len(),
                    len: 0,
                    is_decoy: !decoy_prefix.is_empty() && name.starts_with(decoy_prefix),
                });
            } else if cur.is_some() {
                db.seq.extend(
                    line.bytes()
                        .filter(|b| !b.is_ascii_whitespace())
                        .map(|b| b.to_ascii_uppercase()),
                );
            }
        }
        if let Some(mut p) = cur.take() {
            p.len = db.seq.len() - p.start;
            db.proteins.push(p);
        }
        Ok(db)
    }

    /// The residue slice of one protein.
    #[inline]
    pub fn sequence(&self, i: usize) -> &[u8] {
        let p = &self.proteins[i];
        &self.seq[p.start..p.start + p.len]
    }

    /// Number of proteins already flagged as decoys.
    pub fn n_decoys(&self) -> usize {
        self.proteins.iter().filter(|p| p.is_decoy).count()
    }

    /// Append a decoy for every target protein. No-op for [`DecoyStrategy::None`], and no-op if the
    /// database already carries decoys (a pre-concatenated `revCat` FASTA) — re-reversing those
    /// would produce duplicate targets and corrupt the FDR estimate.
    pub fn add_decoys(&mut self, strategy: DecoyStrategy, prefix: &str) -> usize {
        if strategy == DecoyStrategy::None || self.n_decoys() > 0 {
            return 0;
        }
        let n_target = self.proteins.len();
        for i in 0..n_target {
            let (start, len) = (self.proteins[i].start, self.proteins[i].len);
            let decoy_start = self.seq.len();
            // Reverse the residues (MS-GF+ `-tda 1`): same composition and length, so the decoy
            // peptide mass distribution matches the target one.
            for k in (0..len).rev() {
                self.seq.push(self.seq[start + k]);
            }
            self.proteins.push(Protein {
                name: format!("{prefix}{}", self.proteins[i].name),
                desc: self.proteins[i].desc.clone(),
                start: decoy_start,
                len,
                is_decoy: true,
            });
        }
        n_target
    }

    /// Background amino-acid frequencies over the **target** proteins, counting only the 20
    /// standard residues (mirrors `DBScanner.setAminoAcidProbabilities`, which normalises over the
    /// residues it can score). These weight the de novo graph's edges, so the generating function
    /// reflects the composition of the database actually being searched rather than a uniform 1/20.
    ///
    /// Falls back to uniform `0.05` if the database has no scorable residue.
    pub fn aa_probabilities(&self) -> Vec<(u8, f64)> {
        let mut counts = [0u64; 26];
        let mut total = 0u64;
        for (i, p) in self.proteins.iter().enumerate() {
            if p.is_decoy {
                continue; // decoys are permutations of targets — same composition, don't double count
            }
            for &b in self.sequence(i) {
                if msgf_chem::residue_mass(b).is_some() {
                    counts[(b - b'A') as usize] += 1;
                    total += 1;
                }
            }
        }
        let residues: Vec<u8> = (b'A'..=b'Z')
            .filter(|&b| msgf_chem::residue_mass(b).is_some())
            .collect();
        if total == 0 {
            return residues.into_iter().map(|r| (r, 0.05)).collect();
        }
        residues
            .into_iter()
            .map(|r| (r, counts[(r - b'A') as usize] as f64 / total as f64))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_tmp(name: &str, body: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join("msgf-db-tests");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join(name);
        let mut f = File::create(&p).unwrap();
        f.write_all(body.as_bytes()).unwrap();
        p
    }

    #[test]
    fn reads_multiline_fasta() {
        let p = write_tmp(
            "a.fasta",
            ">sp|P1|ONE first protein\nPEPTIDE\nKSAMPLER\n>sp|P2|TWO\nMKRAA\n",
        );
        let db = ProteinDb::read(&p, DEFAULT_DECOY_PREFIX).unwrap();
        assert_eq!(db.proteins.len(), 2);
        assert_eq!(db.proteins[0].name, "sp|P1|ONE");
        assert_eq!(db.proteins[0].desc, "first protein");
        assert_eq!(db.sequence(0), b"PEPTIDEKSAMPLER");
        assert_eq!(db.sequence(1), b"MKRAA");
    }

    #[test]
    fn decoys_are_reversed_and_prefixed() {
        let p = write_tmp("b.fasta", ">P1\nPEPTIDEK\n");
        let mut db = ProteinDb::read(&p, DEFAULT_DECOY_PREFIX).unwrap();
        assert_eq!(
            db.add_decoys(DecoyStrategy::Reverse, DEFAULT_DECOY_PREFIX),
            1
        );
        assert_eq!(db.proteins.len(), 2);
        assert_eq!(db.proteins[1].name, "XXX_P1");
        assert!(db.proteins[1].is_decoy);
        assert_eq!(db.sequence(1), b"KEDITPEP");
    }

    #[test]
    fn existing_decoys_are_detected_not_regenerated() {
        let p = write_tmp("c.fasta", ">P1\nPEPTIDEK\n>XXX_P1\nKEDITPEP\n");
        let mut db = ProteinDb::read(&p, DEFAULT_DECOY_PREFIX).unwrap();
        assert_eq!(db.n_decoys(), 1);
        assert_eq!(
            db.add_decoys(DecoyStrategy::Reverse, DEFAULT_DECOY_PREFIX),
            0
        );
        assert_eq!(db.proteins.len(), 2);
    }

    #[test]
    fn composition_ignores_nonstandard_and_decoys() {
        let p = write_tmp("d.fasta", ">P1\nAAAAG\n>XXX_P1\nGAAAA\n");
        let db = ProteinDb::read(&p, DEFAULT_DECOY_PREFIX).unwrap();
        let probs: std::collections::HashMap<u8, f64> = db.aa_probabilities().into_iter().collect();
        assert!((probs[&b'A'] - 0.8).abs() < 1e-12);
        assert!((probs[&b'G'] - 0.2).abs() < 1e-12);
        assert!((probs[&b'W'] - 0.0).abs() < 1e-12);
    }
}
