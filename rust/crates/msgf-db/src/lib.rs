//! msgf-db — the sequence-database layer for MSGF_Rust.
//!
//! Three things a search needs before any spectrum is scored:
//!
//! - [`fasta`] — streaming FASTA reading into a flat [`fasta::ProteinDb`] (all sequences in one
//!   buffer, proteins as offset+length slices), plus the database's amino-acid composition, which
//!   is what weights the de novo graph's edges.
//! - [`decoy`] — MS-GF+-compatible target-decoy construction, validated **byte-for-byte** against
//!   the reference `.revCat.fasta` databases (see `PLAN2.md` §1.1).
//! - [`enzyme`] — enzymes and in-silico digestion, including semi- and non-enzymatic searches.
//!
//! Nothing here depends on a spectrum or a scoring model; the crate has one dependency
//! (`msgf-chem`, to know which residues are scorable).

pub mod decoy;
pub mod enzyme;
pub mod fasta;

pub use decoy::{write_decoy_database, DecoyOptions, LineSep, Output};
pub use enzyme::{digest, DigestParams, DigestedPeptide, Enzyme};
pub use fasta::{DecoyStrategy, Protein, ProteinDb};
