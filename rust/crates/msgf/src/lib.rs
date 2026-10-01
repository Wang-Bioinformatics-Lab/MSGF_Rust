//! # msgf — MSGF_Rust as a single library
//!
//! A Rust implementation of **MS-GF+-style significance scoring** — the generating-function
//! spectral E-value (SpecEValue) for high-resolution tandem MS — plus a database search engine
//! built on top of it. The scorer and generating function are a clean-room implementation from
//! a written specification (`docs/cleanroom/` in the repository); no MS-GF+ code is used.
//!
//! This crate is a facade: it re-exports the workspace's `msgf-*` crates under short module names
//! so a downstream project takes **one** dependency instead of seven. Every item is the same type
//! as in the underlying crate, so mixing the two styles is fine.
//!
//! ```toml
//! [dependencies]
//! msgf = { git = "https://github.com/mwang87/MSGF_Rust" }
//! # scoring only, without the search engine and its rayon dependency:
//! msgf = { git = "https://github.com/mwang87/MSGF_Rust", default-features = false }
//! ```
//!
//! ## The two entry points
//!
//! **Rescore** — you already have peptide-spectrum matches and want MS-GF+ scores for them:
//!
//! ```no_run
//! use msgf::{genfunc, io, scorer};
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // The bundled MassIVE-KB-trained model; `scorer::ScoreModel::from_file(path)` loads another.
//! let model = scorer::bundled::score_model()?;
//! let spectra = io::read_mgf_file("run.mgf")?;
//! let spectrum = &spectra[0];
//! let peaks: Vec<(f64, f64)> = spectrum.peaks.iter().map(|p| (p.mz, p.intensity)).collect();
//! let ms2 = scorer::Ms2 {
//!     peaks: &peaks,
//!     precursor_mz: spectrum.precursor_mz.unwrap_or_default(),
//!     charge: spectrum.charge.unwrap_or(2),
//! };
//! let prep = scorer::prepare(&model, &ms2).expect("plausible precursor");
//!
//! // RawScore of one candidate (residues with modification deltas, N-terminal context decided
//! // by the caller), then the null distribution and its tail.
//! let residues: Vec<scorer::Residue> = b"SAMPLER"
//!     .iter()
//!     .map(|&letter| scorer::Residue { letter, delta: 0.0 })
//!     .collect();
//! let cand = scorer::Candidate { residues: &residues, n_term_credit: true };
//! let null = genfunc::NullModel::uniform(); // 20 residues at 0.05, trypsin, isotope range (0, 1)
//! let raw = scorer::match_and_terminal_score(&prep, &cand, &null.cleavage).unwrap();
//! let tail = genfunc::score_distribution(&prep, &null, Some(raw)).expect("reachable sink");
//! println!("RawScore {raw}  DeNovoScore {}  SpecEValue {:e}", tail.best_possible(), tail.tail_mass(raw));
//! # Ok(()) }
//! ```
//!
//! **Search** — you have a FASTA and want identifications with q-values (requires the default
//! `search` feature). See [`search`] for the full three-step example.
//!
//! ## Crate map
//!
//! | Module | Crate | What it holds |
//! |---|---|---|
//! | [`chem`] | `msgf-chem` | masses, residues, peptides, fragment ions, tolerance, mass-grid scaling |
//! | [`io`] | `msgf-io` | `Spectrum`/`Peak` and the MGF reader |
//! | [`scorer`] | `msgf-scorer` | `.param` model read/write, the bundled default model (`scorer::bundled`), spectrum preparation (`prepare`) → RawScore |
//! | [`genfunc`] | `msgf-genfunc` | null model + score-distribution DP → DeNovoScore / SpecEValue |
//! | [`db`] | `msgf-db` | FASTA, target-decoy construction, digestion |
//! | [`fdr`] | `msgf-fdr` | MS-GF+-compatible PSM- and peptide-level q-values |
//! | [`search`] | `msgf-search` | candidate index and the search driver |
//!
//! The scoring path is a linear chain — `io → scorer → genfunc` — with `chem` underneath all of it.
//! `docs/cleanroom/SPEC.md` in the repository specifies the scoring; `LICENSING.md` records its
//! clean-room provenance.

pub use msgf_chem as chem;
pub use msgf_genfunc as genfunc;
pub use msgf_io as io;
pub use msgf_scorer as scorer;

#[cfg(feature = "search")]
pub use msgf_db as db;
#[cfg(feature = "search")]
pub use msgf_fdr as fdr;
#[cfg(feature = "search")]
pub use msgf_search as search;

/// The workspace version, so a consumer can report which MSGF_Rust it linked.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// The most-used types, for `use msgf::prelude::*;`.
pub mod prelude {
    pub use crate::chem::{peptide::Residue, residue_mass, Tolerance, Unit};
    pub use crate::genfunc::{score_distribution, NullModel, NullTail};
    pub use crate::io::{MgfReader, Peak, Spectrum};
    pub use crate::scorer::{
        bundled, match_and_terminal_score, prepare, Candidate, Cleavage, Ms2, PreparedSpectrum,
        ScoreModel, ScoringModel,
    };

    #[cfg(feature = "search")]
    pub use crate::db::{enzyme::DigestParams, fasta::ProteinDb, Enzyme};
    #[cfg(feature = "search")]
    pub use crate::fdr::TargetDecoyAnalysis;
    #[cfg(feature = "search")]
    pub use crate::search::{PeptideIndex, Psm, SearchEngine, SearchParams};
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn re_exports_are_the_underlying_types() {
        // The facade must not wrap or shadow — these are type identities, checked at compile time.
        let _: fn(u8) -> Option<f64> = chem::residue_mass;
        let _: chem::Tolerance = chem::Tolerance::ppm(10.0);
        assert_eq!(genfunc::NullModel::uniform().alphabet.len(), 20);
        assert!(!VERSION.is_empty());
    }

    #[test]
    fn prelude_covers_the_pipeline() {
        use crate::prelude::*;
        let _: Spectrum = Spectrum::default();
        let _: NullModel = NullModel::uniform();
        let _: Tolerance = Tolerance::da(0.5);
    }

    #[cfg(feature = "search")]
    #[test]
    fn search_surface_is_reachable() {
        use crate::prelude::*;
        let p = DigestParams::default();
        assert_eq!(p.enzyme.name, "Tryp");
        let _: SearchParams = SearchParams::default();
        assert_eq!(fdr::peptide_key("K.SAMPLER.A"), "SAMPLER");
    }
}
