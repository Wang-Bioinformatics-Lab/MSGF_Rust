//! # msgf — MSGF_Rust as a single library
//!
//! A Rust reimplementation of **MS-GF+ significance scoring** — the generating-function spectral
//! E-value (SpecEValue) for high-resolution tandem MS — validated to be bit-exact against the
//! reference Java MS-GF+, plus a database search engine built on top of it.
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
//! use msgf::{chem, genfunc, io, scorer};
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // The bundled MassIVE-KB-trained model; `scorer::read_param_file(path)` loads another.
//! let model = scorer::bundled::model()?;
//! let spectra = io::read_mgf_file("run.mgf")?;
//! let spectrum = &spectra[0];
//!
//! let charge = spectrum.charge.unwrap_or(2);
//! let mz = spectrum.precursor_mz.unwrap_or_default();
//! let parent_mass = mz as f32 * charge as f32 - charge as f32 * chem::mass::PROTON as f32;
//! let peaks: Vec<(f32, f32)> = spectrum.peaks.iter().map(|p| (p.mz as f32, p.intensity as f32)).collect();
//!
//! let ranked = scorer::preprocess::preprocess(&model, charge, parent_mass, &peaks);
//! let scored = scorer::scored_spectrum::ScoredSpectrum::from_ranked_peaks(
//!     &model, charge, parent_mass, ranked);
//!
//! let residues = chem::peptide::parse("K.SAMPLER.A").expect("valid peptide");
//! let raw = scored.raw_score(
//!     &chem::peptide::nominal_prefix_masses(&residues),
//!     &chem::peptide::accurate_prefix_masses(&residues),
//!     chem::peptide::num_mods(&residues) as i32,
//! );
//! let _ = (raw, genfunc::graph::standard_aa());
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
//! | [`scorer`] | `msgf-scorer` | `.param` model read/write, the bundled default model (`scorer::bundled`), preprocessing, `ScoredSpectrum` → RawScore |
//! | [`genfunc`] | `msgf-genfunc` | de novo graph + score-distribution DP → DeNovoScore / SpecEValue |
//! | [`db`] | `msgf-db` | FASTA, target-decoy construction, digestion |
//! | [`fdr`] | `msgf-fdr` | MS-GF+-compatible PSM- and peptide-level q-values |
//! | [`search`] | `msgf-search` | candidate index and the search driver |
//!
//! The scoring path is a linear chain — `io → scorer → genfunc` — with `chem` underneath all of it.
//! `PLAN.md` in the repository has the design and `CLAUDE.md` the fidelity contract: integer scores
//! must match MS-GF+ **exactly**, SpecEValue within `|log10(rust/java)| ≤ 0.05`.

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
    pub use crate::genfunc::{compute, compute_into, DpScratch, GenFunc, ScoreDist};
    pub use crate::io::{MgfReader, Peak, Spectrum};
    pub use crate::scorer::{bundled, scored_spectrum::ScoredSpectrum, ScoringModel};

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
        assert_eq!(genfunc::graph::standard_aa().len(), 20);
        assert!(!VERSION.is_empty());
    }

    #[test]
    fn prelude_covers_the_pipeline() {
        use crate::prelude::*;
        let _: Spectrum = Spectrum::default();
        let _: DpScratch = DpScratch::default();
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
