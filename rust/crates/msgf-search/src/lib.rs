//! msgf-search — the MS-GF+ database search engine.
//!
//! ```text
//!   FASTA ──digest──► candidates ──precursor window──► RawScore ──► SpecEValue ──► q-values
//!  (msgf-db)          (index)                        (msgf-scorer) (msgf-genfunc)  (msgf-fdr)
//! ```
//!
//! A search is three steps:
//!
//! ```no_run
//! use msgf_db::{fasta::ProteinDb, enzyme::DigestParams, fasta::DEFAULT_DECOY_PREFIX};
//! use msgf_search::{index::PeptideIndex, mods::ModSet, SearchEngine, SearchParams, assign_q_values};
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let model = msgf_scorer::read_param_file("HCD_HighRes_Tryp.param")?;
//! let db = ProteinDb::read("human.revCat.fasta", DEFAULT_DECOY_PREFIX)?;
//! let (digest, mods) = (DigestParams::default(), ModSet::default());
//!
//! let index = PeptideIndex::build(&db, &digest, &mods);              // 1. build the candidate index
//! let engine = SearchEngine::new(&model, &db, &index, &mods, &digest, SearchParams::default());
//! let spectra = msgf_io::read_mgf_file("run.mgf")?;
//! let mut psms = engine.run(&spectra);                                // 2. search (parallel)
//! assign_q_values(&mut psms);                                         // 3. global target-decoy FDR
//! # Ok(()) }
//! ```
//!
//! Step 3 is deliberately separate and serial: FDR is a property of the whole result set, so it is
//! an epilogue to the parallel search (`plans/PLAN2.md` §TD-3).
//!
//! The generating function is built **once per `(spectrum, charge)`** and shared by every candidate
//! in the precursor window — see [`search`] for the details and for the two documented divergences
//! from MS-GF+ (N-terminal enzymes, E-value scaling).

pub mod index;
pub mod mods;
pub mod report;
pub mod search;

pub use index::{Candidate, PeptideIndex};
pub use mods::{ModSet, ModSpec};
pub use search::{Psm, SearchEngine, SearchParams, SearchScratch};

use msgf_fdr::{PsmRecord, TargetDecoyAnalysis};

/// Fill in `q_value` and `pep_q_value` on every PSM using target-decoy analysis.
///
/// FDR is global, so this must see the **whole** result set at once — call it after every spectrum
/// has been searched, never per spectrum. Decoy status is already decided per match from all of its
/// protein occurrences (`plans/PLAN2.md` §1.3).
///
/// A search over a database with no decoys yields q-values of 0; that is not an FDR estimate, and
/// [`has_decoys`] lets a caller detect and say so.
pub fn assign_q_values(psms: &mut [Psm]) {
    if psms.is_empty() {
        return;
    }
    let records: Vec<PsmRecord> = psms
        .iter()
        .map(|p| PsmRecord {
            score: p.spec_evalue as f32,
            peptide: p.peptide_key.clone(),
            is_decoy: p.is_decoy,
        })
        .collect();
    let tda = TargetDecoyAnalysis::new(&records, 1.0);
    for p in psms.iter_mut() {
        let score = p.spec_evalue as f32;
        p.q_value = tda.psm_q_value(score);
        p.pep_q_value = tda.pep_q_value(&p.peptide_key, score);
    }
}

/// Whether the result set contains any decoy match — i.e. whether the q-values mean anything.
pub fn has_decoys(psms: &[Psm]) -> bool {
    psms.iter().any(|p| p.is_decoy)
}

/// Number of **target** PSMs at or below a q-value threshold — the "IDs at 1% FDR" figure a search
/// is judged by.
pub fn n_targets_below(psms: &[Psm], q_threshold: f32) -> usize {
    psms.iter()
        .filter(|p| !p.is_decoy && p.q_value <= q_threshold)
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn psm(spec_evalue: f64, is_decoy: bool, peptide: &str) -> Psm {
        Psm {
            spec_index: 0,
            scan: String::new(),
            title: String::new(),
            precursor_mz: 0.0,
            charge: 2,
            isotope_error: 0,
            precursor_error_ppm: 0.0,
            peptide: peptide.to_string(),
            peptide_key: peptide.to_string(),
            proteins: vec![if is_decoy { "XXX_P".into() } else { "P".into() }],
            is_decoy,
            raw_score: 0,
            denovo_score: 0,
            spec_evalue,
            evalue: 0.0,
            q_value: f32::NAN,
            pep_q_value: f32::NAN,
        }
    }

    #[test]
    fn q_values_are_assigned_and_monotone() {
        let mut p = vec![
            psm(1e-10, false, "A"),
            psm(1e-9, false, "B"),
            psm(1e-8, false, "C"),
            psm(1e-7, false, "D"),
            psm(1e-6, true, "E"),
        ];
        assign_q_values(&mut p);
        assert!(p.iter().all(|x| x.q_value.is_finite()));
        assert_eq!(p[0].q_value, 0.0);
        assert!(has_decoys(&p));
        assert_eq!(n_targets_below(&p, 0.01), 4);
    }

    #[test]
    fn no_decoys_is_detectable() {
        let mut p = vec![psm(1e-10, false, "A"), psm(1e-9, false, "B")];
        assign_q_values(&mut p);
        assert!(!has_decoys(&p));
        assert!(p.iter().all(|x| x.q_value == 0.0));
    }

    #[test]
    fn empty_result_set_is_fine() {
        let mut p: Vec<Psm> = Vec::new();
        assign_q_values(&mut p);
        assert!(!has_decoys(&p));
    }
}
