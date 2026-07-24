//! msgf-train — a clean-room trainer that **produces** a fragment-scoring model (`.param`).
//!
//! This is the last piece of the model-ownership plan (`PLAN1.md` step 5): `msgf-scorer` can read
//! and write the format, and this crate fills it with numbers counted from a corpus of confident
//! peptide-spectrum matches, so MSGF_Rust can ship a model of its own rather than a UC-licensed one.
//!
//! **Clean-room boundary.** Nothing here is transcribed from MS-GF+'s
//! `ScoringParameterGeneratorWithErrors`. The container is the documented format
//! (`docs/param-format.md`); the *statistics* are defined in [`counts`] from the semantics the
//! scorer gives each table (`ScoringModel::score_from_table` consumes `ln(ion/noise)`, so the
//! trainer produces exactly that ratio's numerator and denominator). Constants like the mass
//! tolerance, segment count and rank ceiling are configuration of the identity being trained
//! ([`TrainConfig`]), not trained values.
//!
//! Training is a counting pass — no optimiser, no randomness — so the same corpus and config
//! always produce the same bytes.
//!
//! ```no_run
//! # use msgf_train::{TrainConfig, corpus, counts};
//! let cfg = TrainConfig::high_res_hcd_tryptic();
//! let (mut psms, mut stats) = (Vec::new(), corpus::CorpusStats::default());
//! corpus::read_annotated_mgf("library.mgf", &corpus::CorpusFilter::default(), &mut psms, &mut stats).unwrap();
//! let (model, _report, _scheme, _n) = counts::train(&psms, &cfg);
//! msgf_scorer::write_param_file("HCD_HighRes_Tryp.param", &model).unwrap();
//! ```

pub mod corpus;
pub mod counts;
pub mod ions;
pub mod partition;

use msgf_chem::Tolerance;

/// Everything that is *chosen* rather than counted: the model identity, the instrument
/// configuration the statistics are conditioned on, and the counting knobs.
#[derive(Debug, Clone)]
pub struct TrainConfig {
    // ---- identity (written into the model header)
    pub version: i32,
    pub activation: String,
    pub instrument: String,
    pub enzyme: Option<String>,
    /// `None` is MS-GF+'s "Automatic".
    pub protocol: Option<String>,

    // ---- instrument configuration the tables are conditioned on
    /// Fragment mass tolerance used for peak matching, and stored in the model.
    pub mme: Tolerance,
    pub apply_deconvolution: bool,
    pub deconvolution_error_tolerance: f32,
    /// How many m/z segments each precursor is split into (ion statistics differ low vs. high m/z).
    pub num_segments: i32,
    /// Highest intensity rank scored; ranks beyond it reuse the last bin, and bin `max_rank` is
    /// "ion absent".
    pub max_rank: i32,
    /// Mass-error histogram resolution: `2·esf + 1` bins over ±1 Da.
    pub error_scaling_factor: i32,
    /// Highest fragment charge considered as a candidate ion type.
    pub max_fragment_charge: i32,

    // ---- partition scheme
    pub charge_min: i32,
    pub charge_max: i32,
    pub min_psms_per_partition: usize,
    pub max_partitions_per_charge: usize,

    // ---- selection thresholds
    /// An ion type is scored in a partition when it is observed at ≥ this fraction of sites.
    pub ion_freq_threshold: f32,
    pub max_ions_per_partition: usize,
    pub precursor_freq_threshold: f32,
    /// How far a precursor offset must stand above the median frequency of the scanned window to
    /// count as a real precursor artefact rather than ordinary fragment density.
    pub precursor_contrast: f32,
    /// Cap on precursor offsets kept per (charge, reduced charge).
    pub max_precursor_offsets_per_charge: usize,
    /// Emit chemistry-derived precursor offsets for any charge where counting found none — the
    /// case when the corpus is library spectra whose precursor region was stripped before deposit.
    pub precursor_defaults: bool,
    /// Nominal-grid range scanned for precursor offsets.
    pub precursor_offset_lo: i32,
    pub precursor_offset_hi: i32,

    // ---- smoothing (keeps every ln(ion/noise) finite)
    /// Add-λ on rank-distribution bins, in units of *sites per partition*.
    pub smoothing: f64,
    /// Rank-pooling width as a fraction of the rank (0 disables). High ranks are observed rarely,
    /// and a ratio of two sparse bins is noise, not signal.
    pub rank_smoothing: f64,
    /// Add-λ on error-distribution bins.
    pub error_smoothing: f64,
}

impl TrainConfig {
    /// The identity this project validates against first: HCD, high-resolution, tryptic.
    ///
    /// The instrument settings match what a high-resolution `.param` is built for — 0.5 Da
    /// fragment matching window, isotope deconvolution at 0.02 Da, two mass segments, rank ceiling
    /// 150, mass-error histogram at 0.01 Da resolution.
    pub fn high_res_hcd_tryptic() -> Self {
        Self {
            version: 1,
            activation: "HCD".into(),
            instrument: "HighRes".into(),
            enzyme: Some("Tryp".into()),
            protocol: None,
            mme: Tolerance::da(0.5),
            apply_deconvolution: true,
            deconvolution_error_tolerance: 0.02,
            num_segments: 2,
            max_rank: 150,
            error_scaling_factor: 100,
            max_fragment_charge: 2,
            charge_min: 2,
            charge_max: 8,
            min_psms_per_partition: 400,
            max_partitions_per_charge: 30,
            ion_freq_threshold: 0.15,
            max_ions_per_partition: 6,
            precursor_freq_threshold: 0.15,
            precursor_contrast: 2.0,
            max_precursor_offsets_per_charge: 8,
            precursor_defaults: true,
            precursor_offset_lo: -70,
            precursor_offset_hi: 10,
            smoothing: 0.005,
            rank_smoothing: 0.1,
            error_smoothing: 0.01,
        }
    }
}
