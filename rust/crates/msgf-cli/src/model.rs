//! Which scoring model a run uses.
//!
//! `--param` is optional: with no model given, every subcommand falls back to the model this
//! project ships (`msgf_scorer::bundled`) — HCD / high-resolution / tryptic, counted by
//! `msgf-train` from the CC0 MassIVE-KB corpus. That default exists so a plain `msgf search` works
//! out of the box **and** so nothing on the shipping path is UC-licensed; MS-GF+'s own `.param`
//! files remain usable via `--param`, they are just not distributed with this tool.

use std::path::Path;

use msgf_scorer::ScoringModel;

/// Load the model for a run: the file if one was given, otherwise the bundled default.
///
/// Returns the model and a human-readable description of where it came from — printed to stderr so
/// a result table can always be traced back to the model that produced it.
pub fn load(param: Option<&Path>) -> Result<(ScoringModel, String), String> {
    match param {
        Some(p) => {
            let m = msgf_scorer::read_param_file(p)
                .map_err(|e| format!("reading model {}: {e:?}", p.display()))?;
            Ok((m, p.display().to_string()))
        }
        None => {
            let m = msgf_scorer::bundled::model()
                .map_err(|e| format!("decoding the bundled model: {e:?}"))?;
            Ok((
                m,
                format!(
                    "{} (bundled default: {})",
                    msgf_scorer::bundled::NAME,
                    msgf_scorer::bundled::PROVENANCE
                ),
            ))
        }
    }
}

/// Report the model in use on stderr, so stdout stays a clean data stream.
pub fn announce(source: &str, model: &ScoringModel) {
    eprintln!(
        "model: {source} [{}/{}/{}/{}]",
        model.activation,
        model.instrument,
        model.enzyme.as_deref().unwrap_or("none"),
        model.protocol_name()
    );
}
