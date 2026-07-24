//! Result output. The column set mirrors MS-GF+'s TSV (`MzIDToTsv`) so a run can be diffed
//! against a Java one, and so downstream tooling that already parses MS-GF+ output keeps working.
//!
//! One row per reported match. `Protein` lists every occurrence of the peptide, semicolon-separated
//! — MS-GF+'s `-unroll 1` instead emits one row per occurrence, which [`write_tsv_unrolled`] does.

use crate::search::Psm;
use std::io::{self, Write};

/// The MS-GF+ TSV header.
pub const HEADER: &str = "#SpecFile\tSpecID\tScanNum\tFragMethod\tPrecursor\tIsotopeError\t\
PrecursorError(ppm)\tCharge\tPeptide\tProtein\tDeNovoScore\tMSGFScore\tSpecEValue\tEValue\t\
QValue\tPepQValue";

/// Write results as MS-GF+-compatible TSV, one row per match.
pub fn write_tsv(w: &mut impl Write, spec_file: &str, psms: &[Psm]) -> io::Result<()> {
    write_rows(w, spec_file, psms, false)
}

/// Write results with one row per protein occurrence (MS-GF+ `-unroll 1`).
pub fn write_tsv_unrolled(w: &mut impl Write, spec_file: &str, psms: &[Psm]) -> io::Result<()> {
    write_rows(w, spec_file, psms, true)
}

fn write_rows(w: &mut impl Write, spec_file: &str, psms: &[Psm], unroll: bool) -> io::Result<()> {
    writeln!(w, "{HEADER}")?;
    for p in psms {
        let joined = p.proteins.join(";");
        let cells: Vec<&str> = if unroll && !p.proteins.is_empty() {
            p.proteins.iter().map(String::as_str).collect()
        } else {
            vec![joined.as_str()]
        };
        for protein in cells {
            writeln!(
                w,
                "{spec_file}\tindex={}\t{}\t{}\t{:.5}\t{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{:.6E}\t{:.6E}\t{}\t{}",
                p.spec_index,
                p.scan,
                fragmentation_method(&p.title),
                p.precursor_mz,
                p.isotope_error,
                p.precursor_error_ppm,
                p.charge,
                p.peptide,
                protein,
                p.denovo_score,
                p.raw_score,
                p.spec_evalue,
                p.evalue,
                fmt_q(p.q_value),
                fmt_q(p.pep_q_value),
            )?;
        }
    }
    Ok(())
}

/// Q-values are `NaN` until [`crate::assign_q_values`] runs; render that as `NA` rather than
/// printing `NaN` into a numeric column.
fn fmt_q(q: f32) -> String {
    if q.is_nan() {
        "NA".to_string()
    } else {
        format!("{q}")
    }
}

/// MS-GF+ records the activation method per spectrum; MGF rarely carries it, so this reports what
/// the title says when it says anything and `N/A` otherwise.
fn fragmentation_method(title: &str) -> &str {
    for m in ["HCD", "CID", "ETD", "ETHCD", "UVPD"] {
        if title.to_ascii_uppercase().contains(m) {
            return m;
        }
    }
    "N/A"
}

/// A one-line summary of a finished search, for stderr.
pub fn summary(psms: &[Psm], has_decoys: bool) -> String {
    let n = psms.len();
    let targets = psms.iter().filter(|p| !p.is_decoy).count();
    if !has_decoys {
        return format!(
            "{n} PSM(s); no decoys in the database, so q-values are not an FDR estimate"
        );
    }
    let at_1 = crate::n_targets_below(psms, 0.01);
    let at_5 = crate::n_targets_below(psms, 0.05);
    format!("{n} PSM(s), {targets} target; {at_1} target PSMs at 1% FDR, {at_5} at 5%")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn psm(peptide: &str, proteins: &[&str]) -> Psm {
        Psm {
            spec_index: 0,
            scan: "42".into(),
            title: "scan=42 HCD".into(),
            precursor_mz: 500.25,
            charge: 2,
            isotope_error: 0,
            precursor_error_ppm: 1.5,
            peptide: peptide.into(),
            peptide_key: peptide.into(),
            proteins: proteins.iter().map(|s| s.to_string()).collect(),
            is_decoy: false,
            raw_score: 30,
            denovo_score: 55,
            spec_evalue: 1e-10,
            evalue: 1e-4,
            q_value: 0.0,
            pep_q_value: 0.0,
        }
    }

    #[test]
    fn header_matches_msgfplus_columns() {
        assert!(HEADER.starts_with("#SpecFile\tSpecID\tScanNum"));
        assert!(HEADER.ends_with("QValue\tPepQValue"));
        assert_eq!(HEADER.split('\t').count(), 16);
    }

    #[test]
    fn one_row_per_match_joins_proteins() {
        let mut out = Vec::new();
        write_tsv(&mut out, "run.mgf", &[psm("K.SAMPLER.A", &["P1", "P2"])]).unwrap();
        let s = String::from_utf8(out).unwrap();
        let rows: Vec<&str> = s.lines().skip(1).collect();
        assert_eq!(rows.len(), 1);
        assert!(rows[0].contains("P1;P2"), "{}", rows[0]);
        assert_eq!(rows[0].split('\t').count(), 16);
    }

    #[test]
    fn unrolled_emits_one_row_per_protein() {
        let mut out = Vec::new();
        write_tsv_unrolled(&mut out, "run.mgf", &[psm("K.SAMPLER.A", &["P1", "P2"])]).unwrap();
        let s = String::from_utf8(out).unwrap();
        let rows: Vec<&str> = s.lines().skip(1).collect();
        assert_eq!(rows.len(), 2);
        assert!(rows[0].contains("\tP1\t"));
        assert!(rows[1].contains("\tP2\t"));
    }

    #[test]
    fn unassigned_q_values_render_as_na() {
        let mut p = psm("K.SAMPLER.A", &["P1"]);
        p.q_value = f32::NAN;
        p.pep_q_value = f32::NAN;
        let mut out = Vec::new();
        write_tsv(&mut out, "run.mgf", &[p]).unwrap();
        let s = String::from_utf8(out).unwrap();
        assert!(s.lines().nth(1).unwrap().ends_with("NA\tNA"));
    }

    #[test]
    fn fragmentation_method_from_title() {
        assert_eq!(fragmentation_method("scan=1 HCD"), "HCD");
        assert_eq!(fragmentation_method("nothing here"), "N/A");
    }
}
