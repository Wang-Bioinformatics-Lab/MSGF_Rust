//! Target-decoy FASTA construction, byte-compatible with MS-GF+'s `-tda 1`
//! (`msdbsearch/ReverseDB.java`; see `plans/PLAN2.md` §1.1).
//!
//! The output is a **concatenated** database: the whole target block, then one decoy per target.
//! The byte-level rules below were confirmed by reconstructing both reference databases shipped
//! with MS-GF+ (`Tryp_Pig_Bov.revCat.fasta`, `iprg2013_human.revCat.fasta`) byte-for-byte:
//!
//! 1. Input is consumed **line by line**, stripping `\r\n` or `\n`. Every line of the target block
//!    is re-emitted with the writer's line separator — so the target block's *wrapping* is
//!    preserved but its line endings are normalised. (This is why the two reference files differ:
//!    the pig/bovine one was generated on Windows and kept CRLF, the human one on Linux and did
//!    not. [`DecoyOptions::line_sep`] selects which.)
//! 2. A decoy header is `>` + prefix + the original header text (the **whole** header, not just
//!    the accession).
//! 3. A decoy sequence is every sequence line of the record concatenated **verbatim** — internal
//!    whitespace included — then reversed as a whole, then trimmed, then written as **one
//!    unwrapped line**. (The `iprg2013_human` GFP record carries trailing spaces mid-sequence and
//!    is the record that pins this rule.)
//!
//! The prefix follows MS-GF+'s normalisation: any trailing `_` is stripped and exactly one is
//! re-appended, so `XXX`, `XXX_`, and `XXX__` all yield `XXX_`.

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// MS-GF+'s default decoy prefix.
pub const DEFAULT_PREFIX: &str = "XXX";

/// Line separator for the generated file.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum LineSep {
    /// `\n` — what a Linux-generated MS-GF+ database uses (`iprg2013_human.revCat.fasta`).
    #[default]
    Lf,
    /// `\r\n` — what a Windows-generated one uses (`Tryp_Pig_Bov.revCat.fasta`).
    Crlf,
}

impl LineSep {
    #[inline]
    fn bytes(self) -> &'static [u8] {
        match self {
            LineSep::Lf => b"\n",
            LineSep::Crlf => b"\r\n",
        }
    }
}

/// What to write.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Output {
    /// Targets followed by decoys — what MS-GF+ `-tda 1` searches.
    Concatenated,
    /// Decoys only.
    DecoyOnly,
}

/// Options for [`write_decoy_database`].
#[derive(Debug, Clone)]
pub struct DecoyOptions {
    /// Accession prefix; normalised to end in exactly one `_`.
    pub prefix: String,
    pub line_sep: LineSep,
    pub output: Output,
}

impl Default for DecoyOptions {
    fn default() -> Self {
        Self {
            prefix: DEFAULT_PREFIX.to_string(),
            line_sep: LineSep::default(),
            output: Output::Concatenated,
        }
    }
}

/// Normalise a decoy prefix the way MS-GF+ does: strip trailing `_`, then append exactly one.
pub fn normalize_prefix(prefix: &str) -> String {
    format!("{}_", prefix.trim_end_matches('_'))
}

/// Write a target-decoy database from `input` to `output`. Returns the number of decoys written.
///
/// Streams the input: memory use is one protein record, not the whole database.
pub fn write_decoy_database(input: &Path, output: &Path, opts: &DecoyOptions) -> io::Result<usize> {
    let sep = opts.line_sep.bytes();
    let prefix = normalize_prefix(&opts.prefix);
    let mut reader = BufReader::new(File::open(input)?);
    let mut writer = BufWriter::new(File::create(output)?);

    // Pass 1: echo the target block (line-by-line) while collecting each record's header and its
    // raw concatenated sequence, so the decoy block can be appended after it.
    let mut headers: Vec<String> = Vec::new();
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    let mut line = Vec::new();
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        // Strip the terminator (`\r\n` or `\n`); everything else is content.
        while line.last() == Some(&b'\n') || line.last() == Some(&b'\r') {
            line.pop();
        }
        if opts.output == Output::Concatenated {
            writer.write_all(&line)?;
            writer.write_all(sep)?;
        }
        if line.first() == Some(&b'>') {
            headers.push(String::from_utf8_lossy(&line[1..]).into_owned());
            seqs.push(Vec::new());
        } else if let Some(cur) = seqs.last_mut() {
            cur.extend_from_slice(&line); // verbatim — internal whitespace is part of the sequence
        }
    }

    // Pass 2: the decoy block.
    for (hdr, seq) in headers.iter().zip(&seqs) {
        writer.write_all(b">")?;
        writer.write_all(prefix.as_bytes())?;
        writer.write_all(hdr.as_bytes())?;
        writer.write_all(sep)?;
        let mut rev: Vec<u8> = seq.iter().rev().copied().collect();
        trim_ascii_whitespace(&mut rev);
        writer.write_all(&rev)?;
        writer.write_all(sep)?;
    }
    writer.flush()?;
    Ok(headers.len())
}

/// In-place equivalent of Java's `String.trim()` for a byte buffer.
fn trim_ascii_whitespace(v: &mut Vec<u8>) {
    while v.last().is_some_and(|b| b.is_ascii_whitespace()) {
        v.pop();
    }
    let lead = v.iter().take_while(|b| b.is_ascii_whitespace()).count();
    v.drain(..lead);
}

/// The two load-time sanity checks MS-GF+ applies to a concatenated database
/// (`MSGFPlus.java:238-252`). Returns an error message when the database looks unusable.
///
/// - the fraction of **unique** protein accessions must be at least 0.5;
/// - the decoy fraction must land in `[0.4, 0.6]`.
pub fn validate_concatenated(accessions: &[String], decoy_prefix: &str) -> Result<(), String> {
    if accessions.is_empty() {
        return Err("database contains no proteins".into());
    }
    let n = accessions.len() as f64;
    let unique: std::collections::HashSet<&str> = accessions.iter().map(String::as_str).collect();
    let unique_ratio = unique.len() as f64 / n;
    if unique_ratio < 0.5 {
        return Err(format!(
            "only {:.1}% of protein accessions are unique (MS-GF+ requires >= 50%); \
             the database is probably duplicated",
            unique_ratio * 100.0
        ));
    }
    let prefix = normalize_prefix(decoy_prefix);
    let decoy_fraction = accessions.iter().filter(|a| a.starts_with(&prefix)).count() as f64 / n;
    if !(0.4..=0.6).contains(&decoy_fraction) {
        return Err(format!(
            "decoys are {:.1}% of the database, outside the required 40-60% \
             (expected a concatenated target-decoy FASTA with prefix `{prefix}`)",
            decoy_fraction * 100.0
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    fn tmp(name: &str) -> std::path::PathBuf {
        // Unit tests do not get CARGO_TARGET_TMPDIR (integration tests do), so use the system
        // temp dir with a crate-scoped subdirectory.
        let dir = std::env::temp_dir().join("msgf-db-tests");
        std::fs::create_dir_all(&dir).unwrap();
        dir.join(name)
    }

    fn run(body: &[u8], opts: &DecoyOptions, name: &str) -> Vec<u8> {
        let inp = tmp(&format!("{name}.in.fasta"));
        let out = tmp(&format!("{name}.out.fasta"));
        File::create(&inp).unwrap().write_all(body).unwrap();
        write_decoy_database(&inp, &out, opts).unwrap();
        let mut v = Vec::new();
        File::open(&out).unwrap().read_to_end(&mut v).unwrap();
        v
    }

    #[test]
    fn prefix_is_normalized_like_msgfplus() {
        assert_eq!(normalize_prefix("XXX"), "XXX_");
        assert_eq!(normalize_prefix("XXX_"), "XXX_");
        assert_eq!(normalize_prefix("REV__"), "REV_");
    }

    #[test]
    fn concatenates_target_then_reversed_decoy() {
        let got = run(b">P1 desc\nPEPTIDEK\n", &DecoyOptions::default(), "d1");
        assert_eq!(got, b">P1 desc\nPEPTIDEK\n>XXX_P1 desc\nKEDITPEP\n");
    }

    #[test]
    fn wrapping_is_preserved_in_the_target_block_but_decoys_are_one_line() {
        let got = run(b">P1\nPEPT\nIDEK\n", &DecoyOptions::default(), "d2");
        assert_eq!(got, b">P1\nPEPT\nIDEK\n>XXX_P1\nKEDITPEP\n");
    }

    #[test]
    fn crlf_input_is_normalized_to_the_chosen_separator() {
        let got = run(b">P1\r\nPEPTIDEK\r\n", &DecoyOptions::default(), "d3");
        assert_eq!(got, b">P1\nPEPTIDEK\n>XXX_P1\nKEDITPEP\n");
        let opts = DecoyOptions {
            line_sep: LineSep::Crlf,
            ..Default::default()
        };
        let got = run(b">P1\r\nPEPTIDEK\r\n", &opts, "d4");
        assert_eq!(got, b">P1\r\nPEPTIDEK\r\n>XXX_P1\r\nKEDITPEP\r\n");
    }

    #[test]
    fn internal_whitespace_survives_reversal_then_is_trimmed() {
        // The rule the iprg2013_human GFP record pins: sequence lines concatenate verbatim, the
        // whole string reverses, and only the ends are trimmed.
        let got = run(b">P1\nPEPT \nIDEK\n", &DecoyOptions::default(), "d5");
        assert_eq!(got, b">P1\nPEPT \nIDEK\n>XXX_P1\nKEDI TPEP\n");
        // A trailing space on the final line becomes a leading one after reversal, and is trimmed.
        let got = run(b">P1\nPEPTIDEK \n", &DecoyOptions::default(), "d6");
        assert_eq!(got, b">P1\nPEPTIDEK \n>XXX_P1\nKEDITPEP\n");
    }

    #[test]
    fn decoy_only_output_skips_the_target_block() {
        let opts = DecoyOptions {
            output: Output::DecoyOnly,
            ..Default::default()
        };
        let got = run(b">P1\nPEPTIDEK\n", &opts, "d7");
        assert_eq!(got, b">XXX_P1\nKEDITPEP\n");
    }

    #[test]
    fn validation_catches_bad_databases() {
        let ok: Vec<String> = vec!["A".into(), "B".into(), "XXX_A".into(), "XXX_B".into()];
        assert!(validate_concatenated(&ok, "XXX").is_ok());
        // no decoys at all
        assert!(validate_concatenated(&["A".to_string(), "B".to_string()], "XXX").is_err());
        // A 50% unique ratio is exactly on MS-GF+'s boundary and is accepted.
        let boundary: Vec<String> = vec!["A".into(), "A".into(), "XXX_A".into(), "XXX_A".into()];
        assert!(validate_concatenated(&boundary, "XXX").is_ok());
        // Heavily duplicated accessions (25% unique) are rejected, with the decoy fraction still 50%.
        let dup: Vec<String> = vec![
            "A".into(),
            "A".into(),
            "A".into(),
            "A".into(),
            "XXX_A".into(),
            "XXX_A".into(),
            "XXX_A".into(),
            "XXX_A".into(),
        ];
        assert!(validate_concatenated(&dup, "XXX").is_err());
        assert!(validate_concatenated(&[], "XXX").is_err());
    }
}
