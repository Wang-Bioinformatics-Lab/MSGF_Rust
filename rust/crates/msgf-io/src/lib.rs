//! msgf-io — spectrum types and format readers for MSGF_Rust.
//!
//! This first cut implements a streaming MGF reader. mzML (via the `mzdata` crate) comes next.
//! Parsing is validated byte-for-byte against `validation/golden/spectra/` in
//! `tests/golden_spectra.rs`, including a canonical peak-list hash.

use std::io::{self, BufRead};
use std::path::Path;

/// A single m/z + intensity pair.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Peak {
    pub mz: f64,
    pub intensity: f64,
}

/// One MS/MS spectrum.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Spectrum {
    /// 0-based position in the source file.
    pub index: usize,
    pub title: Option<String>,
    pub scan: Option<String>,
    pub charge: Option<i32>,
    pub precursor_mz: Option<f64>,
    pub peaks: Vec<Peak>,
}

impl Spectrum {
    pub fn n_peaks(&self) -> usize {
        self.peaks.len()
    }

    /// Canonical textual representation of the peak list: one line `"{mz:.5} {intensity:.5}"`
    /// per peak in file order, joined by `\n`. This is the exact form hashed by the golden
    /// (`peaks_sha1`), so any reader that parses identically reproduces the hash.
    pub fn canonical_peak_string(&self) -> String {
        let mut s = String::with_capacity(self.peaks.len() * 20);
        for (i, p) in self.peaks.iter().enumerate() {
            if i > 0 {
                s.push('\n');
            }
            s.push_str(&format!("{:.5} {:.5}", p.mz, p.intensity));
        }
        s
    }

    /// Total ion current (Σ intensities).
    pub fn tic(&self) -> f64 {
        self.peaks.iter().map(|p| p.intensity).sum()
    }
}

/// Streaming MGF reader over any buffered source.
pub struct MgfReader<R: BufRead> {
    reader: R,
    index: usize,
    line: String,
}

impl<R: BufRead> MgfReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            index: 0,
            line: String::new(),
        }
    }

    fn read_line(&mut self) -> io::Result<Option<&str>> {
        self.line.clear();
        let n = self.reader.read_line(&mut self.line)?;
        if n == 0 {
            Ok(None)
        } else {
            Ok(Some(self.line.trim_end_matches(['\r', '\n'])))
        }
    }

    fn read_next(&mut self) -> io::Result<Option<Spectrum>> {
        // scan forward to BEGIN IONS
        loop {
            match self.read_line()? {
                None => return Ok(None),
                Some(l) if l.trim() == "BEGIN IONS" => break,
                _ => continue,
            }
        }
        let mut spec = Spectrum {
            index: self.index,
            ..Default::default()
        };
        loop {
            let line = match self.read_line()? {
                None => break, // EOF inside a spectrum: emit what we have
                Some(l) => l,
            };
            if !apply_body_line(&mut spec, line) {
                break;
            }
        }
        self.index += 1;
        Ok(Some(spec))
    }
}

/// Apply one line (line terminator already removed) inside a `BEGIN IONS` block to `spec`.
/// Returns `false` at `END IONS`. Shared by the streaming and the parallel reader so both parse
/// identically.
#[inline]
fn apply_body_line(spec: &mut Spectrum, line: &str) -> bool {
    let t = line.trim();
    if t == "END IONS" {
        return false;
    }
    if t.is_empty() {
        return true;
    }
    if let Some(v) = t.strip_prefix("TITLE=") {
        spec.title = Some(v.to_string());
    } else if let Some(v) = t.strip_prefix("SCANS=") {
        spec.scan = Some(v.to_string());
    } else if let Some(v) = t.strip_prefix("PEPMASS=") {
        spec.precursor_mz = v.split_whitespace().next().and_then(|x| x.parse().ok());
    } else if let Some(v) = t.strip_prefix("CHARGE=") {
        spec.charge = v.trim_end_matches(['+', '-']).parse().ok();
    } else if t.as_bytes()[0].is_ascii_digit() {
        let mut it = t.split_whitespace();
        if let Some(mz) = it.next().and_then(|x| x.parse::<f64>().ok()) {
            let intensity = it.next().and_then(|x| x.parse::<f64>().ok()).unwrap_or(0.0);
            spec.peaks.push(Peak { mz, intensity });
        }
    }
    // any other header (RTINSECONDS=, etc.) is ignored
    true
}

impl<R: BufRead> Iterator for MgfReader<R> {
    type Item = io::Result<Spectrum>;
    fn next(&mut self) -> Option<Self::Item> {
        self.read_next().transpose()
    }
}

/// Read all spectra from an MGF file. The file is read whole; spectrum blocks are located by one
/// sequential pass with the streaming reader's state machine and then parsed in parallel (rayon),
/// with exactly the streaming reader's line rules, so the result is identical to
/// `MgfReader::new(BufReader::new(file)).collect()`. A file that is not valid UTF-8 falls back to
/// the streaming reader (which reports the error the same way it always did).
pub fn read_mgf_file<P: AsRef<Path>>(path: P) -> io::Result<Vec<Spectrum>> {
    let bytes = std::fs::read(path)?;
    match std::str::from_utf8(&bytes) {
        Ok(text) => Ok(parse_mgf_parallel(text)),
        Err(_) => MgfReader::new(io::Cursor::new(bytes)).collect(),
    }
}

/// The lines of `text` as `BufRead::read_line` yields them (split after each `\n`, last line
/// possibly unterminated), each with its `\r`/`\n` terminator trimmed, plus its start offset.
fn lines_with_offsets(text: &str) -> impl Iterator<Item = (usize, &str)> {
    let b = text.as_bytes();
    let mut pos = 0usize;
    std::iter::from_fn(move || {
        if pos >= b.len() {
            return None;
        }
        let start = pos;
        let end = match b[pos..].iter().position(|&c| c == b'\n') {
            Some(i) => pos + i + 1,
            None => b.len(),
        };
        pos = end;
        Some((start, text[start..end].trim_end_matches(['\r', '\n'])))
    })
}

/// Whether a (terminator-trimmed) line could be a block marker; a line starting with an ASCII
/// digit (a peak) never trims to `BEGIN IONS` / `END IONS`, which skips the trim for most lines.
#[inline]
fn is_marker(line: &str, marker: &str) -> bool {
    !line.as_bytes().first().is_some_and(|c| c.is_ascii_digit()) && line.trim() == marker
}

fn parse_mgf_parallel(text: &str) -> Vec<Spectrum> {
    use rayon::prelude::*;
    // Body byte ranges, following the streaming reader: outside a block skip to a line that trims
    // to BEGIN IONS; inside, stop at a line that trims to END IONS or at EOF.
    let mut blocks: Vec<(usize, usize)> = Vec::new();
    let mut open: Option<usize> = None;
    for (start, line) in lines_with_offsets(text) {
        match open {
            None => {
                if is_marker(line, "BEGIN IONS") {
                    open = Some(
                        start
                            + text[start..]
                                .find('\n')
                                .map_or(text.len() - start, |i| i + 1),
                    );
                }
            }
            Some(b0) => {
                if is_marker(line, "END IONS") {
                    blocks.push((b0, start));
                    open = None;
                }
            }
        }
    }
    if let Some(b0) = open {
        blocks.push((b0, text.len()));
    }
    blocks
        .par_iter()
        .enumerate()
        .map(|(index, &(b0, b1))| {
            let mut spec = Spectrum {
                index,
                ..Default::default()
            };
            for (_, line) in lines_with_offsets(&text[b0.min(b1)..b1]) {
                apply_body_line(&mut spec, line);
            }
            spec
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const SAMPLE: &str = "\
BEGIN IONS
TITLE=scan=144
PEPMASS=319.7108
SCANS=144
CHARGE=2+
105.01811 73.3644
126.12766 554.4465
END IONS
BEGIN IONS
TITLE=scan=145
PEPMASS=400.5 1000.0
CHARGE=3+
200.0 10.0
END IONS
";

    #[test]
    fn parses_sample() {
        let specs: Vec<_> = MgfReader::new(Cursor::new(SAMPLE))
            .map(|s| s.unwrap())
            .collect();
        assert_eq!(specs.len(), 2);
        assert_eq!(specs[0].index, 0);
        assert_eq!(specs[0].scan.as_deref(), Some("144"));
        assert_eq!(specs[0].charge, Some(2));
        assert_eq!(specs[0].precursor_mz, Some(319.7108));
        assert_eq!(specs[0].n_peaks(), 2);
        assert_eq!(specs[1].charge, Some(3));
        assert_eq!(specs[1].precursor_mz, Some(400.5));
    }

    #[test]
    fn canonical_string_format() {
        let s = Spectrum {
            peaks: vec![
                Peak {
                    mz: 105.01811,
                    intensity: 73.3644,
                },
                Peak {
                    mz: 126.12766,
                    intensity: 554.4465,
                },
            ],
            ..Default::default()
        };
        assert_eq!(
            s.canonical_peak_string(),
            "105.01811 73.36440\n126.12766 554.44650"
        );
    }

    /// The parallel whole-file reader must equal the streaming reader on awkward inputs.
    #[test]
    fn parallel_reader_matches_streaming_reader() {
        let cases: &[&str] = &[
            SAMPLE,
            "",
            "no blocks here\n",
            "junk\r\nBEGIN IONS\r\nTITLE=a\r\nPEPMASS=1.5\r\n10 1\r\nEND IONS\r\njunk\r\n",
            // missing END IONS: the next BEGIN IONS is just an ignored line inside the block
            "BEGIN IONS\nSCANS=1\n100 1\nBEGIN IONS\nSCANS=2\n200 2\nEND IONS\n",
            // EOF inside a block, unterminated last line
            "BEGIN IONS\nSCANS=3\n300 3",
            "  BEGIN IONS  \n\n  TITLE=x \n 5 6 \n\t END IONS\nBEGIN IONS\nEND IONS\n",
            "END IONS\nBEGIN IONS\nCHARGE=2-\n1e2 3\n7\nEND IONS",
            "BEGIN IONS\nEND IONS\nBEGIN IONS\n",
        ];
        for (i, text) in cases.iter().enumerate() {
            let streamed: Vec<Spectrum> = MgfReader::new(Cursor::new(text.as_bytes()))
                .collect::<io::Result<_>>()
                .unwrap();
            assert_eq!(parse_mgf_parallel(text), streamed, "case {i}");
        }
    }
}
