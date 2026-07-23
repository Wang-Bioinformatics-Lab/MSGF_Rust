//! msgf-scorer — reader for MS-GF+ binary scoring models (`.param`).
//!
//! The `.param` files are written by `NewRankScorer.writeParameters` (Java `DataOutputStream`,
//! big-endian) and hold a trained rank-scoring model per (activation, instrument, enzyme,
//! protocol). This module decodes that format faithfully. The stream ends with an
//! `0x7FFFFFFF` sentinel that MS-GF+ itself checks; we validate it too, so a misaligned parse
//! is caught rather than silently accepted.
//!
//! Validated against `validation/golden/models/*.model.golden.json` (derived from the
//! authoritative `writeParametersPlainText` dump) in `tests/golden_model.rs`.
//!
//! Per-node spectrum scoring (RawScore) is built on top of this model next.

use msgf_chem::Tolerance;
use std::fs;
use std::io;
use std::path::Path;

/// Terminator written after all model data (`Integer.MAX_VALUE`).
const TERMINATOR: i32 = i32::MAX;

/// A scoring partition: (precursor charge, parent mass boundary, mass segment index).
#[derive(Debug, Clone, PartialEq)]
pub struct Partition {
    pub charge: i32,
    pub parent_mass: f32,
    pub seg: i32,
}

/// A fragment-ion offset-frequency entry (one theoretical ion type in a partition).
#[derive(Debug, Clone, PartialEq)]
pub struct FragOff {
    pub is_prefix: bool,
    pub charge: i32,
    pub offset: f32,
    pub frequency: f32,
    /// Name as MS-GF+ builds it: `"P_{charge}_{round(offset)}"` / `"S_{charge}_{round(offset)}"`.
    pub name: String,
}

/// A precursor offset-frequency entry.
#[derive(Debug, Clone, PartialEq)]
pub struct PrecursorOff {
    pub charge: i32,
    pub reduced_charge: i32,
    pub offset: f32,
    pub tol_ppm: bool,
    pub tol_val: f32,
    pub frequency: f32,
}

/// Rank-distribution table for one partition: each ion type (fragment ions in read order, then
/// `noise`) maps to `max_rank + 1` frequencies indexed by observed peak rank.
#[derive(Debug, Clone, PartialEq)]
pub struct RankDist {
    pub partition_index: usize,
    pub ions: Vec<(String, Vec<f32>)>,
}

/// Mass-error distribution for one partition.
#[derive(Debug, Clone, PartialEq)]
pub struct ErrorDist {
    pub signal: Vec<f32>,
    pub noise: Vec<f32>,
    pub ion_existence: [f32; 4],
}

/// A fully decoded MS-GF+ scoring model.
#[derive(Debug, Clone)]
pub struct ScoringModel {
    pub version: i32,
    pub activation: String,
    pub instrument: String,
    pub enzyme: Option<String>,
    pub protocol: Option<String>, // None == Automatic
    pub mme: Tolerance,
    pub apply_deconvolution: bool,
    pub deconvolution_error_tolerance: f32,
    pub charge_histogram: Vec<(i32, i32)>,
    pub num_segments: i32,
    /// Partitions in MS-GF+ `TreeSet` order: (charge, seg, parent_mass).
    pub partitions: Vec<Partition>,
    pub precursor_off: Vec<PrecursorOff>,
    /// Fragment offset frequencies, parallel to `partitions`.
    pub frag_off: Vec<Vec<FragOff>>,
    pub max_rank: i32,
    /// Rank distributions (only for partitions that have ≥1 fragment ion type).
    pub rank_dist: Vec<RankDist>,
    pub error_scaling_factor: i32,
    /// Error distributions, parallel to `partitions` (empty if `error_scaling_factor == 0`).
    pub error_dist: Vec<ErrorDist>,
}

/// Failure decoding a `.param` stream.
#[derive(Debug)]
pub enum ParamError {
    /// Ran off the end of the buffer at `pos`.
    UnexpectedEof {
        pos: usize,
        need: usize,
    },
    /// The trailing sentinel was wrong — the parse desynced somewhere.
    BadTerminator {
        got: i32,
        pos: usize,
    },
    Io(io::Error),
}

impl std::fmt::Display for ParamError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParamError::UnexpectedEof { pos, need } => {
                write!(f, "unexpected EOF at {pos} (need {need} bytes)")
            }
            ParamError::BadTerminator { got, pos } => {
                write!(f, "bad terminator {got:#x} at {pos} (parse desynced)")
            }
            ParamError::Io(e) => write!(f, "io: {e}"),
        }
    }
}
impl std::error::Error for ParamError {}
impl From<io::Error> for ParamError {
    fn from(e: io::Error) -> Self {
        ParamError::Io(e)
    }
}

/// Big-endian reader mirroring Java `DataInputStream`.
struct Reader<'a> {
    b: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    fn new(b: &'a [u8]) -> Self {
        Self { b, pos: 0 }
    }
    fn take(&mut self, n: usize) -> Result<&'a [u8], ParamError> {
        if self.pos + n > self.b.len() {
            return Err(ParamError::UnexpectedEof {
                pos: self.pos,
                need: n,
            });
        }
        let s = &self.b[self.pos..self.pos + n];
        self.pos += n;
        Ok(s)
    }
    fn u8(&mut self) -> Result<u8, ParamError> {
        Ok(self.take(1)?[0])
    }
    fn bool(&mut self) -> Result<bool, ParamError> {
        Ok(self.u8()? != 0)
    }
    fn i32(&mut self) -> Result<i32, ParamError> {
        Ok(i32::from_be_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn f32(&mut self) -> Result<f32, ParamError> {
        Ok(f32::from_be_bytes(self.take(4)?.try_into().unwrap()))
    }
    /// Java `writeByte(len)` + `writeChars` (len UTF-16BE chars).
    fn jstring(&mut self) -> Result<String, ParamError> {
        let len = self.u8()? as usize;
        let mut s = String::with_capacity(len);
        for _ in 0..len {
            let hi = self.u8()? as u16;
            let lo = self.u8()? as u16;
            s.push(char::from_u32(((hi << 8) | lo) as u32).unwrap_or('\u{FFFD}'));
        }
        Ok(s)
    }
}

/// Java `Math.round` for the ion-name suffix: `floor(x + 0.5)`.
fn java_round(x: f32) -> i64 {
    (x + 0.5).floor() as i64
}

/// Decode a `.param` model from raw bytes.
pub fn read_param(bytes: &[u8]) -> Result<ScoringModel, ParamError> {
    let mut r = Reader::new(bytes);

    let version = r.i32()?;
    let activation = r.jstring()?;
    let instrument = r.jstring()?;
    let enzyme = {
        // length-0 => absent; jstring already handles len byte
        let start = r.pos;
        let s = r.jstring()?;
        if r.pos == start + 1 {
            None
        } else {
            Some(s)
        } // len byte was 0
    };
    let protocol = {
        let start = r.pos;
        let s = r.jstring()?;
        if r.pos == start + 1 {
            None
        } else {
            Some(s)
        }
    };

    let mme_ppm = r.bool()?;
    let mme_val = r.f32()?;
    let mme = if mme_ppm {
        Tolerance::ppm(mme_val as f64)
    } else {
        Tolerance::da(mme_val as f64)
    };

    let apply_deconvolution = r.bool()?;
    let deconvolution_error_tolerance = r.f32()?;

    // charge histogram
    let n = r.i32()? as usize;
    let mut charge_histogram = Vec::with_capacity(n);
    for _ in 0..n {
        let charge = r.i32()?;
        let count = r.i32()?;
        charge_histogram.push((charge, count));
    }

    // partitions
    let n = r.i32()? as usize;
    let num_segments = r.i32()?;
    let mut partitions = Vec::with_capacity(n);
    for _ in 0..n {
        let charge = r.i32()?;
        let parent_mass = r.f32()?;
        let seg = r.i32()?;
        partitions.push(Partition {
            charge,
            parent_mass,
            seg,
        });
    }
    // MS-GF+ stores these in a TreeSet: order by (charge, seg, parent_mass), unique.
    partitions.sort_by(|a, b| {
        a.charge
            .cmp(&b.charge)
            .then(a.seg.cmp(&b.seg))
            .then(a.parent_mass.partial_cmp(&b.parent_mass).unwrap())
    });
    partitions
        .dedup_by(|a, b| a.charge == b.charge && a.seg == b.seg && a.parent_mass == b.parent_mass);

    // precursor offset frequencies
    let n = r.i32()? as usize;
    let mut precursor_off = Vec::with_capacity(n);
    for _ in 0..n {
        let charge = r.i32()?;
        let reduced_charge = r.i32()?;
        let offset = r.f32()?;
        let tol_ppm = r.bool()?;
        let tol_val = r.f32()?;
        let frequency = r.f32()?;
        precursor_off.push(PrecursorOff {
            charge,
            reduced_charge,
            offset,
            tol_ppm,
            tol_val,
            frequency,
        });
    }

    // fragment offset frequencies — one block per partition, in sorted order
    let mut frag_off: Vec<Vec<FragOff>> = Vec::with_capacity(partitions.len());
    for _ in 0..partitions.len() {
        let size = r.i32()? as usize;
        let mut block = Vec::with_capacity(size);
        for _ in 0..size {
            let is_prefix = r.bool()?;
            let charge = r.i32()?;
            let offset = r.f32()?;
            let frequency = r.f32()?;
            let tag = if is_prefix { 'P' } else { 'S' };
            let name = format!("{tag}_{charge}_{}", java_round(offset));
            block.push(FragOff {
                is_prefix,
                charge,
                offset,
                frequency,
                name,
            });
        }
        frag_off.push(block);
    }

    // rank distributions — per partition with ≥1 ion type, ions then NOISE
    let max_rank = r.i32()?;
    let ncols = (max_rank + 1) as usize;
    let mut rank_dist = Vec::new();
    for (pi, block) in frag_off.iter().enumerate() {
        if block.is_empty() {
            continue; // getIonTypes empty => MS-GF+ skips this partition
        }
        let mut ions = Vec::with_capacity(block.len() + 1);
        for fo in block {
            let mut freqs = Vec::with_capacity(ncols);
            for _ in 0..ncols {
                freqs.push(r.f32()?);
            }
            ions.push((fo.name.clone(), freqs));
        }
        let mut noise = Vec::with_capacity(ncols);
        for _ in 0..ncols {
            noise.push(r.f32()?);
        }
        ions.push(("noise".to_string(), noise));
        rank_dist.push(RankDist {
            partition_index: pi,
            ions,
        });
    }

    // error distributions
    let error_scaling_factor = r.i32()?;
    let mut error_dist = Vec::new();
    if error_scaling_factor > 0 {
        let width = (error_scaling_factor * 2 + 1) as usize;
        for _ in 0..partitions.len() {
            let mut signal = Vec::with_capacity(width);
            for _ in 0..width {
                signal.push(r.f32()?);
            }
            let mut noise = Vec::with_capacity(width);
            for _ in 0..width {
                noise.push(r.f32()?);
            }
            let mut ion_existence = [0.0f32; 4];
            for slot in &mut ion_existence {
                let v = r.f32()?;
                *slot = if v == 0.0 { 0.001 } else { v };
            }
            error_dist.push(ErrorDist {
                signal,
                noise,
                ion_existence,
            });
        }
    }

    // sentinel — proves the whole parse stayed aligned
    let term = r.i32()?;
    if term != TERMINATOR {
        return Err(ParamError::BadTerminator {
            got: term,
            pos: r.pos - 4,
        });
    }

    Ok(ScoringModel {
        version,
        activation,
        instrument,
        enzyme,
        protocol,
        mme,
        apply_deconvolution,
        deconvolution_error_tolerance,
        charge_histogram,
        num_segments,
        partitions,
        precursor_off,
        frag_off,
        max_rank,
        rank_dist,
        error_scaling_factor,
        error_dist,
    })
}

/// Read and decode a `.param` file from disk.
pub fn read_param_file<P: AsRef<Path>>(path: P) -> Result<ScoringModel, ParamError> {
    read_param(&fs::read(path)?)
}

impl RankDist {
    /// Frequency row for an ion by name (`None` if the ion is not scored in this partition).
    fn row(&self, name: &str) -> Option<&[f32]> {
        self.ions
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, v)| v.as_slice())
    }
    /// The noise row (always stored last).
    fn noise_row(&self) -> &[f32] {
        &self
            .ions
            .last()
            .expect("rank distribution has a noise row")
            .1
    }
}

impl ScoringModel {
    /// Protocol name, mapping the absent case to MS-GF+'s `"Automatic"` default.
    pub fn protocol_name(&self) -> &str {
        self.protocol.as_deref().unwrap_or("Automatic")
    }

    /// Rank distribution for a partition index (only partitions with ≥1 ion have one).
    pub fn rank_dist_for(&self, partition_index: usize) -> Option<&RankDist> {
        self.rank_dist
            .iter()
            .find(|r| r.partition_index == partition_index)
    }

    /// Log-likelihood score for observing `ion`'s peak at 1-based `rank` (rank 1 = most intense)
    /// in partition `partition_index`. Mirrors `NewRankScorer.getNodeScore`.
    pub fn node_score(&self, partition_index: usize, ion: &FragOff, rank: i32) -> f32 {
        let idx = if rank > self.max_rank {
            (self.max_rank - 1) as usize
        } else {
            (rank - 1) as usize
        };
        self.score_from_table(partition_index, ion, idx)
    }

    /// Log-likelihood score for `ion`'s peak being absent (the `maxRank` bin).
    /// Mirrors `NewRankScorer.getMissingIonScore`.
    pub fn missing_ion_score(&self, partition_index: usize, ion: &FragOff) -> f32 {
        self.score_from_table(partition_index, ion, self.max_rank as usize)
    }

    /// `log( ionFreq[idx] / (noiseFreq[idx] * min(ionCharge, numSegments)) )`, computed in the
    /// same float order as Java (`NewRankScorer.getScoreFromTable`, `isError = false`).
    fn score_from_table(&self, partition_index: usize, ion: &FragOff, idx: usize) -> f32 {
        let rd = self
            .rank_dist_for(partition_index)
            .expect("partition has a rank distribution");
        let ion_freq = rd.row(&ion.name).expect("ion is scored in this partition")[idx];
        let noise = rd.noise_row()[idx] * ion.charge.min(self.num_segments) as f32;
        ((ion_freq / noise) as f64).ln() as f32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn java_round_matches() {
        assert_eq!(java_round(19.01839), 19);
        assert_eq!(java_round(1.9918417), 2);
        assert_eq!(java_round(-26.98709), -27);
        assert_eq!(java_round(1.007825), 1);
    }
}
