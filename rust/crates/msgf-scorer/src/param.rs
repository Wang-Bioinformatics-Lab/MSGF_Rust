//! The `.param` scoring-model file as data, and its decoder.
//!
//! Written from the format description in `docs/param-format.md` (big-endian positional stream,
//! sections §1–§8, the two read-side transforms). [`crate::write_param`] is the inverse. For
//! scoring, decode the same bytes with [`crate::ScoreModel::from_bytes`], which derives the
//! log-likelihood tables.

use msgf_chem::{Tolerance, Unit};
use std::fmt;
use std::path::Path;

/// The trailing sentinel of every `.param` stream (§8).
pub const TERMINATOR: i32 = 0x7FFF_FFFF;

/// One `(charge, parent-mass boundary, segment)` partition (§3).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Partition {
    pub charge: i32,
    pub parent_mass: f32,
    pub seg: i32,
}

/// One precursor-offset entry (§4).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrecursorOff {
    pub charge: i32,
    pub reduced_charge: i32,
    pub offset: f32,
    pub tol_ppm: bool,
    pub tol_val: f32,
    pub frequency: f32,
}

/// One fragment ion type of a partition (§5).
#[derive(Debug, Clone, PartialEq)]
pub struct FragOff {
    pub is_prefix: bool,
    pub charge: i32,
    pub offset: f32,
    pub frequency: f32,
    /// `{P|S}_{charge}_{floor(offset + 0.5)}` — derived on read, never stored.
    pub name: String,
}

impl FragOff {
    /// The derived ion name, `{P|S}_{charge}_{floor(offset + 0.5)}` (the addition in f32).
    pub fn derive_name(is_prefix: bool, charge: i32, offset: f32) -> String {
        format!(
            "{}_{}_{}",
            if is_prefix { 'P' } else { 'S' },
            charge,
            (offset + 0.5f32).floor() as i64
        )
    }
}

/// The rank rows of one partition with a non-empty ion block (§6): one `(name, row)` per ion
/// type in block order, then `("noise", row)`. Each row has `max_rank + 1` columns.
#[derive(Debug, Clone, PartialEq)]
pub struct RankDist {
    pub partition_index: usize,
    pub ions: Vec<(String, Vec<f32>)>,
}

/// The error / ion-existence tables of one partition (§7).
#[derive(Debug, Clone, PartialEq)]
pub struct ErrorDist {
    /// `2E + 1` values.
    pub signal: Vec<f32>,
    /// `2E + 1` values.
    pub noise: Vec<f32>,
    /// Zeros are floored to 0.001 on read.
    pub ion_existence: [f32; 4],
}

/// A `.param` model, field for field.
#[derive(Debug, Clone, PartialEq)]
pub struct ScoringModel {
    pub version: i32,
    pub activation: String,
    pub instrument: String,
    pub enzyme: Option<String>,
    /// `None` = the "Automatic" protocol.
    pub protocol: Option<String>,
    /// Fragment mass tolerance.
    pub mme: Tolerance,
    /// Isotope-cluster reduction on/off.
    pub apply_deconvolution: bool,
    /// Isotope-cluster reduction tolerance (Da).
    pub deconvolution_error_tolerance: f32,
    pub charge_histogram: Vec<(i32, i32)>,
    pub num_segments: i32,
    /// Sorted by `(charge, seg, parent_mass)`, unique.
    pub partitions: Vec<Partition>,
    pub precursor_off: Vec<PrecursorOff>,
    /// One block per partition, parallel to `partitions`.
    pub frag_off: Vec<Vec<FragOff>>,
    pub max_rank: i32,
    /// One entry per partition with a non-empty `frag_off` block, in partition order.
    pub rank_dist: Vec<RankDist>,
    /// `E`; 0 means §7 is absent.
    pub error_scaling_factor: i32,
    /// One per partition when `E > 0`, else empty.
    pub error_dist: Vec<ErrorDist>,
}

impl ScoringModel {
    /// The protocol name: the stored one, or `"Automatic"` when absent.
    pub fn protocol_name(&self) -> &str {
        self.protocol.as_deref().unwrap_or("Automatic")
    }

    /// The rank row of ion `ion` in partition `pi` and that partition's noise row, if scored.
    fn rows(&self, pi: usize, ion: &FragOff) -> Option<(&[f32], &[f32])> {
        let rd = self.rank_dist.iter().find(|r| r.partition_index == pi)?;
        let (_, noise) = rd.ions.last()?;
        let (_, row) = rd.ions.iter().find(|(n, _)| *n == ion.name)?;
        Some((row, noise))
    }

    /// `ln(ion[c] / (noise[c] * min(ion charge, num_segments)))` at column `c` (§6 note).
    fn table_score(&self, pi: usize, ion: &FragOff, col: usize) -> f32 {
        let Some((row, noise)) = self.rows(pi, ion) else {
            return 0.0;
        };
        let div = ion.charge.min(self.num_segments) as f32;
        ((row[col] / (noise[col] * div)) as f64).ln() as f32
    }

    /// Node score of ion `ion` of partition `pi` whose matched peak has intensity rank `rank`
    /// (1 = most intense; ranks beyond `max_rank - 1` share the last present column).
    pub fn node_score(&self, pi: usize, ion: &FragOff, rank: u32) -> f32 {
        let col = (rank.clamp(1, self.max_rank.max(1) as u32) - 1) as usize;
        self.table_score(pi, ion, col)
    }

    /// Score of ion `ion` of partition `pi` when no peak matches (the `max_rank` column).
    pub fn missing_ion_score(&self, pi: usize, ion: &FragOff) -> f32 {
        self.table_score(pi, ion, self.max_rank as usize)
    }
}

/// Why a `.param` stream could not be decoded.
#[derive(Debug)]
pub enum ParamError {
    Io(std::io::Error),
    /// The stream ended while reading the named field.
    Truncated(&'static str),
    /// A count or value that cannot be valid.
    Invalid(String),
    /// The trailing sentinel was wrong or bytes remain after it (the parse desynchronised).
    Sentinel(String),
}

impl fmt::Display for ParamError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParamError::Io(e) => write!(f, "model read failed: {e}"),
            ParamError::Truncated(what) => write!(f, "model truncated while reading {what}"),
            ParamError::Invalid(s) => write!(f, "invalid model: {s}"),
            ParamError::Sentinel(s) => write!(f, "model sentinel check failed: {s}"),
        }
    }
}

impl std::error::Error for ParamError {}

struct Reader<'a> {
    b: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    fn take(&mut self, n: usize, what: &'static str) -> Result<&'a [u8], ParamError> {
        if self.pos + n > self.b.len() {
            return Err(ParamError::Truncated(what));
        }
        let s = &self.b[self.pos..self.pos + n];
        self.pos += n;
        Ok(s)
    }
    fn i32(&mut self, what: &'static str) -> Result<i32, ParamError> {
        Ok(i32::from_be_bytes(self.take(4, what)?.try_into().unwrap()))
    }
    fn f32(&mut self, what: &'static str) -> Result<f32, ParamError> {
        Ok(f32::from_be_bytes(self.take(4, what)?.try_into().unwrap()))
    }
    fn bool(&mut self, what: &'static str) -> Result<bool, ParamError> {
        Ok(self.take(1, what)?[0] != 0)
    }
    fn count(&mut self, what: &'static str) -> Result<usize, ParamError> {
        let n = self.i32(what)?;
        if n < 0 || n as usize > self.b.len() {
            return Err(ParamError::Invalid(format!("{what} count {n}")));
        }
        Ok(n as usize)
    }
    fn f32s(&mut self, n: usize, what: &'static str) -> Result<Vec<f32>, ParamError> {
        (0..n).map(|_| self.f32(what)).collect()
    }
    /// Length-prefixed UTF-16BE string; length 0 = absent.
    fn string(&mut self, what: &'static str) -> Result<Option<String>, ParamError> {
        let n = self.take(1, what)?[0] as usize;
        if n == 0 {
            return Ok(None);
        }
        let raw = self.take(2 * n, what)?;
        let units: Vec<u16> = raw
            .chunks(2)
            .map(|c| u16::from_be_bytes([c[0], c[1]]))
            .collect();
        Ok(Some(String::from_utf16_lossy(&units)))
    }
}

/// Read and decode a `.param` file.
pub fn read_param_file(path: impl AsRef<Path>) -> Result<ScoringModel, ParamError> {
    read_param(&std::fs::read(path).map_err(ParamError::Io)?)
}

/// Decode a `.param` byte stream (`docs/param-format.md`).
pub fn read_param(bytes: &[u8]) -> Result<ScoringModel, ParamError> {
    let mut r = Reader { b: bytes, pos: 0 };

    // §1 header / identity
    let version = r.i32("version")?;
    let activation = r.string("activation")?.unwrap_or_default();
    let instrument = r.string("instrument")?.unwrap_or_default();
    let enzyme = r.string("enzyme")?;
    let protocol = r.string("protocol")?;
    let mme_ppm = r.bool("tolerance unit")?;
    let mme_value = r.f32("tolerance")?;
    let apply_deconvolution = r.bool("deconvolution flag")?;
    let deconvolution_error_tolerance = r.f32("deconvolution tolerance")?;
    let mme = Tolerance {
        value: mme_value as f64,
        unit: if mme_ppm { Unit::Ppm } else { Unit::Da },
    };

    // §2 charge histogram
    let n_hist = r.count("charge histogram")?;
    let mut charge_histogram = Vec::with_capacity(n_hist);
    for _ in 0..n_hist {
        let c = r.i32("charge histogram")?;
        let n = r.i32("charge histogram")?;
        charge_histogram.push((c, n));
    }

    // §3 partitions (sorted by (charge, seg, parent_mass), unique)
    let n_part = r.count("partitions")?;
    let num_segments = r.i32("segment count")?;
    let mut partitions = Vec::with_capacity(n_part);
    for _ in 0..n_part {
        let charge = r.i32("partition charge")?;
        let parent_mass = r.f32("partition boundary")?;
        let seg = r.i32("partition segment")?;
        partitions.push(Partition {
            charge,
            parent_mass,
            seg,
        });
    }
    partitions.sort_by(|a, b| {
        (a.charge, a.seg).cmp(&(b.charge, b.seg)).then(
            a.parent_mass
                .partial_cmp(&b.parent_mass)
                .unwrap_or(std::cmp::Ordering::Equal),
        )
    });
    partitions.dedup_by(|a, b| {
        a.charge == b.charge && a.seg == b.seg && a.parent_mass.to_bits() == b.parent_mass.to_bits()
    });

    // §4 precursor offsets
    let n_off = r.count("precursor offsets")?;
    let mut precursor_off = Vec::with_capacity(n_off);
    for _ in 0..n_off {
        precursor_off.push(PrecursorOff {
            charge: r.i32("precursor offset")?,
            reduced_charge: r.i32("precursor offset")?,
            offset: r.f32("precursor offset")?,
            tol_ppm: r.bool("precursor offset")?,
            tol_val: r.f32("precursor offset")?,
            frequency: r.f32("precursor offset")?,
        });
    }

    // §5 fragment offsets, one block per partition
    let mut frag_off = Vec::with_capacity(partitions.len());
    for _ in 0..partitions.len() {
        let n = r.count("ion types")?;
        let mut block = Vec::with_capacity(n);
        for _ in 0..n {
            let is_prefix = r.bool("ion type")?;
            let charge = r.i32("ion type")?;
            let offset = r.f32("ion type")?;
            let frequency = r.f32("ion type")?;
            block.push(FragOff {
                is_prefix,
                charge,
                offset,
                frequency,
                name: FragOff::derive_name(is_prefix, charge, offset),
            });
        }
        frag_off.push(block);
    }

    // §6 rank distributions, partitions with a non-empty block only; noise row last
    let max_rank = r.i32("max rank")?;
    if max_rank < 0 {
        return Err(ParamError::Invalid(format!("max rank {max_rank}")));
    }
    let cols = max_rank as usize + 1;
    let mut rank_dist = Vec::new();
    for (pi, block) in frag_off.iter().enumerate() {
        if block.is_empty() {
            continue;
        }
        let mut ions = Vec::with_capacity(block.len() + 1);
        for fo in block {
            ions.push((fo.name.clone(), r.f32s(cols, "rank row")?));
        }
        ions.push(("noise".to_string(), r.f32s(cols, "noise row")?));
        rank_dist.push(RankDist {
            partition_index: pi,
            ions,
        });
    }

    // §7 error / ion-existence tables, every partition, present iff E > 0
    let error_scaling_factor = r.i32("error scale")?;
    if error_scaling_factor < 0 {
        return Err(ParamError::Invalid(format!(
            "error scale {error_scaling_factor}"
        )));
    }
    let mut error_dist = Vec::new();
    if error_scaling_factor > 0 {
        let w = 2 * error_scaling_factor as usize + 1;
        for _ in 0..partitions.len() {
            let signal = r.f32s(w, "error signal")?;
            let noise = r.f32s(w, "error noise")?;
            let mut ion_existence = [0f32; 4];
            for e in ion_existence.iter_mut() {
                let v = r.f32("ion existence")?;
                *e = if v == 0.0 { 0.001 } else { v };
            }
            error_dist.push(ErrorDist {
                signal,
                noise,
                ion_existence,
            });
        }
    }

    // §8 sentinel
    let sentinel = r.i32("sentinel")?;
    if sentinel != TERMINATOR {
        return Err(ParamError::Sentinel(format!("found {sentinel:#x}")));
    }
    if r.pos != bytes.len() {
        return Err(ParamError::Sentinel(format!(
            "{} trailing bytes",
            bytes.len() - r.pos
        )));
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
