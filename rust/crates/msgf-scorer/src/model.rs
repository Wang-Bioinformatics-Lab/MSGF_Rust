//! Decoding of the `.param` scoring model (big-endian positional stream; layout documented in the
//! model package's `param-format` interface document) and the derived per-partition tables.

use std::fmt;
use std::path::Path;

/// Why a model could not be decoded.
#[derive(Debug)]
pub enum ModelError {
    Io(std::io::Error),
    /// The stream ended before the named field.
    Truncated(&'static str),
    /// A count or value that cannot be valid.
    Invalid(String),
    /// The trailing sentinel was not `0x7FFFFFFF` (the parse desynchronised) or bytes remain.
    Sentinel(String),
    /// A model without error tables (`E = 0`): edge scoring is undefined, not supported.
    NoErrorTables,
}

impl fmt::Display for ModelError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ModelError::Io(e) => write!(f, "model read failed: {e}"),
            ModelError::Truncated(what) => write!(f, "model truncated while reading {what}"),
            ModelError::Invalid(s) => write!(f, "invalid model: {s}"),
            ModelError::Sentinel(s) => write!(f, "model sentinel check failed: {s}"),
            ModelError::NoErrorTables => {
                write!(f, "model has no error tables (E = 0); not supported")
            }
        }
    }
}

impl std::error::Error for ModelError {}

/// One fragment ion type of a partition.
#[derive(Clone, Debug)]
pub struct IonType {
    pub prefix: bool,
    pub charge: i32,
    pub offset: f32,
    pub frequency: f32,
}

impl IonType {
    /// Join key standing for the derived name `{P|S}_{charge}_{floor(offset + 0.5)}`.
    fn name_key(&self) -> (bool, i32, i64) {
        (
            self.prefix,
            self.charge,
            (self.offset + 0.5f32).floor() as i64,
        )
    }
}

/// One (charge, parent-mass boundary, segment) partition with its derived tables.
#[derive(Clone, Debug)]
pub struct Partition {
    pub charge: i32,
    pub boundary: f32,
    pub segment: i32,
    pub ions: Vec<IonType>,
    /// Per ion type, `Rmax + 1` rank log-likelihoods (already the row of the first type sharing
    /// this type's name).
    pub(crate) rank_ll: Vec<Vec<f32>>,
    /// `ln(signal/noise)` per error bin (`2E + 1` values), f32.
    pub(crate) error_ll: Vec<f32>,
    /// The four ion-existence values (zeros floored to 0.001).
    pub(crate) existence: [f32; 4],
}

/// A decoded scoring model.
#[derive(Clone, Debug)]
pub struct ScoreModel {
    pub version: i32,
    pub activation: String,
    pub instrument: String,
    pub enzyme: Option<String>,
    pub protocol: Option<String>,
    pub tol_is_ppm: bool,
    pub tol_value: f32,
    pub reduce_isotopes: bool,
    pub reduce_tol: f32,
    pub segments: i32,
    pub max_rank: i32,
    pub error_scale: i32,
    /// Sorted by `(charge, segment, boundary)`, unique.
    pub partitions: Vec<Partition>,
    /// `(charge, reduced charge, offset)` of each precursor-offset entry, in file order.
    pub precursor_offsets: Vec<(i32, i32, f32)>,
}

struct Reader<'a> {
    b: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    fn take(&mut self, n: usize, what: &'static str) -> Result<&'a [u8], ModelError> {
        if self.pos + n > self.b.len() {
            return Err(ModelError::Truncated(what));
        }
        let s = &self.b[self.pos..self.pos + n];
        self.pos += n;
        Ok(s)
    }
    fn i32(&mut self, what: &'static str) -> Result<i32, ModelError> {
        Ok(i32::from_be_bytes(self.take(4, what)?.try_into().unwrap()))
    }
    fn f32(&mut self, what: &'static str) -> Result<f32, ModelError> {
        Ok(f32::from_be_bytes(self.take(4, what)?.try_into().unwrap()))
    }
    fn bool(&mut self, what: &'static str) -> Result<bool, ModelError> {
        Ok(self.take(1, what)?[0] != 0)
    }
    fn count(&mut self, what: &'static str) -> Result<usize, ModelError> {
        let n = self.i32(what)?;
        if n < 0 || n as usize > self.b.len() {
            return Err(ModelError::Invalid(format!("{what} count {n}")));
        }
        Ok(n as usize)
    }
    /// Length-prefixed UTF-16BE string; length 0 = absent.
    fn string(&mut self, what: &'static str) -> Result<Option<String>, ModelError> {
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

impl ScoreModel {
    /// Read and decode a `.param` file.
    pub fn from_file(path: impl AsRef<Path>) -> Result<ScoreModel, ModelError> {
        let bytes = std::fs::read(path).map_err(ModelError::Io)?;
        ScoreModel::from_bytes(&bytes)
    }

    /// Decode a `.param` byte stream. Rejects a bad sentinel, trailing bytes and `E = 0`.
    pub fn from_bytes(bytes: &[u8]) -> Result<ScoreModel, ModelError> {
        let mut r = Reader { b: bytes, pos: 0 };
        let version = r.i32("version")?;
        let activation = r.string("activation")?.unwrap_or_default();
        let instrument = r.string("instrument")?.unwrap_or_default();
        let enzyme = r.string("enzyme")?;
        let protocol = r.string("protocol")?;
        let tol_is_ppm = r.bool("tolerance unit")?;
        let tol_value = r.f32("tolerance")?;
        let reduce_isotopes = r.bool("reduce flag")?;
        let reduce_tol = r.f32("reduce tolerance")?;

        let n_hist = r.count("charge histogram")?;
        for _ in 0..n_hist {
            r.i32("charge histogram")?;
            r.i32("charge histogram")?;
        }

        let n_part = r.count("partitions")?;
        let segments = r.i32("segment count")?;
        if segments < 1 {
            return Err(ModelError::Invalid(format!("segment count {segments}")));
        }
        let mut keys = Vec::with_capacity(n_part);
        for _ in 0..n_part {
            let z = r.i32("partition charge")?;
            let b = r.f32("partition boundary")?;
            let s = r.i32("partition segment")?;
            keys.push((z, s, b));
        }
        keys.sort_by(|a, b| {
            (a.0, a.1)
                .cmp(&(b.0, b.1))
                .then(a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal))
        });
        keys.dedup_by(|a, b| a.0 == b.0 && a.1 == b.1 && a.2.to_bits() == b.2.to_bits());

        let n_off = r.count("precursor offsets")?;
        let mut precursor_offsets = Vec::with_capacity(n_off);
        for _ in 0..n_off {
            let z = r.i32("precursor offset")?;
            let rz = r.i32("precursor offset")?;
            let off = r.f32("precursor offset")?;
            r.bool("precursor offset")?;
            r.f32("precursor offset")?;
            r.f32("precursor offset")?;
            precursor_offsets.push((z, rz, off));
        }

        let mut ion_lists = Vec::with_capacity(keys.len());
        for _ in 0..keys.len() {
            let n = r.count("ion types")?;
            let mut v = Vec::with_capacity(n);
            for _ in 0..n {
                let prefix = r.bool("ion type")?;
                let charge = r.i32("ion type")?;
                let offset = r.f32("ion type")?;
                let frequency = r.f32("ion type")?;
                if charge < 1 {
                    return Err(ModelError::Invalid(format!("ion charge {charge}")));
                }
                v.push(IonType {
                    prefix,
                    charge,
                    offset,
                    frequency,
                });
            }
            ion_lists.push(v);
        }

        let max_rank = r.i32("max rank")?;
        if max_rank < 1 {
            return Err(ModelError::Invalid(format!("max rank {max_rank}")));
        }
        let cols = max_rank as usize + 1;
        // Raw rank rows per partition: (ion rows, noise row).
        let mut rank_rows: Vec<Option<(Vec<Vec<f32>>, Vec<f32>)>> = Vec::with_capacity(keys.len());
        for ions in &ion_lists {
            if ions.is_empty() {
                rank_rows.push(None);
                continue;
            }
            let mut rows = Vec::with_capacity(ions.len());
            for _ in ions {
                rows.push(
                    (0..cols)
                        .map(|_| r.f32("rank row"))
                        .collect::<Result<Vec<_>, _>>()?,
                );
            }
            let noise = (0..cols)
                .map(|_| r.f32("noise row"))
                .collect::<Result<Vec<_>, _>>()?;
            rank_rows.push(Some((rows, noise)));
        }

        let error_scale = r.i32("error scale")?;
        if error_scale < 0 {
            return Err(ModelError::Invalid(format!("error scale {error_scale}")));
        }
        let mut err_tables = Vec::with_capacity(keys.len());
        if error_scale > 0 {
            let w = 2 * error_scale as usize + 1;
            for _ in 0..keys.len() {
                let signal = (0..w)
                    .map(|_| r.f32("error signal"))
                    .collect::<Result<Vec<_>, _>>()?;
                let noise = (0..w)
                    .map(|_| r.f32("error noise"))
                    .collect::<Result<Vec<_>, _>>()?;
                let mut ex = [0f32; 4];
                for e in ex.iter_mut() {
                    let v = r.f32("existence")?;
                    *e = if v == 0.0 { 0.001 } else { v };
                }
                err_tables.push((signal, noise, ex));
            }
        }
        let sentinel = r.i32("sentinel")?;
        if sentinel != 0x7FFF_FFFF {
            return Err(ModelError::Sentinel(format!("found {sentinel:#x}")));
        }
        if r.pos != bytes.len() {
            return Err(ModelError::Sentinel(format!(
                "{} trailing bytes",
                bytes.len() - r.pos
            )));
        }
        if error_scale == 0 {
            return Err(ModelError::NoErrorTables);
        }

        let mut partitions = Vec::with_capacity(keys.len());
        for (i, ((z, s, b), ions)) in keys.into_iter().zip(ion_lists).enumerate() {
            let mut rank_ll = Vec::with_capacity(ions.len());
            if let Some((rows, noise)) = &rank_rows[i] {
                for ion in &ions {
                    // Row of the first type with the same derived name.
                    let first = ions
                        .iter()
                        .position(|o| o.name_key() == ion.name_key())
                        .unwrap();
                    let row = &rows[first];
                    let div = ion.charge.min(segments) as f32;
                    rank_ll.push(
                        (0..cols)
                            .map(|b| ((row[b] / (noise[b] * div)) as f64).ln() as f32)
                            .collect(),
                    );
                }
            }
            let (signal, noise, existence) = &err_tables[i];
            let error_ll = signal
                .iter()
                .zip(noise)
                .map(|(s, n)| (*s as f64 / *n as f64).ln() as f32)
                .collect();
            partitions.push(Partition {
                charge: z,
                boundary: b,
                segment: s,
                ions,
                rank_ll,
                error_ll,
                existence: *existence,
            });
        }

        Ok(ScoreModel {
            version,
            activation,
            instrument,
            enzyme,
            protocol,
            tol_is_ppm,
            tol_value,
            reduce_isotopes,
            reduce_tol,
            segments,
            max_rank,
            error_scale,
            partitions,
            precursor_offsets,
        })
    }

    /// Fragment tolerance (Da) at theoretical m/z `x`.
    #[inline]
    pub(crate) fn tolerance_at(&self, x: f32) -> f32 {
        if self.tol_is_ppm {
            (x as f64 * self.tol_value as f64 * 1e-6) as f32
        } else {
            self.tol_value
        }
    }

    /// The partition floor: index of the last partition whose `(charge, segment, boundary)` is
    /// `<=` `(z, s, mass)`, lexicographically.
    fn floor_index(&self, z: i32, s: i32, mass: f32) -> Option<usize> {
        let n = self.partitions.partition_point(|p| {
            (p.charge, p.segment) < (z, s)
                || ((p.charge, p.segment) == (z, s) && p.boundary <= mass)
        });
        n.checked_sub(1)
    }

    /// The partition used for segment `s` of a spectrum with charge `z` and parent mass `mass`.
    pub(crate) fn partition_for(&self, z: i32, s: i32, mass: f32) -> Option<usize> {
        match self.floor_index(z, s, mass) {
            None => {
                let z0 = self.partitions.first()?.charge;
                self.floor_index(z0, s, mass)
            }
            Some(i) if self.partitions[i].charge == z => Some(i),
            Some(i) => self.floor_index(self.partitions[i].charge, s, mass),
        }
    }
}
