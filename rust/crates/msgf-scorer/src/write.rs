//! Encoder for the `.param` scoring-model format — the inverse of [`crate::read_param`].
//!
//! It serialises any in-memory [`ScoringModel`] into the on-disk byte format documented in
//! `docs/param-format.md`, which is how `msgf-train` emits the models this project ships. The
//! `.param` *file format* is an interface; a model we count ourselves carries no upstream licence.
//!
//! # Fidelity
//!
//! `read_param(write_param(m)) == m` for every model the reader produces. The encoding is
//! big-endian scalars, a length byte + UTF-16BE code units for strings, per-partition parallel
//! sections in the sorted partition order, and the trailing `0x7FFFFFFF` sentinel.

use crate::ScoringModel;
use msgf_chem::Unit;
use std::io;
use std::path::Path;

/// Big-endian sink (the write side of [`crate::read_param`]).
struct Writer {
    b: Vec<u8>,
}

impl Writer {
    fn new() -> Self {
        Self { b: Vec::new() }
    }
    fn u8(&mut self, v: u8) {
        self.b.push(v);
    }
    fn bool(&mut self, v: bool) {
        self.b.push(v as u8);
    }
    fn i32(&mut self, v: i32) {
        self.b.extend_from_slice(&v.to_be_bytes());
    }
    fn f32(&mut self, v: f32) {
        self.b.extend_from_slice(&v.to_be_bytes());
    }
    /// Length byte + UTF-16BE code units (the format's string encoding).
    fn jstring(&mut self, s: &str) {
        let units: Vec<u16> = s.encode_utf16().collect();
        self.u8(units.len() as u8);
        for u in units {
            self.b.extend_from_slice(&u.to_be_bytes());
        }
    }
    /// Optional string: `None` is a single `0` length byte (how the reader detects absence).
    fn jstring_opt(&mut self, s: &Option<String>) {
        match s {
            Some(x) => self.jstring(x),
            None => self.u8(0),
        }
    }
}

/// Serialise a [`ScoringModel`] into `.param` bytes. Exact inverse of [`crate::read_param`].
///
/// `crate::read_param(&write_param(m))` reproduces `m` for any model obtained from the reader.
/// Note the reader re-derives [`crate::FragOff::name`] from `(is_prefix, charge, offset)` and floors
/// zero ion-existence entries to `0.001`, so these fields are recomputed on the next read rather
/// than carried in the bytes — the round-trip is exact at the [`ScoringModel`] level, and
/// byte-for-byte for any file whose ion-existence entries are all non-zero.
pub fn write_param(m: &ScoringModel) -> Vec<u8> {
    let mut w = Writer::new();

    // header / identity
    w.i32(m.version);
    w.jstring(&m.activation);
    w.jstring(&m.instrument);
    w.jstring_opt(&m.enzyme);
    w.jstring_opt(&m.protocol);
    w.bool(m.mme.unit == Unit::Ppm);
    w.f32(m.mme.value as f32);
    w.bool(m.apply_deconvolution);
    w.f32(m.deconvolution_error_tolerance);

    // charge histogram
    w.i32(m.charge_histogram.len() as i32);
    for &(charge, count) in &m.charge_histogram {
        w.i32(charge);
        w.i32(count);
    }

    // partitions (already in the reader's sorted order: charge, seg, parent_mass)
    w.i32(m.partitions.len() as i32);
    w.i32(m.num_segments);
    for p in &m.partitions {
        w.i32(p.charge);
        w.f32(p.parent_mass);
        w.i32(p.seg);
    }

    // precursor offset frequencies
    w.i32(m.precursor_off.len() as i32);
    for po in &m.precursor_off {
        w.i32(po.charge);
        w.i32(po.reduced_charge);
        w.f32(po.offset);
        w.bool(po.tol_ppm);
        w.f32(po.tol_val);
        w.f32(po.frequency);
    }

    // fragment offset frequencies — one block per partition, in partition order
    for block in &m.frag_off {
        w.i32(block.len() as i32);
        for fo in block {
            w.bool(fo.is_prefix);
            w.i32(fo.charge);
            w.f32(fo.offset);
            w.f32(fo.frequency);
            // fo.name is derived by the reader, not stored.
        }
    }

    // rank distributions — max_rank, then for each non-empty partition its ion rows then noise.
    // `rank_dist` is already in reader order (increasing partition index; ions in block order,
    // noise stored last), so a straight iteration reproduces the stream the reader expects.
    w.i32(m.max_rank);
    for rd in &m.rank_dist {
        for (_name, freqs) in &rd.ions {
            for &f in freqs {
                w.f32(f);
            }
        }
    }

    // error / isotope distributions (present iff error_scaling_factor > 0), one per partition
    w.i32(m.error_scaling_factor);
    if m.error_scaling_factor > 0 {
        for ed in &m.error_dist {
            for &f in &ed.signal {
                w.f32(f);
            }
            for &f in &ed.noise {
                w.f32(f);
            }
            for &f in &ed.ion_existence {
                w.f32(f);
            }
        }
    }

    // sentinel — the reader validates this to prove the stream stayed aligned
    w.i32(crate::TERMINATOR);
    w.b
}

/// Encode a [`ScoringModel`] and write it to `path`.
pub fn write_param_file<P: AsRef<Path>>(path: P, m: &ScoringModel) -> io::Result<()> {
    std::fs::write(path, write_param(m))
}
