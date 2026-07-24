//! Encoder for the MS-GF+ `.param` scoring-model format — the inverse of [`crate::read_param`].
//!
//! # Why this exists
//!
//! [`crate::read_param`] can only *consume* the trained `.param` models shipped with MS-GF+, which
//! are Copyright UC Regents under a non-commercial license (see `validation/README.md`,
//! `docs/models.md`). To ship a permissively-licensed MSGF_Rust we must be able to **produce** a
//! model of our own — trained from openly-licensed data (e.g. MassIVE-KB, CC0). This module is the
//! first half of that: it serialises any in-memory [`ScoringModel`] back into the on-disk byte
//! format, so a future `msgf-train` crate can emit a real `.param` that the existing reader,
//! scorer, and generating function accept unchanged.
//!
//! # Licensing boundary
//!
//! This encoder is **clean-room**: it is written from the format documented in `docs/param-format.md`
//! and the structure of [`crate::read_param`] in this repo — *not* transcribed from MS-GF+'s Java
//! `NewRankScorer.writeParameters`. The `.param` *file format* is an interface (uncopyrightable);
//! the encumbered artifacts are the trained *numbers* in UC's shipped `.param` files, which live
//! only under `validation/data/` (gitignored, fetched, test-only) and are never vendored here. A
//! [`ScoringModel`] we construct ourselves and write with this module carries no upstream license.
//!
//! # Fidelity
//!
//! `read_param(write_param(m)) == m` for every model produced by [`crate::read_param`] — the
//! round-trip is validated against all four high-res UC models in `tests/roundtrip_write.rs`. The
//! encoding mirrors Java `DataOutputStream`: big-endian scalars, `writeByte(len)` + UTF-16BE chars
//! for strings, per-partition parallel arrays in the reader's sorted order, and the trailing
//! `0x7FFFFFFF` (`Integer.MAX_VALUE`) sentinel.

use crate::ScoringModel;
use msgf_chem::Unit;
use std::io;
use std::path::Path;

/// Big-endian sink mirroring Java `DataOutputStream` (the write-side of [`crate::read_param`]'s
/// `Reader`).
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
    /// Java `writeByte(len)` + `writeChars` (len UTF-16BE code units). Inverse of `Reader::jstring`.
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

    // partitions (already in the reader's TreeSet order: charge, seg, parent_mass)
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
