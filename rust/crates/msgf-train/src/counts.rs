//! The counting pass — corpus in, `ScoringModel` out.
//!
//! Training a fragment-scoring model is histogramming, not optimisation: every number in the
//! emitted `.param` is a ratio of counts, so a run is fully reproducible (no RNG, no learning
//! rate). The four trained sections and how each is counted:
//!
//! | Section | Counted as |
//! |---|---|
//! | §4 precursor offsets | fraction of spectra with a peak at a (charge-reduced) precursor offset, on the nominal-mass grid |
//! | §5 fragment offsets | matched sites / all sites, per candidate ion type per partition; kept above a frequency threshold |
//! | §6 rank distributions | per partition: `sites whose matched peak had rank r` ÷ `spectra`, with the `max_rank` bin holding absent sites |
//! | §7 error / isotope dists | main-ion edge mass error (signal) and ion-existence, over true cleavage pairs |
//!
//! **Signal vs. noise.** The `noise` row of §6 and the `Noise` row of §7 are counted the same way
//! as their signal counterparts, but at the node positions of a *decoy* peptide — a deterministic
//! shuffle of the identified peptide (C-terminal residue fixed, so the parent mass and partition
//! are unchanged). Decoy positions that collide with a true node position are skipped, so the
//! noise population is "what this spectrum offers a wrong peptide at the same mass".
//!
//! **Normalisation.** Rank rows are per *spectrum*, not per site: `row[r]` is the average number
//! of sites per spectrum whose matched peak had rank `r`, so a row sums to the average number of
//! scored sites per spectrum and the ratio `ion[r] / noise[r]` in
//! `ScoringModel::score_from_table` is a like-for-like likelihood ratio. Every bin is floored by
//! add-λ smoothing, which keeps `ln(ion/noise)` finite and makes never-observed ranks score ~0.

use crate::corpus::TrainingPsm;
use crate::ions::Candidate;
use crate::partition::PartitionScheme;
use crate::TrainConfig;
use msgf_chem::peptide::{self, Residue};
use msgf_chem::{round_half_up, scaling};
use msgf_scorer::preprocess::preprocess;
use msgf_scorer::scored_spectrum::{peak_by_mass, RankedPeak};
use msgf_scorer::{ErrorDist, FragOff, Partition, PrecursorOff, RankDist, ScoringModel};
use rayon::prelude::*;

/// `Composition.ChargeCarrierMass()` as MS-GF+'s preprocessing uses it.
const CHARGE_CARRIER: f32 = 1.00727649_f64 as f32;

/// Raw counts for one training run, indexed `[partition][candidate][bin]`.
pub struct Counts {
    n_cand: usize,
    n_rank: usize,
    n_err: usize,
    /// Cleavage sites evaluated (per partition, per candidate).
    sites: Vec<u64>,
    /// Sites whose theoretical m/z matched a peak.
    matched: Vec<u64>,
    /// Rank histogram; bin `max_rank` = "ion absent".
    rank: Vec<u32>,
    dsites: Vec<u64>,
    drank: Vec<u32>,
    /// Mass-error histogram over true edges, per candidate acting as the main ion.
    err_sig: Vec<u32>,
    /// Same over decoy edges.
    err_noise: Vec<u32>,
    /// The four (cur-present, prev-present) combinations over true edges.
    ion_exist: Vec<u64>,
    /// Spectra visited (for the report).
    pub spectra: u64,
}

impl Counts {
    fn new(n_part: usize, n_cand: usize, n_rank: usize, n_err: usize) -> Self {
        Self {
            n_cand,
            n_rank,
            n_err,
            sites: vec![0; n_part * n_cand],
            matched: vec![0; n_part * n_cand],
            rank: vec![0; n_part * n_cand * n_rank],
            dsites: vec![0; n_part * n_cand],
            drank: vec![0; n_part * n_cand * n_rank],
            err_sig: vec![0; n_part * n_cand * n_err],
            err_noise: vec![0; n_part * n_cand * n_err],
            ion_exist: vec![0; n_part * n_cand * 4],
            spectra: 0,
        }
    }

    #[inline]
    fn pc(&self, p: usize, c: usize) -> usize {
        p * self.n_cand + c
    }
    #[inline]
    fn rk(&self, p: usize, c: usize, b: usize) -> usize {
        (p * self.n_cand + c) * self.n_rank + b
    }
    #[inline]
    fn er(&self, p: usize, c: usize, b: usize) -> usize {
        (p * self.n_cand + c) * self.n_err + b
    }

    fn merge(mut self, o: Counts) -> Self {
        for (a, b) in self.sites.iter_mut().zip(&o.sites) {
            *a += b;
        }
        for (a, b) in self.matched.iter_mut().zip(&o.matched) {
            *a += b;
        }
        for (a, b) in self.rank.iter_mut().zip(&o.rank) {
            *a += b;
        }
        for (a, b) in self.dsites.iter_mut().zip(&o.dsites) {
            *a += b;
        }
        for (a, b) in self.drank.iter_mut().zip(&o.drank) {
            *a += b;
        }
        for (a, b) in self.err_sig.iter_mut().zip(&o.err_sig) {
            *a += b;
        }
        for (a, b) in self.err_noise.iter_mut().zip(&o.err_noise) {
            *a += b;
        }
        for (a, b) in self.ion_exist.iter_mut().zip(&o.ion_exist) {
            *a += b;
        }
        self.spectra += o.spectra;
        self
    }
}

/// Deterministic decoy peptide: shuffle all but the C-terminal residue (so the parent mass, and
/// therefore the partition, is unchanged). Seeded from the sequence — no RNG state, so two runs
/// over the same corpus produce byte-identical models.
fn decoy_residues(res: &[Residue]) -> Vec<Residue> {
    let mut out = res.to_vec();
    let n = out.len();
    if n < 4 {
        return out;
    }
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for r in res {
        h ^= r.aa as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    let mut s = h | 1;
    for i in (1..n - 1).rev() {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        out.swap(i, (s % (i as u64 + 1)) as usize);
    }
    out
}

/// Scratch buffers reused across PSMs (one per worker thread).
struct Scratch {
    peaks: Vec<RankedPeak>,
    /// Matched main-ion mass per `[candidate][site]`, `-1` when the ion is absent.
    mm: Vec<f32>,
    dmm: Vec<f32>,
    nominal: Vec<i32>,
    accurate: Vec<f64>,
    dnominal: Vec<i32>,
    daccurate: Vec<f64>,
    /// Sorted target node masses (nominal) per direction, for decoy-collision rejection.
    tgt_prefix: Vec<i32>,
    tgt_suffix: Vec<i32>,
}

impl Scratch {
    fn new() -> Self {
        Self {
            peaks: Vec::new(),
            mm: Vec::new(),
            dmm: Vec::new(),
            nominal: Vec::new(),
            accurate: Vec::new(),
            dnominal: Vec::new(),
            daccurate: Vec::new(),
            tgt_prefix: Vec::new(),
            tgt_suffix: Vec::new(),
        }
    }
}

/// Count one PSM into `acc`.
fn accumulate(
    psm: &TrainingPsm,
    cfg: &TrainConfig,
    scheme: &PartitionScheme,
    cands: &[Candidate],
    bootstrap: &ScoringModel,
    acc: &mut Counts,
    s: &mut Scratch,
) {
    s.peaks = preprocess(bootstrap, psm.charge, psm.parent_mass, &psm.peaks);
    if s.peaks.is_empty() {
        return;
    }
    s.nominal = peptide::nominal_prefix_masses(&psm.residues);
    s.accurate = peptide::accurate_prefix_masses(&psm.residues);
    let n = s.nominal.len();
    if n < 3 {
        return;
    }
    let sites = n - 1; // cleavage after residue 0..n-2
    let pep_nom = s.nominal[n - 1];

    let decoy = decoy_residues(&psm.residues);
    s.dnominal = peptide::nominal_prefix_masses(&decoy);
    s.daccurate = peptide::accurate_prefix_masses(&decoy);

    // Target node positions, for rejecting decoy positions that land on a true one.
    s.tgt_prefix.clear();
    s.tgt_suffix.clear();
    for j in 0..sites {
        s.tgt_prefix.push(s.nominal[j]);
        s.tgt_suffix.push(pep_nom - s.nominal[j]);
    }
    s.tgt_prefix.sort_unstable();
    s.tgt_suffix.sort_unstable();

    // The partition serving each mass segment is constant for the spectrum.
    let mut part_for_seg = [usize::MAX; 8];
    let nseg = cfg.num_segments.min(8);
    for seg in 0..nseg {
        part_for_seg[seg as usize] = scheme
            .index_of(psm.charge, psm.parent_mass, seg)
            .unwrap_or(usize::MAX);
    }
    if part_for_seg[..nseg as usize]
        .iter()
        .all(|&p| p == usize::MAX)
    {
        return;
    }
    acc.spectra += 1;

    let nc = cands.len();
    s.mm.clear();
    s.mm.resize(nc * sites, -1.0);
    s.dmm.clear();
    s.dmm.resize(nc * sites, -1.0);

    let max_rank = cfg.max_rank;
    let absent_bin = max_rank as usize;

    for (ci, cand) in cands.iter().enumerate() {
        for j in 0..sites {
            // ---- target site
            let node_nom = if cand.is_prefix {
                s.nominal[j]
            } else {
                pep_nom - s.nominal[j]
            };
            let node_mass = scaling::nominal_to_mass(node_nom);
            let theo = node_mass / cand.charge as f32 + cand.offset;
            if theo > 0.0 {
                let seg = scheme.segment_num(theo, psm.parent_mass);
                let part = part_for_seg[seg as usize];
                if part != usize::MAX {
                    let tol = cfg.mme.window_da(theo as f64) as f32;
                    let i = acc.pc(part, ci);
                    acc.sites[i] += 1;
                    match peak_by_mass(&s.peaks, theo, tol) {
                        Some(p) => {
                            acc.matched[i] += 1;
                            let bin = if p.rank > max_rank {
                                (max_rank - 1) as usize
                            } else {
                                (p.rank - 1) as usize
                            };
                            let k = acc.rk(part, ci, bin);
                            acc.rank[k] += 1;
                            s.mm[ci * sites + j] = (p.mz - cand.offset) * cand.charge as f32;
                        }
                        None => {
                            let k = acc.rk(part, ci, absent_bin);
                            acc.rank[k] += 1;
                        }
                    }
                }
            }

            // ---- decoy site (skip positions that coincide with a true node)
            let dnode_nom = if cand.is_prefix {
                s.dnominal[j]
            } else {
                pep_nom - s.dnominal[j]
            };
            // A decoy position must not sit on a true cleavage node of *either* series: a decoy
            // prefix mass can coincide with a target suffix mass, and that peak is a real y ion,
            // not noise.
            if collides(&s.tgt_prefix, dnode_nom) || collides(&s.tgt_suffix, dnode_nom) {
                continue;
            }
            let dnode_mass = scaling::nominal_to_mass(dnode_nom);
            let dtheo = dnode_mass / cand.charge as f32 + cand.offset;
            if dtheo <= 0.0 {
                continue;
            }
            let seg = scheme.segment_num(dtheo, psm.parent_mass);
            let part = part_for_seg[seg as usize];
            if part == usize::MAX {
                continue;
            }
            let tol = cfg.mme.window_da(dtheo as f64) as f32;
            let di = acc.pc(part, ci);
            acc.dsites[di] += 1;
            match peak_by_mass(&s.peaks, dtheo, tol) {
                Some(p) => {
                    let bin = if p.rank > max_rank {
                        (max_rank - 1) as usize
                    } else {
                        (p.rank - 1) as usize
                    };
                    let k = acc.rk(part, ci, bin);
                    acc.drank[k] += 1;
                    s.dmm[ci * sites + j] = (p.mz - cand.offset) * cand.charge as f32;
                }
                None => {
                    let k = acc.rk(part, ci, absent_bin);
                    acc.drank[k] += 1;
                }
            }
        }
    }

    // ---- edges: mass error + ion existence, per candidate acting as the main ion.
    // The scorer scores every edge in the *last* segment's partition, so they are counted there.
    let epart = part_for_seg[(nseg - 1) as usize];
    if epart == usize::MAX {
        return;
    }
    let esf = cfg.error_scaling_factor;
    for (ci, cand) in cands.iter().enumerate() {
        for j in 0..sites {
            for (mm, acc_masses, sig) in [(&s.mm, &s.accurate, true), (&s.dmm, &s.daccurate, false)]
            {
                let cur = mm[ci * sites + j];
                let (prev, theo) = if cand.is_prefix {
                    let prev = if j == 0 { 0.0 } else { mm[ci * sites + j - 1] };
                    let base = if j == 0 { 0.0 } else { acc_masses[j - 1] };
                    (prev, (acc_masses[j] - base) as f32)
                } else {
                    let prev = if j + 1 < sites {
                        mm[ci * sites + j + 1]
                    } else {
                        0.0
                    };
                    (prev, (acc_masses[j + 1] - acc_masses[j]) as f32)
                };
                let index = (cur >= 0.0) as usize + 2 * (prev >= 0.0) as usize;
                if sig {
                    acc.ion_exist[(epart * acc.n_cand + ci) * 4 + index] += 1;
                }
                if index == 3 {
                    let e = cur - prev - theo;
                    let bin = (round_half_up(e * esf as f32).clamp(-esf, esf) + esf) as usize;
                    if sig {
                        let k = acc.er(epart, ci, bin);
                        acc.err_sig[k] += 1;
                    } else {
                        let k = acc.er(epart, ci, bin);
                        acc.err_noise[k] += 1;
                    }
                }
            }
        }
    }
}

/// Is `x` within one nominal bin of any (sorted) true node position?
#[inline]
fn collides(sorted: &[i32], x: i32) -> bool {
    match sorted.binary_search(&x) {
        Ok(_) => true,
        Err(i) => {
            (i < sorted.len() && (sorted[i] - x).abs() <= 1)
                || (i > 0 && (x - sorted[i - 1]).abs() <= 1)
        }
    }
}

/// Run the counting sweep over the corpus (parallel over chunks; the result is order-independent
/// because every worker accumulates integer counts that are summed at the end).
pub fn sweep(
    psms: &[TrainingPsm],
    cfg: &TrainConfig,
    scheme: &PartitionScheme,
    cands: &[Candidate],
    bootstrap: &ScoringModel,
) -> Counts {
    let n_part = scheme.partitions.len();
    let n_cand = cands.len();
    let n_rank = (cfg.max_rank + 1) as usize;
    let n_err = (cfg.error_scaling_factor * 2 + 1) as usize;
    psms.par_chunks(512)
        .map(|chunk| {
            let mut acc = Counts::new(n_part, n_cand, n_rank, n_err);
            let mut s = Scratch::new();
            for psm in chunk {
                accumulate(psm, cfg, scheme, cands, bootstrap, &mut acc, &mut s);
            }
            acc
        })
        .reduce(
            || Counts::new(n_part, n_cand, n_rank, n_err),
            |a, b| a.merge(b),
        )
}

/// §4 — precursor offset frequencies, counted on the **raw** spectrum (the scorer filters
/// precursor peaks before deconvolution, so the offsets live in raw m/z space).
///
/// Offsets are enumerated on the nominal-mass grid: `offset = k / NOMINAL` for integer `k`, which
/// is the grid every real `.param` model's precursor offsets fall on.
pub fn precursor_offsets(psms: &[TrainingPsm], cfg: &TrainConfig) -> Vec<PrecursorOff> {
    let lo = cfg.precursor_offset_lo;
    let hi = cfg.precursor_offset_hi;
    let width = (hi - lo + 1) as usize;
    let max_charge = psms.iter().map(|p| p.charge).max().unwrap_or(0);
    let mut out = Vec::new();

    for charge in cfg.charge_min..=max_charge {
        let group: Vec<&TrainingPsm> = psms.iter().filter(|p| p.charge == charge).collect();
        if group.len() < cfg.min_psms_per_partition {
            continue;
        }
        let before = out.len();
        let mut at_precursor = 0.0f32;
        for rc in 0..charge {
            let cc = charge - rc;
            if cc == 0 {
                continue;
            }
            let hits: Vec<u64> = group
                .par_iter()
                .fold(
                    || vec![0u64; width],
                    |mut h, psm| {
                        let base = (psm.parent_mass + cc as f32 * CHARGE_CARRIER) / cc as f32;
                        for k in lo..=hi {
                            let mz = base + k as f32 / scaling::NOMINAL;
                            if mz <= 0.0 {
                                continue;
                            }
                            let tol = cfg.mme.window_da(mz as f64) as f32;
                            if has_peak(&psm.peaks, mz, tol) {
                                h[(k - lo) as usize] += 1;
                            }
                        }
                        h
                    },
                )
                .reduce(
                    || vec![0u64; width],
                    |mut a, b| {
                        for (x, y) in a.iter_mut().zip(&b) {
                            *x += y;
                        }
                        a
                    },
                );
            // An offset only means something if it stands *out of* the local peak density: for a
            // charge-reduced species the scanned window can sit inside the fragment-rich part of
            // the spectrum, where every offset "hits" and none of them is a precursor artefact.
            // So require both the absolute threshold and a clear margin over the window's median.
            let mut freqs: Vec<f32> = (lo..=hi)
                .map(|k| hits[(k - lo) as usize] as f32 / group.len() as f32)
                .collect();
            let mut sorted = freqs.clone();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let median = sorted[sorted.len() / 2];
            let floor = cfg
                .precursor_freq_threshold
                .max(median * cfg.precursor_contrast);
            let mut picked: Vec<(i32, f32)> = (lo..=hi)
                .zip(freqs.drain(..))
                .filter(|&(_, f)| f >= floor)
                .collect();
            picked.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
            picked.truncate(cfg.max_precursor_offsets_per_charge);
            picked.sort_by_key(|&(k, _)| k);
            if rc == 0 {
                at_precursor = hits[(0 - lo) as usize] as f32 / group.len() as f32;
            }
            for (k, f) in picked {
                out.push(PrecursorOff {
                    charge,
                    reduced_charge: rc,
                    offset: k as f32 / scaling::NOMINAL,
                    tol_ppm: matches!(cfg.mme.unit, msgf_chem::Unit::Ppm),
                    tol_val: cfg.mme.value as f32,
                    frequency: f,
                });
            }
        }
        if cfg.precursor_defaults && out.len() == before {
            out.extend(chemistry_precursor_offsets(charge, at_precursor, cfg));
        }
    }
    out
}

/// Precursor offsets that follow from chemistry rather than from counting: the precursor peak
/// itself (with a one-bin margin either side) and the precursor water loss.
///
/// These exist because a *library* corpus cannot teach them — consensus library spectra have had
/// the precursor region deleted before deposition, so the counted frequency is 0 and the trained
/// model would filter nothing, while a real query spectrum still carries its precursor. Only the
/// `(reduced_charge, offset)` pair is used at scoring time (`preprocess` ignores the frequency),
/// so this adds a filtering *rule*, not a trained number; the frequency written is the one we
/// actually measured, whatever it was.
fn chemistry_precursor_offsets(charge: i32, measured: f32, cfg: &TrainConfig) -> Vec<PrecursorOff> {
    let mut bins = vec![-1i32, 0, 1];
    // Water loss, in nominal bins of the charge-reduced m/z; take both neighbours when it falls
    // between bins so the ±mme window covers it either way.
    let w = -(msgf_chem::mass::WATER as f32) / charge as f32 * scaling::NOMINAL;
    bins.push(w.floor() as i32);
    bins.push(w.ceil() as i32);
    bins.sort_unstable();
    bins.dedup();
    bins.into_iter()
        .map(|k| PrecursorOff {
            charge,
            reduced_charge: 0,
            offset: k as f32 / scaling::NOMINAL,
            tol_ppm: matches!(cfg.mme.unit, msgf_chem::Unit::Ppm),
            tol_val: cfg.mme.value as f32,
            frequency: measured,
        })
        .collect()
}

/// Any raw peak within `±tol` of `mz` (raw peaks are m/z-sorted).
fn has_peak(peaks: &[(f32, f32)], mz: f32, tol: f32) -> bool {
    let lo = mz - tol;
    let i = peaks.partition_point(|p| p.0 < lo);
    i < peaks.len() && peaks[i].0 <= mz + tol
}

/// What a training run selected, for the run report.
#[derive(Debug, Clone)]
pub struct PartitionReport {
    pub charge: i32,
    pub parent_mass: f32,
    pub seg: i32,
    pub spectra: u64,
    /// `(label, name, frequency)` of every ion type kept.
    pub ions: Vec<(String, String, f32)>,
    pub main_ion: Option<String>,
}

/// Assemble the counted tables into a `ScoringModel`.
pub fn build_model(
    counts: &Counts,
    cfg: &TrainConfig,
    scheme: &PartitionScheme,
    cands: &[Candidate],
    precursor_off: Vec<PrecursorOff>,
    charge_histogram: Vec<(i32, i32)>,
) -> (ScoringModel, Vec<PartitionReport>) {
    let n_rank = (cfg.max_rank + 1) as usize;
    let n_err = (cfg.error_scaling_factor * 2 + 1) as usize;
    let lambda = cfg.smoothing;

    // ---- §5: pick the ion types worth scoring in each partition.
    let mut selected: Vec<Vec<usize>> = Vec::with_capacity(scheme.partitions.len());
    let mut frag_off: Vec<Vec<FragOff>> = Vec::with_capacity(scheme.partitions.len());
    for p in 0..scheme.partitions.len() {
        let observed: Vec<(usize, f32)> = (0..cands.len())
            .filter_map(|ci| {
                let i = counts.pc(p, ci);
                let sites = counts.sites[i];
                (sites > 0).then(|| (ci, counts.matched[i] as f32 / sites as f32))
            })
            .collect();
        let mut scored: Vec<(usize, f32)> = observed
            .iter()
            .copied()
            .filter(|&(_, f)| f >= cfg.ion_freq_threshold)
            .collect();
        // A partition that clears the threshold with nothing would contribute a flat 0 to every
        // node it serves — worse than scoring its weak best ion, and it costs the spectrum a whole
        // m/z segment. So a populated partition always keeps at least its most frequent ion type.
        if scored.is_empty() {
            if let Some(best) = observed
                .iter()
                .copied()
                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            {
                scored.push(best);
            }
        }
        scored.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        // Two candidates can round to the same `.param` ion name (charge-2 losses do); the name is
        // the model's lookup key, so keep only the more frequent of any such pair.
        let mut names: Vec<&str> = Vec::new();
        scored.retain(|(ci, _)| {
            let n = cands[*ci].name.as_str();
            if names.contains(&n) {
                false
            } else {
                names.push(n);
                true
            }
        });
        scored.truncate(cfg.max_ions_per_partition);
        frag_off.push(
            scored
                .iter()
                .map(|(ci, f)| FragOff {
                    is_prefix: cands[*ci].is_prefix,
                    charge: cands[*ci].charge,
                    offset: cands[*ci].offset,
                    frequency: *f,
                    name: cands[*ci].name.clone(),
                })
                .collect(),
        );
        selected.push(scored.into_iter().map(|(ci, _)| ci).collect());
    }

    // ---- §6: rank distributions (ion rows + the shared decoy-derived noise row).
    let mut rank_dist = Vec::new();
    for (p, sel) in selected.iter().enumerate() {
        if sel.is_empty() {
            continue; // §6 skips partitions with no scored ion type
        }
        let n_spec = scheme.group_spectra[p].max(1) as f64;
        let mut ions = Vec::with_capacity(sel.len() + 1);
        for &ci in sel {
            let raw: Vec<f64> = (0..n_rank)
                .map(|b| counts.rank[counts.rk(p, ci, b)] as f64)
                .collect();
            let row = rank_row(&raw, cfg.rank_smoothing, lambda, n_spec);
            ions.push((cands[ci].name.clone(), row));
        }
        // Noise: decoy-position ranks, averaged over the singly-charged scored ion types so the
        // row keeps one ion-type's worth of sites (the score divides it by min(charge, segments)).
        let noise_src: Vec<usize> = {
            let ones: Vec<usize> = sel
                .iter()
                .copied()
                .filter(|&ci| cands[ci].charge == 1)
                .collect();
            if ones.is_empty() {
                sel.clone()
            } else {
                ones
            }
        };
        let k = noise_src.len() as f64;
        let raw: Vec<f64> = (0..n_rank)
            .map(|b| {
                let sum: u64 = noise_src
                    .iter()
                    .map(|&ci| counts.drank[counts.rk(p, ci, b)] as u64)
                    .sum();
                sum as f64 / k
            })
            .collect();
        let noise = rank_row(&raw, cfg.rank_smoothing, lambda, n_spec);
        ions.push(("noise".to_string(), noise));
        rank_dist.push(RankDist {
            partition_index: p,
            ions,
        });
    }

    // ---- §7: error distributions. Edges are scored in the last segment's partition, so each
    // (charge, parent_mass) group's counts live there; every segment of the group gets a copy.
    let last_seg = cfg.num_segments - 1;
    let mut error_dist = Vec::with_capacity(scheme.partitions.len());
    let mut main_ions: Vec<Option<usize>> = vec![None; scheme.partitions.len()];
    for (p, part) in scheme.partitions.iter().enumerate() {
        let src = scheme
            .partitions
            .iter()
            .position(|q| {
                q.charge == part.charge && q.parent_mass == part.parent_mass && q.seg == last_seg
            })
            .unwrap_or(p);
        // `NewScoredSpectrum.determineIonTypes`: the main ion is the one with the greatest
        // frequency *summed over the group's segments* — not the top ion of one segment. The
        // error distribution describes that ion's edges, so it must be chosen the same way.
        let main = {
            let mut best: Option<(usize, f32)> = None;
            for (q, qpart) in scheme.partitions.iter().enumerate() {
                if qpart.charge != part.charge || qpart.parent_mass != part.parent_mass {
                    continue;
                }
                for &ci in &selected[q] {
                    let f: f32 = scheme
                        .partitions
                        .iter()
                        .enumerate()
                        .filter(|(_, r)| {
                            r.charge == part.charge && r.parent_mass == part.parent_mass
                        })
                        .map(|(r, _)| {
                            frag_off[r]
                                .iter()
                                .find(|fo| fo.name == cands[ci].name)
                                .map(|fo| fo.frequency)
                                .unwrap_or(0.0)
                        })
                        .sum();
                    if best.is_none_or(|(_, bf)| f > bf) {
                        best = Some((ci, f));
                    }
                }
            }
            best.map(|(ci, _)| ci)
        };
        main_ions[p] = main;
        let (signal, noise, ion_existence) = match main {
            Some(ci) => {
                let sig: Vec<u64> = (0..n_err)
                    .map(|b| counts.err_sig[counts.er(src, ci, b)] as u64)
                    .collect();
                let noi: Vec<u64> = (0..n_err)
                    .map(|b| counts.err_noise[counts.er(src, ci, b)] as u64)
                    .collect();
                let ie: Vec<u64> = (0..4)
                    .map(|b| counts.ion_exist[(src * counts.n_cand + ci) * 4 + b])
                    .collect();
                (
                    normalize(&sig, cfg.error_smoothing),
                    normalize(&noi, cfg.error_smoothing),
                    normalize(&ie, 0.5),
                )
            }
            None => (
                vec![1.0 / n_err as f32; n_err],
                vec![1.0 / n_err as f32; n_err],
                vec![0.25; 4],
            ),
        };
        error_dist.push(ErrorDist {
            signal,
            noise,
            ion_existence: [
                ion_existence[0].max(1e-4),
                ion_existence[1].max(1e-4),
                ion_existence[2].max(1e-4),
                ion_existence[3].max(1e-4),
            ],
        });
    }

    let report: Vec<PartitionReport> = scheme
        .partitions
        .iter()
        .enumerate()
        .map(|(p, part)| PartitionReport {
            charge: part.charge,
            parent_mass: part.parent_mass,
            seg: part.seg,
            spectra: scheme.group_spectra[p],
            ions: selected[p]
                .iter()
                .map(|&ci| {
                    (
                        cands[ci].label.clone(),
                        cands[ci].name.clone(),
                        frag_off[p]
                            .iter()
                            .find(|f| f.name == cands[ci].name)
                            .map(|f| f.frequency)
                            .unwrap_or(0.0),
                    )
                })
                .collect(),
            main_ion: main_ions[p].map(|ci| cands[ci].label.clone()),
        })
        .collect();

    let model = ScoringModel {
        version: cfg.version,
        activation: cfg.activation.clone(),
        instrument: cfg.instrument.clone(),
        enzyme: cfg.enzyme.clone(),
        protocol: cfg.protocol.clone(),
        mme: cfg.mme,
        apply_deconvolution: cfg.apply_deconvolution,
        deconvolution_error_tolerance: cfg.deconvolution_error_tolerance,
        charge_histogram,
        num_segments: cfg.num_segments,
        partitions: scheme.partitions.clone(),
        precursor_off,
        frag_off,
        max_rank: cfg.max_rank,
        rank_dist,
        error_scaling_factor: cfg.error_scaling_factor,
        error_dist,
    };
    (model, report)
}

/// Turn a rank histogram into a `.param` rank row.
///
/// The counts thin out badly at high rank — few spectra even *have* a 120th peak — and the score
/// is a ratio of two such bins, so raw counts there produce large scores out of one or two
/// observations. Rank distributions are smooth in rank, so neighbouring ranks are pooled with a
/// window proportional to the rank (±`frac·r`): dense low ranks keep their resolution, sparse high
/// ranks get averaged. Both the ion rows and the noise row are pooled identically, so a region
/// where they agree still scores 0. The "ion absent" bin (last) is never pooled — it is not a
/// rank, and it is always well populated.
///
/// The row is then per-*spectrum*: `row[r]` = average number of sites per spectrum whose matched
/// peak had rank `r`, with add-λ keeping every bin non-zero.
fn rank_row(counts: &[f64], frac: f64, lambda: f64, n_spectra: f64) -> Vec<f32> {
    let n = counts.len();
    let present = n - 1;
    let mut out = Vec::with_capacity(n);
    for b in 0..present {
        let w = (b as f64 * frac).round() as usize;
        let lo = b.saturating_sub(w);
        let hi = (b + w + 1).min(present);
        let mean = counts[lo..hi].iter().sum::<f64>() / (hi - lo) as f64;
        out.push(((mean + lambda) / n_spectra) as f32);
    }
    out.push(((counts[present] + lambda) / n_spectra) as f32);
    out
}

/// Add-λ smoothed probability vector (sums to 1; never zero, so `ln` of a ratio stays finite).
fn normalize(counts: &[u64], lambda: f64) -> Vec<f32> {
    let total: f64 = counts.iter().map(|&c| c as f64).sum::<f64>() + lambda * counts.len() as f64;
    if total <= 0.0 {
        return vec![1.0 / counts.len() as f32; counts.len()];
    }
    counts
        .iter()
        .map(|&c| ((c as f64 + lambda) / total) as f32)
        .collect()
}

/// The minimal model `preprocess` needs before anything is trained: tolerances, deconvolution
/// settings and (once counted) the precursor offsets it filters on.
pub fn bootstrap_model(cfg: &TrainConfig, precursor_off: Vec<PrecursorOff>) -> ScoringModel {
    ScoringModel {
        version: cfg.version,
        activation: cfg.activation.clone(),
        instrument: cfg.instrument.clone(),
        enzyme: cfg.enzyme.clone(),
        protocol: cfg.protocol.clone(),
        mme: cfg.mme,
        apply_deconvolution: cfg.apply_deconvolution,
        deconvolution_error_tolerance: cfg.deconvolution_error_tolerance,
        charge_histogram: Vec::new(),
        num_segments: cfg.num_segments,
        partitions: Vec::new(),
        precursor_off,
        frag_off: Vec::new(),
        max_rank: cfg.max_rank,
        rank_dist: Vec::new(),
        error_scaling_factor: cfg.error_scaling_factor,
        error_dist: Vec::new(),
    }
}

/// Convenience: the whole pipeline, corpus → model.
pub fn train(
    psms: &[TrainingPsm],
    cfg: &TrainConfig,
) -> (ScoringModel, Vec<PartitionReport>, PartitionScheme, u64) {
    let precursor_off = precursor_offsets(psms, cfg);
    let bootstrap = bootstrap_model(cfg, precursor_off.clone());
    let scheme = PartitionScheme::build(
        psms,
        cfg.num_segments,
        cfg.min_psms_per_partition,
        cfg.max_partitions_per_charge,
    );
    let cands = crate::ions::candidates(cfg.max_fragment_charge);
    let counts = sweep(psms, cfg, &scheme, &cands, &bootstrap);

    let mut hist: Vec<(i32, i32)> = Vec::new();
    for p in psms {
        match hist.iter_mut().find(|(c, _)| *c == p.charge) {
            Some(e) => e.1 += 1,
            None => hist.push((p.charge, 1)),
        }
    }
    hist.sort_by_key(|(c, _)| *c);

    let spectra = counts.spectra;
    let (model, report) = build_model(&counts, cfg, &scheme, &cands, precursor_off, hist);
    (model, report, scheme, spectra)
}

/// Partitions in the model whose `Partition` list is needed by callers building reports.
pub fn partitions_of(model: &ScoringModel) -> &[Partition] {
    &model.partitions
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decoy_keeps_mass_and_cterm() {
        let res = peptide::parse("SAMPLERPEPTIDEK").unwrap();
        let d = decoy_residues(&res);
        assert_eq!(d.len(), res.len());
        assert_eq!(d.last().unwrap().aa, b'K');
        let m0: f64 = res
            .iter()
            .map(|r| msgf_chem::residue_mass(r.aa).unwrap())
            .sum();
        let m1: f64 = d
            .iter()
            .map(|r| msgf_chem::residue_mass(r.aa).unwrap())
            .sum();
        assert!((m0 - m1).abs() < 1e-9);
        assert_ne!(
            res.iter().map(|r| r.aa).collect::<Vec<_>>(),
            d.iter().map(|r| r.aa).collect::<Vec<_>>()
        );
    }

    #[test]
    fn collision_window_is_one_nominal_bin() {
        let sorted = vec![100, 200, 300];
        assert!(collides(&sorted, 200));
        assert!(collides(&sorted, 201));
        assert!(collides(&sorted, 199));
        assert!(!collides(&sorted, 202));
    }

    #[test]
    fn normalize_sums_to_one_and_floors_zeros() {
        let p = normalize(&[10, 0, 30], 0.5);
        assert!((p.iter().sum::<f32>() - 1.0).abs() < 1e-6);
        assert!(p[1] > 0.0);
    }
}
