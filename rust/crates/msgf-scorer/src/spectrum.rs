//! Spectrum preparation (peak order, precursor suppression, ranks, isotope-cluster reduction) and
//! the per-spectrum scoring context: segment partitions, site scores, anchor masses, edge terms.

use super::model::ScoreModel;
use super::{
    nominal_to_real, round_nominal, round_score, CARRIER, ISO_STEP1, ISO_STEP2, NOMINAL_SCALE,
    PROTON, WATER,
};
use std::cmp::Ordering;

/// One MS/MS spectrum at one charge, as input.
#[derive(Clone, Copy, Debug)]
pub struct Ms2<'a> {
    /// `(m/z, intensity)` as read from the source.
    pub peaks: &'a [(f64, f64)],
    pub precursor_mz: f64,
    pub charge: i32,
}

/// The fixed anchor ion type used by the edge term.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Anchor {
    pub prefix: bool,
    pub charge: i32,
    pub offset: f32,
}

/// A spectrum prepared for scoring at one charge. Immutable; cheap to share across threads.
#[derive(Clone, Debug)]
pub struct PreparedSpectrum<'m> {
    pub(crate) model: &'m ScoreModel,
    pub(crate) charge: i32,
    pub(crate) parent_mass: f32,
    pub(crate) n0: i32,
    /// Final peak list: m/z, intensity, rank; sorted by m/z.
    mz: Vec<f32>,
    intensity: Vec<f32>,
    rank: Vec<u32>,
    /// Partition index per segment.
    pub(crate) seg_partition: Vec<usize>,
    pub(crate) rho: f32,
    pub(crate) anchor: Option<Anchor>,
    /// R(eps_i) for i = 0, 1, 2 (neither / current-only / previous-only anchored).
    pub(crate) edge_const: [i32; 3],
    /// R(eps_3 + eta_j) for j in [0, 2E].
    pub(crate) edge_both: Vec<i32>,
    /// Pre(k), Suf(k), alpha(k) for k in [0, cache_max].
    pre: Vec<f32>,
    suf: Vec<f32>,
    alpha: Vec<f32>,
}

fn cmp_f32(a: f32, b: f32) -> Ordering {
    a.partial_cmp(&b).unwrap_or(Ordering::Equal)
}

/// Prepare a spectrum at its charge, caching node tables up to `N0 + 2`. `None` if the candidate
/// nominal mass `N0` is outside `[50, 10000]` (the spectrum is not scored at this charge).
pub fn prepare<'m>(model: &'m ScoreModel, ms2: &Ms2) -> Option<PreparedSpectrum<'m>> {
    prepare_with_cache(model, ms2, 2)
}

/// As [`prepare`], caching the candidate-independent node tables over `[0, N0 + above_n0]`
/// (use `-isotope_lo` when the isotope range goes below 0); other nodes are computed on demand.
pub fn prepare_with_cache<'m>(
    model: &'m ScoreModel,
    ms2: &Ms2,
    above_n0: i32,
) -> Option<PreparedSpectrum<'m>> {
    let z = ms2.charge;
    let zf = z as f32;
    let parent_mass = (ms2.precursor_mz as f32 * zf) - (zf * PROTON as f32);
    let n0 = round_nominal((parent_mass - WATER as f32) * NOMINAL_SCALE);
    if !(50..=10000).contains(&n0) {
        return None;
    }

    let peaks: Vec<(f32, f32)> = ms2
        .peaks
        .iter()
        .map(|&(m, i)| (m as f32, i as f32))
        .collect();
    let (mz, intensity, rank) = preprocess_peaks(&model.preprocess_view(), z, parent_mass, peaks);
    let n = mz.len();

    let segs = model.segments;
    let seg_partition: Vec<usize> = (0..segs)
        .map(|s| model.partition_for(z, s, parent_mass))
        .collect::<Option<_>>()?;

    // Peak density and edge tables from the last segment's partition.
    let rho =
        (n.max(1) as f32) / ((parent_mass - WATER as f32) / (model.tol_value * 2.0f32)).max(1.0f32);
    let edge_part = &model.partitions[*seg_partition.last().unwrap()];
    let one_m = 1.0f32 - rho;
    let q = [one_m * one_m, rho * one_m, rho * one_m, rho * rho];
    let eps: Vec<f32> = (0..4)
        .map(|i| (edge_part.existence[i] as f64 / q[i] as f64).ln() as f32)
        .collect();
    let edge_const = [
        round_score(eps[0]),
        round_score(eps[1]),
        round_score(eps[2]),
    ];
    let edge_both = edge_part
        .error_ll
        .iter()
        .map(|&eta| round_score(eps[3] + eta))
        .collect();

    // Anchor type: largest summed frequency over the segment partitions, grouped by
    // (derived name, exact offset bits); ties -> first appearance.
    let mut groups: Vec<((bool, i32, i64, u32), f32, Anchor)> = Vec::new();
    for &pi in &seg_partition {
        for ion in &model.partitions[pi].ions {
            let key = (
                ion.prefix,
                ion.charge,
                (ion.offset + 0.5f32).floor() as i64,
                ion.offset.to_bits(),
            );
            match groups.iter_mut().find(|g| g.0 == key) {
                Some(g) => g.1 += ion.frequency,
                None => groups.push((
                    key,
                    ion.frequency,
                    Anchor {
                        prefix: ion.prefix,
                        charge: ion.charge,
                        offset: ion.offset,
                    },
                )),
            }
        }
    }
    let mut anchor: Option<(f32, Anchor)> = None;
    for g in &groups {
        if anchor.map_or(true, |(best, _)| g.1 > best) {
            anchor = Some((g.1, g.2));
        }
    }

    let mut prep = PreparedSpectrum {
        model,
        charge: z,
        parent_mass,
        n0,
        mz,
        intensity,
        rank,
        seg_partition,
        rho,
        anchor: anchor.map(|a| a.1),
        edge_const,
        edge_both,
        pre: Vec::new(),
        suf: Vec::new(),
        alpha: Vec::new(),
    };
    let top = (n0 + above_n0).max(0);
    prep.pre = prep.site_table(top, true);
    prep.suf = prep.site_table(top, false);
    prep.alpha = prep.alpha_table(top);
    Some(prep)
}

/// The model settings spectrum preparation reads: fragment tolerance, isotope-cluster reduction
/// and the precursor-offset windows. Borrowed from a [`ScoreModel`] or a [`PreprocessParams`].
pub(crate) struct PreprocessView<'a> {
    pub tol_is_ppm: bool,
    pub tol_value: f32,
    pub reduce_isotopes: bool,
    pub reduce_tol: f32,
    pub precursor_offsets: &'a [(i32, i32, f32)],
}

impl PreprocessView<'_> {
    /// Fragment tolerance (Da) at theoretical m/z `x` (same arithmetic as the scoring model's).
    #[inline]
    fn tolerance_at(&self, x: f32) -> f32 {
        if self.tol_is_ppm {
            (x as f64 * self.tol_value as f64 * 1e-6) as f32
        } else {
            self.tol_value
        }
    }
}

/// Spectrum-preparation settings on their own, for callers that have no scoring tables yet (the
/// trainer prepares its corpus with a model it has not finished counting).
#[derive(Clone, Debug, PartialEq)]
pub struct PreprocessParams {
    pub tol_is_ppm: bool,
    pub tol_value: f32,
    pub reduce_isotopes: bool,
    pub reduce_tol: f32,
    /// `(charge, reduced charge, offset)` of each precursor-offset entry, in file order.
    pub precursor_offsets: Vec<(i32, i32, f32)>,
}

impl PreprocessParams {
    /// The preparation settings of a `.param` model.
    pub fn from_param(m: &crate::ScoringModel) -> PreprocessParams {
        PreprocessParams {
            tol_is_ppm: m.mme.unit == msgf_chem::Unit::Ppm,
            tol_value: m.mme.value as f32,
            reduce_isotopes: m.apply_deconvolution,
            reduce_tol: m.deconvolution_error_tolerance,
            precursor_offsets: m
                .precursor_off
                .iter()
                .map(|p| (p.charge, p.reduced_charge, p.offset))
                .collect(),
        }
    }

    fn view(&self) -> PreprocessView<'_> {
        PreprocessView {
            tol_is_ppm: self.tol_is_ppm,
            tol_value: self.tol_value,
            reduce_isotopes: self.reduce_isotopes,
            reduce_tol: self.reduce_tol,
            precursor_offsets: &self.precursor_offsets,
        }
    }
}

impl ScoreModel {
    pub(crate) fn preprocess_view(&self) -> PreprocessView<'_> {
        PreprocessView {
            tol_is_ppm: self.tol_is_ppm,
            tol_value: self.tol_value,
            reduce_isotopes: self.reduce_isotopes,
            reduce_tol: self.reduce_tol,
            precursor_offsets: &self.precursor_offsets,
        }
    }

    /// This model's spectrum-preparation settings.
    pub fn preprocess_params(&self) -> PreprocessParams {
        PreprocessParams {
            tol_is_ppm: self.tol_is_ppm,
            tol_value: self.tol_value,
            reduce_isotopes: self.reduce_isotopes,
            reduce_tol: self.reduce_tol,
            precursor_offsets: self.precursor_offsets.clone(),
        }
    }
}

/// One peak of a prepared spectrum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RankedPeak {
    pub mz: f32,
    pub intensity: f32,
    /// Intensity rank, 1 = most intense.
    pub rank: i32,
}

/// Prepare a peak list (`(m/z, intensity)` as f32) at charge `charge` and neutral parent mass
/// `parent_mass`: the same peak order, precursor suppression, ranking and isotope-cluster
/// reduction [`prepare`] applies, returned as the final m/z-sorted list.
pub fn preprocess(
    params: &PreprocessParams,
    charge: i32,
    parent_mass: f32,
    peaks: &[(f32, f32)],
) -> Vec<RankedPeak> {
    let (mz, intensity, rank) =
        preprocess_peaks(&params.view(), charge, parent_mass, peaks.to_vec());
    (0..mz.len())
        .map(|i| RankedPeak {
            mz: mz[i],
            intensity: intensity[i],
            rank: rank[i] as i32,
        })
        .collect()
}

/// The matched peak for theoretical m/z `x` with tolerance `tol` (Da) in an m/z-sorted list:
/// greatest intensity in the inclusive window `[x - tol, x + tol]`, ties to the last (highest
/// m/z) — the lookup rule [`PreparedSpectrum`] scores with.
pub fn peak_by_mass(peaks: &[RankedPeak], x: f32, tol: f32) -> Option<&RankedPeak> {
    let (lo, hi) = (x - tol, x + tol);
    let mut i = peaks.partition_point(|p| p.mz < lo);
    let mut best: Option<usize> = None;
    while i < peaks.len() && peaks[i].mz <= hi {
        if best.map_or(true, |b| peaks[i].intensity >= peaks[b].intensity) {
            best = Some(i);
        }
        i += 1;
    }
    best.map(|b| &peaks[b])
}

/// Peak order, precursor suppression, ranks and isotope-cluster reduction; returns the final
/// list (m/z, intensity, rank) sorted by m/z.
fn preprocess_peaks(
    model: &PreprocessView,
    z: i32,
    parent_mass: f32,
    mut peaks: Vec<(f32, f32)>,
) -> (Vec<f32>, Vec<f32>, Vec<u32>) {
    // Working order: input order if already non-decreasing in m/z, else stable (m/z, intensity).
    if peaks.windows(2).any(|w| w[1].0 < w[0].0) {
        peaks.sort_by(|a, b| cmp_f32(a.0, b.0).then(cmp_f32(a.1, b.1)));
    }
    let mut mz: Vec<f32> = peaks.iter().map(|p| p.0).collect();
    let mut inten: Vec<f32> = peaks.iter().map(|p| p.1).collect();

    suppress_precursor(model, z, parent_mass, &mz, &mut inten);

    // Ranks: intensity desc, then m/z desc, then working order.
    let n = mz.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        cmp_f32(inten[b], inten[a])
            .then(cmp_f32(mz[b], mz[a]))
            .then(a.cmp(&b))
    });
    let mut rank = vec![0u32; n];
    for (r, &i) in order.iter().enumerate() {
        rank[i] = r as u32 + 1;
    }

    if model.reduce_isotopes {
        reduce_isotope_clusters(&mut mz, z, model.reduce_tol);
    }
    // Final list: stable sort by (m/z, intensity) (the no-reduction path is already sorted, and
    // a stable sort leaves it unchanged).
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| cmp_f32(mz[a], mz[b]).then(cmp_f32(inten[a], inten[b])));
    (
        idx.iter().map(|&i| mz[i]).collect(),
        idx.iter().map(|&i| inten[i]).collect(),
        idx.iter().map(|&i| rank[i]).collect(),
    )
}

/// Zero every peak inside a precursor-offset window (inclusive).
fn suppress_precursor(model: &PreprocessView, z: i32, mass: f32, mz: &[f32], inten: &mut [f32]) {
    let offs = model.precursor_offsets;
    if offs.is_empty() {
        return;
    }
    let key = offs
        .iter()
        .map(|o| o.0)
        .filter(|&c| c <= z)
        .max()
        .unwrap_or_else(|| offs.iter().map(|o| o.0).min().unwrap());
    for &(oz, reduced, delta) in offs {
        if oz != key {
            continue;
        }
        let c = z - reduced;
        if c == 0 {
            continue;
        }
        let cf = c as f32;
        let x = ((mass + cf * CARRIER as f32) / cf) + delta;
        let t = model.tolerance_at(x);
        let (lo, hi) = (x - t, x + t);
        for (m, i) in mz.iter().zip(inten.iter_mut()) {
            if *m >= lo && *m <= hi {
                *i = 0.0;
            }
        }
    }
}

/// Rewrite multiply-charged isotope clusters to singly-charged m/z, in place (working order).
fn reduce_isotope_clusters(a: &mut [f32], z: i32, tol: f32) {
    let n = a.len();
    let carrier = CARRIER as f32;
    let step1 = ISO_STEP1 as f32;
    let step2 = ISO_STEP2 as f32;
    let mut consumed = vec![false; n];
    let rewrite = |y: f32, u: f32| u * y - (u - 1.0f32) * carrier;
    for i in 0..n {
        if consumed[i] {
            continue;
        }
        let base = a[i];
        let mut u = 2;
        'charges: while u < z && u < 4 {
            let uf = u as f32;
            let s1 = step1 / uf;
            for j in i + 1..n {
                let d = (a[j] - base) - s1;
                if -tol < d && d < tol {
                    consumed[j] = true;
                    a[i] = rewrite(a[i], uf);
                    let base2 = a[j];
                    let s2 = step2 / uf;
                    for k in j + 1..n {
                        let d2 = (a[k] - base2) - s2;
                        if -tol < d2 && d2 < tol {
                            consumed[k] = true;
                            a[k] = rewrite(a[k], uf);
                            break;
                        } else if d2 > tol {
                            break;
                        }
                    }
                    a[j] = rewrite(a[j], uf);
                    break 'charges;
                } else if d > tol {
                    break;
                }
            }
            u += 1;
        }
    }
}

impl<'m> PreparedSpectrum<'m> {
    /// The final peak list: (m/z, intensity, rank), sorted by m/z.
    pub fn peaks(&self) -> impl Iterator<Item = (f32, f32, u32)> + '_ {
        (0..self.mz.len()).map(move |i| (self.mz[i], self.intensity[i], self.rank[i]))
    }
    /// Candidate nominal mass N0.
    pub fn n0(&self) -> i32 {
        self.n0
    }
    /// Neutral parent mass M (f32).
    pub fn parent_mass(&self) -> f32 {
        self.parent_mass
    }
    pub fn charge(&self) -> i32 {
        self.charge
    }
    /// Peak density rho.
    pub fn rho(&self) -> f32 {
        self.rho
    }
    /// Sorted-partition index chosen for each segment.
    pub fn segment_partitions(&self) -> &[usize] {
        &self.seg_partition
    }
    /// Whether the anchor ion type is a prefix type (`None` if there is no anchor type).
    pub fn anchor_is_prefix(&self) -> Option<bool> {
        self.anchor.map(|a| a.prefix)
    }
    /// The three constant edge scores (neither, current-only, previous-only anchored).
    pub fn edge_constants(&self) -> [i32; 3] {
        self.edge_const
    }

    /// Index of the matched peak for theoretical m/z `x`: greatest intensity in the inclusive
    /// tolerance window, ties to the last (highest m/z).
    #[inline]
    fn lookup(&self, x: f32) -> Option<usize> {
        let t = self.model.tolerance_at(x);
        let (lo, hi) = (x - t, x + t);
        let mut i = self.mz.partition_point(|&m| m < lo);
        let mut best: Option<usize> = None;
        while i < self.mz.len() && self.mz[i] <= hi {
            if best.map_or(true, |b| self.intensity[i] >= self.intensity[b]) {
                best = Some(i);
            }
            i += 1;
        }
        best
    }

    /// [`lookup`](Self::lookup) for a caller whose queries mostly move forward: `cursor` is the
    /// previous query's window start and is walked (either way) to this window's start instead of
    /// binary-searched. Same window, same tie rule, same result.
    #[inline]
    fn lookup_walk(&self, x: f32, cursor: &mut usize) -> Option<usize> {
        let t = self.model.tolerance_at(x);
        let (lo, hi) = (x - t, x + t);
        let n = self.mz.len();
        let mut i = (*cursor).min(n);
        while i < n && self.mz[i] < lo {
            i += 1;
        }
        while i > 0 && !(self.mz[i - 1] < lo) {
            i -= 1;
        }
        *cursor = i;
        let mut best: Option<usize> = None;
        while i < n && self.mz[i] <= hi {
            if best.map_or(true, |b| self.intensity[i] >= self.intensity[b]) {
                best = Some(i);
            }
            i += 1;
        }
        best
    }

    /// `compute_site(k, prefix)` for every k in `[0, top]`, swept ion by ion: each ion's
    /// theoretical m/z rises with k, so its peak lookups walk the m/z list once instead of
    /// binary-searching it per node. Every k still receives its terms in (segment, ion) order
    /// starting from 0.0, so the f32 sums are those of `compute_site`.
    fn site_table(&self, top: i32, prefix: bool) -> Vec<f32> {
        let mut out = vec![0.0f32; top as usize + 1];
        let model = self.model;
        let segs = model.segments;
        let rmax = model.max_rank as u32;
        for (s, &pi) in self.seg_partition.iter().enumerate() {
            let part = &model.partitions[pi];
            for (t, ion) in part.ions.iter().enumerate() {
                if ion.prefix != prefix {
                    continue;
                }
                let ll = &part.rank_ll[t];
                let mut cursor = 0usize;
                for k in 1..=top {
                    let real = nominal_to_real(k);
                    let x = real / ion.charge as f32 + ion.offset;
                    let seg = (((x / self.parent_mass) * segs as f32) as i32).min(segs - 1);
                    if seg != s as i32 {
                        continue;
                    }
                    let bin = match self.lookup_walk(x, &mut cursor) {
                        Some(p) => (self.rank[p].min(rmax) - 1) as usize,
                        None => rmax as usize,
                    };
                    out[k as usize] += ll[bin];
                }
            }
        }
        out
    }

    /// `compute_alpha(k)` for every k in `[0, top]`, with one forward-walking lookup cursor.
    fn alpha_table(&self, top: i32) -> Vec<f32> {
        let mut out = vec![0.0f32; top as usize + 1];
        let Some(a) = self.anchor else {
            for v in out.iter_mut().skip(1) {
                *v = -1.0;
            }
            return out;
        };
        let cf = a.charge as f32;
        let mut cursor = 0usize;
        for k in 1..=top {
            let x = nominal_to_real(k) / cf + a.offset;
            out[k as usize] = match self.lookup_walk(x, &mut cursor) {
                Some(p) => (self.mz[p] - a.offset) * cf,
                None => -1.0,
            };
        }
        out
    }

    /// Site score of nominal node `k` for one polarity (prefix = true).
    fn compute_site(&self, k: i32, prefix: bool) -> f32 {
        if k == 0 {
            return 0.0;
        }
        let model = self.model;
        let segs = model.segments;
        let rmax = model.max_rank as u32;
        let real = nominal_to_real(k);
        let mut sum = 0.0f32;
        for (s, &pi) in self.seg_partition.iter().enumerate() {
            let part = &model.partitions[pi];
            for (t, ion) in part.ions.iter().enumerate() {
                if ion.prefix != prefix {
                    continue;
                }
                let x = real / ion.charge as f32 + ion.offset;
                let seg = (((x / self.parent_mass) * segs as f32) as i32).min(segs - 1);
                if seg != s as i32 {
                    continue;
                }
                let bin = match self.lookup(x) {
                    Some(p) => (self.rank[p].min(rmax) - 1) as usize,
                    None => rmax as usize,
                };
                sum += part.rank_ll[t][bin];
            }
        }
        sum
    }

    fn compute_alpha(&self, k: i32) -> f32 {
        if k == 0 {
            return 0.0;
        }
        let Some(a) = self.anchor else { return -1.0 };
        let cf = a.charge as f32;
        let x = nominal_to_real(k) / cf + a.offset;
        match self.lookup(x) {
            Some(p) => (self.mz[p] - a.offset) * cf,
            None => -1.0,
        }
    }

    /// Prefix site score Pre(k).
    #[inline]
    pub fn pre(&self, k: i32) -> f32 {
        match self.pre.get(k as usize) {
            Some(&v) if k >= 0 => v,
            _ => self.compute_site(k, true),
        }
    }
    /// Suffix site score Suf(k).
    #[inline]
    pub fn suf(&self, k: i32) -> f32 {
        match self.suf.get(k as usize) {
            Some(&v) if k >= 0 => v,
            _ => self.compute_site(k, false),
        }
    }
    /// Anchor mass alpha(k); negative = not anchored.
    #[inline]
    pub fn alpha(&self, k: i32) -> f32 {
        match self.alpha.get(k as usize) {
            Some(&v) if k >= 0 => v,
            _ => self.compute_alpha(k),
        }
    }

    /// Combined vertex score of a cleavage at prefix mass `a` of a peptide of nominal mass `p`.
    #[inline]
    pub fn vertex(&self, p: i32, a: i32) -> i32 {
        round_score(self.pre(a) + self.suf(p - a))
    }

    /// Edge score between anchor masses `cur` = alpha(k) and `prev` = alpha(k'), residue mass `ma`.
    #[inline]
    pub fn edge_between(&self, cur: f32, prev: f32, ma: f32) -> i32 {
        let i = (cur >= 0.0) as usize + 2 * (prev >= 0.0) as usize;
        if i == 3 {
            let e = self.model.error_scale;
            let d = (cur - prev) - ma;
            let j = round_score(d * e as f32).clamp(-e, e) + e;
            self.edge_both[j as usize]
        } else {
            self.edge_const[i]
        }
    }

    /// EdgeScore(k, k', m_a).
    #[inline]
    pub fn edge(&self, k: i32, kp: i32, ma: f32) -> i32 {
        self.edge_between(self.alpha(k), self.alpha(kp), ma)
    }
}
