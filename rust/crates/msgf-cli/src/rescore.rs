//! `msgf rescore` — recompute MS-GF+ **RawScore**, **DeNovoScore** and **SpecEValue** for a list
//! of peptide-spectrum matches.
//!
//! This is *not* a database search: the candidate peptides come from the input PSM list. For a
//! search over a FASTA, see [`crate::search`].
//!
//! The generating function depends only on (spectrum, precursor mass, isotope range, amino-acid
//! alphabet) — not on any one peptide — so it is built **once per `(scan, charge)`** and shared by
//! every PSM against that spectrum, each of which is then a cheap RawScore + tail lookup.
//!
//! Because the whole PSM list is known up front, the driver runs in two passes per spectrum: the
//! RawScore of every PSM sharing a `(scan, charge)` is computed first (it needs only the
//! [`PreparedSpectrum`]), and the **minimum** of those RawScores becomes the pruning threshold for
//! the one generating function they share ([`score_distribution`]). The tail is bit-identical to
//! the full DP at and above that threshold, and no PSM in the group is ever queried below it.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use msgf_genfunc::{score_distribution, NullModel};
use msgf_scorer::{
    match_and_terminal_score, prepare_with_cache, Candidate, Cleavage, Ms2, PreparedSpectrum,
    Residue, ScoreModel, RESIDUE_ORDER,
};

pub const USAGE: &str = "\
msgf rescore — recompute MS-GF+ scores for a PSM list

USAGE:
    msgf rescore --spectra <FILE.mgf> --psms <PSMS.tsv> [OPTIONS]

Recompute MS-GF+ RawScore, DeNovoScore and SpecEValue for each input PSM.

REQUIRED:
    -s, --spectra <FILE>   MS/MS spectra, MGF format (must carry SCANS=, CHARGE=, PEPMASS=)
    -i, --psms    <FILE>   PSMs to rescore, TSV: columns `scan`, `peptide`, optional `charge`

OPTIONS:
    -p, --param   <FILE>   Scoring model (.param). Default: the bundled HCD/HighRes/Tryptic
                           model trained from MassIVE-KB (CC0) — pass a file for another
                           activation/instrument/enzyme, e.g. MS-GF+'s own models.
    -o, --out     <FILE>   Output TSV (default: stdout)
        --ti      <LO,HI>  Isotope-error range, like MS-GF+ -ti (default: 0,1)
        --aa-probs <FILE>  Amino-acid background probabilities, TSV `residue<TAB>prob`
                           (default: uniform 0.05 — MS-GF+ de novo). Use a database's
                           composition to reproduce a real search's SpecEValue.
        --ox-m             Add variable oxidation on M (+15.994915) to the graph alphabet
        --db-size <N>      If set, also emit EValue = SpecEValue * N (candidate count)
        --threads <N>      Worker threads (default: all cores). Output is identical for
                           every thread count.
    -h, --help             Print this help

PEPTIDE FORMAT (in the --psms file):
    Bare sequence `PEPTIDEK`, optional enzyme context `K.PEPTIDEK.A`, and inline modification
    deltas `+d`/`-d` on the preceding residue, e.g. `SM+15.995PEP` or `+42.011SAMPLER`.
    Only the 20 standard residues are accepted; unknown residues skip the PSM.

NOTES:
    The default alphabet (20 residues, uniform 0.05) matches MS-GF+ de novo. To reproduce a
    specific MS-GF+ *search* bit-for-bit, pass all three of that search's settings:
      --param     MS-GF+'s own .param for the acquisition (e.g. HCD_HighRes_Tryp.param).
                  The bundled default is a DIFFERENT trained model and will not reproduce
                  MS-GF+ — a model is the scoring function, so its numbers are its own.
      --aa-probs  the searched database's composition (not the uniform default).
      --ox-m      the same variable mods the search used, in the graph alphabet.
    With all three, RawScore/DeNovoScore match exactly and SpecEValue to f64 accumulation noise.
";

// ---- configuration / argument parsing --------------------------------------------------------

pub struct Config {
    spectra: PathBuf,
    param: Option<PathBuf>,
    psms: PathBuf,
    out: Option<PathBuf>,
    ti: (i32, i32),
    aa_probs: Option<PathBuf>,
    ox_m: bool,
    db_size: Option<f64>,
    threads: Option<usize>,
}

impl Config {
    pub fn parse(args: &[String]) -> Result<Config, String> {
        let (mut spectra, mut param, mut psms, mut out, mut aa_probs) =
            (None, None, None, None, None);
        let mut ti = (0, 1);
        let (mut ox_m, mut db_size) = (false, None);
        let mut threads = None;
        let mut it = args.iter();
        while let Some(a) = it.next() {
            let mut want = |name: &str| -> Result<String, String> {
                it.next()
                    .cloned()
                    .ok_or_else(|| format!("`{name}` needs a value"))
            };
            match a.as_str() {
                "-s" | "--spectra" => spectra = Some(PathBuf::from(want("--spectra")?)),
                "-p" | "--param" => param = Some(PathBuf::from(want("--param")?)),
                "-i" | "--psms" => psms = Some(PathBuf::from(want("--psms")?)),
                "-o" | "--out" => out = Some(PathBuf::from(want("--out")?)),
                "--aa-probs" => aa_probs = Some(PathBuf::from(want("--aa-probs")?)),
                "--ox-m" => ox_m = true,
                "--db-size" => {
                    db_size = Some(
                        want("--db-size")?
                            .parse()
                            .map_err(|_| "--db-size must be a number")?,
                    )
                }
                "--threads" => {
                    threads = Some(
                        want("--threads")?
                            .parse()
                            .map_err(|_| "--threads must be a positive integer")?,
                    )
                }
                "--ti" => {
                    let v = want("--ti")?;
                    let (lo, hi) = v.split_once(',').ok_or("--ti must be LO,HI (e.g. 0,1)")?;
                    ti = (
                        lo.trim()
                            .parse()
                            .map_err(|_| "--ti LO must be an integer")?,
                        hi.trim()
                            .parse()
                            .map_err(|_| "--ti HI must be an integer")?,
                    );
                }
                "-h" | "--help" => {
                    print!("{USAGE}");
                    std::process::exit(0);
                }
                other => return Err(format!("unexpected argument `{other}`")),
            }
        }
        if ti.0 > ti.1 {
            return Err("--ti LO must be <= HI".into());
        }
        Ok(Config {
            spectra: spectra.ok_or("missing --spectra")?,
            param,
            psms: psms.ok_or("missing --psms")?,
            out,
            ti,
            aa_probs,
            ox_m,
            db_size,
            threads,
        })
    }
}

// ---- the rescoring driver --------------------------------------------------------------------

/// Raw spectrum data keyed by scan, indexed once from the MGF.
struct RawSpectrum {
    charge: Option<i32>,
    precursor_mz: f64,
    peaks: Vec<(f64, f64)>,
}

/// One PSM to rescore.
struct Psm {
    scan: String,
    peptide: String,
    charge: Option<i32>,
}

/// Why an input PSM produced no row. Recorded per PSM and replayed in input order so the driver's
/// stderr is identical to the one-pass-per-PSM version's.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Skip {
    NoSpectrum,
    NoCharge,
    NoGenFunc,
    BadPeptide,
}

/// What one input PSM turned into.
#[derive(Clone, Copy, Debug)]
enum Outcome {
    Row {
        charge: i32,
        raw: i32,
        denovo: i32,
        spec: f64,
    },
    Skip(Skip),
}

pub fn run(cfg: &Config) -> Result<(), String> {
    if let Some(n) = cfg.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .map_err(|e| format!("configuring {n} threads: {e}"))?;
    }

    let (model, model_source) = crate::model::load(cfg.param.as_deref())?;
    crate::model::announce(&model_source, &model);
    let spectra = index_spectra(&cfg.spectra)?;
    let psms = read_psms(&cfg.psms)?;
    let null = build_null(cfg.aa_probs.as_deref(), cfg.ox_m, cfg.ti)?;

    // Open the output *before* scoring. Grouping by `(scan, charge)` means rows can only be emitted
    // once every group is done, but an unwritable `--out` must still fail in the first second rather
    // than after a full multi-minute generating-function run.
    let mut writer: Box<dyn Write> = match &cfg.out {
        Some(p) => Box::new(BufWriter::new(
            File::create(p).map_err(|e| format!("creating {}: {e}", p.display()))?,
        )),
        None => Box::new(BufWriter::new(io::stdout())),
    };

    let outcomes = score_all(&model, &spectra, &psms, &null, /* pruned = */ true);

    let mut header = String::from("scan\tpeptide\tcharge\traw_score\tdenovo_score\tspec_evalue");
    if cfg.db_size.is_some() {
        header.push_str("\tevalue");
    }
    writeln!(writer, "{header}").map_err(io_err)?;

    let (mut scored_n, mut skipped_n) = (0usize, 0usize);
    for (psm, outcome) in psms.iter().zip(&outcomes) {
        match *outcome {
            Outcome::Skip(kind) => {
                let why = match kind {
                    Skip::NoSpectrum => "not in spectra file",
                    Skip::NoCharge => "no charge",
                    Skip::NoGenFunc => "could not build generating function",
                    Skip::BadPeptide => "unparseable peptide",
                };
                eprintln!("skip scan {} ({}): {why}", psm.scan, psm.peptide);
                skipped_n += 1;
            }
            Outcome::Row {
                charge,
                raw,
                denovo,
                spec,
            } => {
                // An exact-zero SpecEValue (RawScore above the support) prints as `-0.000000e0`,
                // the established output of this tool; keep it.
                let spec = if spec == 0.0 { -0.0 } else { spec };
                write!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{:.6e}",
                    psm.scan, psm.peptide, charge, raw, denovo, spec
                )
                .map_err(io_err)?;
                if let Some(n) = cfg.db_size {
                    write!(writer, "\t{:.6e}", spec * n).map_err(io_err)?;
                }
                writeln!(writer).map_err(io_err)?;
                scored_n += 1;
            }
        }
    }
    writer.flush().map_err(io_err)?;
    eprintln!("rescored {scored_n} PSM(s); skipped {skipped_n}");
    Ok(())
}

/// Score every PSM, grouped by `(scan, charge)`.
///
/// The PSM list is read whole, so the input order is only an *output* constraint: the driver walks
/// spectra instead, which is what makes tail pruning available. For each `(scan, charge)` group it
///
/// 1. builds the [`PreparedSpectrum`] once and takes the RawScore of every PSM in the group, then
/// 2. builds the group's single generating function pruned to the **minimum** of those RawScores —
///    the lowest score any of them will ever query — and reads each PSM's tail off it.
///
/// Groups are scored in parallel on the current rayon pool (`--threads`); each outcome is stored at
/// its PSM's index, so the result is identical for every thread count.
///
/// Grouping (rather than two passes over the whole list) is what keeps the memory profile honest:
/// one prepared spectrum and one distribution are live per worker thread, where the previous
/// PSM-ordered driver kept one of each for **every** distinct `(scan, charge)` alive in a cache
/// until the run ended. What it costs is one `Vec<usize>` of PSM indices per group plus a 24-byte
/// [`Outcome`] per PSM, held so rows can be emitted in input order.
///
/// `pruned = false` builds the unpruned distribution instead — kept callable so
/// `pruned_matches_unpruned_bitwise` can assert the two agree to the last bit.
fn score_all(
    model: &ScoreModel,
    spectra: &HashMap<String, RawSpectrum>,
    psms: &[Psm],
    null: &NullModel,
    pruned: bool,
) -> Vec<Outcome> {
    // Resolve spectrum + charge per PSM (in input order, so the skip reasons match), and bucket the
    // survivors by key. `groups` maps a key to its slot in `keyed`, which preserves first-appearance
    // order so the run is deterministic.
    let mut outcomes: Vec<Outcome> = Vec::with_capacity(psms.len());
    let mut groups: HashMap<(&str, i32), usize> = HashMap::new();
    let mut keyed: Vec<((&str, i32), Vec<usize>)> = Vec::new();
    for (i, psm) in psms.iter().enumerate() {
        let Some(raw) = spectra.get(&psm.scan) else {
            outcomes.push(Outcome::Skip(Skip::NoSpectrum));
            continue;
        };
        let charge = match psm.charge.or(raw.charge) {
            Some(c) if c > 0 => c,
            _ => {
                outcomes.push(Outcome::Skip(Skip::NoCharge));
                continue;
            }
        };
        // Provisional: a group whose generating function cannot be built reports exactly this for
        // every one of its PSMs, including ones with unparseable peptides (the previous driver
        // checked the spectrum before the peptide, and that precedence is preserved below).
        outcomes.push(Outcome::Skip(Skip::NoGenFunc));
        let key = (psm.scan.as_str(), charge);
        match groups.get(&key) {
            Some(&slot) => keyed[slot].1.push(i),
            None => {
                groups.insert(key, keyed.len());
                keyed.push((key, vec![i]));
            }
        }
    }

    // Groups are independent (each builds its own prepared spectrum and distribution, and the
    // generating function's scratch is thread-local), so they run in parallel on the current rayon
    // pool; results are scattered back by PSM index, so the output does not depend on the pool.
    use rayon::prelude::*;
    let scored: Vec<Vec<(usize, Outcome)>> = keyed
        .par_iter()
        .with_min_len(4)
        .map_init(
            Vec::new,
            |raws: &mut Vec<i32>, &((scan, charge), ref idxs)| {
                score_group(
                    model,
                    &spectra[scan],
                    charge,
                    idxs,
                    psms,
                    null,
                    pruned,
                    raws,
                )
            },
        )
        .collect();
    for group in scored {
        for (i, o) in group {
            outcomes[i] = o;
        }
    }
    outcomes
}

/// One `(scan, charge)` group of [`score_all`]: the RawScore of every PSM, then the group's single
/// generating function pruned to the lowest of them. Empty if the group has no distribution (its
/// PSMs keep the provisional `Skip::NoGenFunc`).
#[allow(clippy::too_many_arguments)]
fn score_group(
    model: &ScoreModel,
    raw_spectrum: &RawSpectrum,
    charge: i32,
    idxs: &[usize],
    psms: &[Psm],
    null: &NullModel,
    pruned: bool,
    raws: &mut Vec<i32>,
) -> Vec<(usize, Outcome)> {
    let Some(prep) = prepare_spec(model, raw_spectrum, charge, null.isotope) else {
        return Vec::new();
    };
    // Pass 1: RawScores. `i32::MIN` marks an unparseable peptide — no real RawScore can reach
    // it, so it neither lowers the threshold nor becomes a row.
    raws.clear();
    let mut threshold = i32::MAX;
    for &i in idxs {
        match raw_score_of(&prep, &null.cleavage, &psms[i].peptide) {
            Some(r) => {
                threshold = threshold.min(r);
                raws.push(r);
            }
            None => raws.push(i32::MIN),
        }
    }
    // Pass 2: the group's generating function, pruned to the lowest score it will be asked for.
    // A group with no scorable PSM leaves `threshold` at `i32::MAX`; the DP clamps the cut to
    // the DeNovoScore, so that is simply the cheapest exact run, and it is still needed to
    // decide whether these PSMs skip as "unparseable" or as "no generating function".
    let Some(tail) = score_distribution(&prep, null, pruned.then_some(threshold)) else {
        return Vec::new();
    };
    let denovo = tail.best_possible();
    idxs.iter()
        .zip(raws.iter())
        .map(|(&i, &raw)| {
            let o = if raw == i32::MIN {
                Outcome::Skip(Skip::BadPeptide)
            } else {
                Outcome::Row {
                    charge,
                    raw,
                    denovo,
                    spec: tail.tail_mass(raw),
                }
            };
            (i, o)
        })
        .collect()
}

/// RawScore = the node + edge match score plus the enzymatic-terminus terms, so the SpecEValue tail
/// is looked up at the score the generating function also uses. `None` if the peptide does not
/// parse.
fn raw_score_of(prep: &PreparedSpectrum, cleavage: &Cleavage, peptide: &str) -> Option<i32> {
    let parsed = msgf_chem::peptide::parse(peptide)?;
    let residues: Vec<Residue> = parsed
        .iter()
        .map(|r| Residue {
            letter: r.aa,
            delta: r.mod_delta,
        })
        .collect();
    let cand = Candidate {
        residues: &residues,
        n_term_credit: n_term_credit(peptide),
    };
    match_and_terminal_score(prep, &cand, cleavage)
}

/// Prepare one `(scan, charge)`: `None` if the precursor is implausible (nominal peptide mass
/// outside 50..=10000) or the model has no partition for it.
fn prepare_spec<'m>(
    model: &'m ScoreModel,
    raw: &RawSpectrum,
    charge: i32,
    ti: (i32, i32),
) -> Option<PreparedSpectrum<'m>> {
    let ms2 = Ms2 {
        peaks: &raw.peaks,
        precursor_mz: raw.precursor_mz,
        charge,
    };
    // Node tables are cached up to the highest isotope sink (N0 - LO); others are computed on use.
    prepare_with_cache(model, &ms2, (-ti.0).max(2))
}

// ---- input parsing ---------------------------------------------------------------------------

/// Index every MGF spectrum by its `SCANS=` value.
fn index_spectra(path: &Path) -> Result<HashMap<String, RawSpectrum>, String> {
    // Open first so a missing file reports as before ("opening ..."); then read in parallel.
    File::open(path).map_err(|e| format!("opening {}: {e}", path.display()))?;
    let mut out = HashMap::new();
    let spectra =
        msgf_io::read_mgf_file(path).map_err(|e| format!("reading {}: {e}", path.display()))?;
    for s in spectra {
        let (Some(scan), Some(mz)) = (s.scan, s.precursor_mz) else {
            continue; // need a scan id and a precursor to score
        };
        out.insert(
            scan,
            RawSpectrum {
                charge: s.charge,
                precursor_mz: mz,
                peaks: s.peaks.iter().map(|p| (p.mz, p.intensity)).collect(),
            },
        );
    }
    if out.is_empty() {
        return Err(format!("no usable spectra in {}", path.display()));
    }
    Ok(out)
}

/// Read the PSM TSV. If the first line looks like a header (contains "peptide"), column order is
/// taken from it (case-insensitive `scan`/`peptide`/`charge`); otherwise columns are assumed to be
/// `scan`, `peptide`, optional `charge` in that order.
fn read_psms(path: &Path) -> Result<Vec<Psm>, String> {
    let file = File::open(path).map_err(|e| format!("opening {}: {e}", path.display()))?;
    let mut lines = BufReader::new(file).lines();

    let first = loop {
        match lines.next() {
            Some(l) => {
                let l = l.map_err(io_err)?;
                if !l.trim().is_empty() {
                    break l;
                }
            }
            None => return Err(format!("{} is empty", path.display())),
        }
    };

    // Locate columns.
    let lower = first.to_ascii_lowercase();
    let (scan_i, pep_i, charge_i, mut pending) = if lower.contains("peptide") {
        let cols: Vec<&str> = first.split('\t').map(str::trim).collect();
        let find = |name: &str| cols.iter().position(|c| c.eq_ignore_ascii_case(name));
        (
            find("scan").ok_or("PSM header has no `scan` column")?,
            find("peptide").ok_or("PSM header has no `peptide` column")?,
            find("charge"),
            None,
        )
    } else {
        // No header: treat the first line as data with fixed column order.
        (0, 1, Some(2), Some(first))
    };

    let mut out = Vec::new();
    let push_row = |line: &str, out: &mut Vec<Psm>| -> Result<(), String> {
        if line.trim().is_empty() {
            return Ok(());
        }
        let f: Vec<&str> = line.split('\t').map(str::trim).collect();
        let scan = f.get(scan_i).filter(|s| !s.is_empty());
        let peptide = f.get(pep_i).filter(|s| !s.is_empty());
        let (Some(scan), Some(peptide)) = (scan, peptide) else {
            return Ok(()); // skip malformed row
        };
        let charge = charge_i
            .and_then(|i| f.get(i))
            .and_then(|c| c.trim_end_matches(['+', '-']).parse::<i32>().ok());
        out.push(Psm {
            scan: scan.to_string(),
            peptide: peptide.to_string(),
            charge,
        });
        Ok(())
    };

    if let Some(first_data) = pending.take() {
        push_row(&first_data, &mut out)?;
    }
    for l in lines {
        push_row(&l.map_err(io_err)?, &mut out)?;
    }
    if out.is_empty() {
        return Err(format!("no PSMs parsed from {}", path.display()));
    }
    Ok(out)
}

/// Build the null model: the 20 standard residues with their background probabilities (uniform
/// 0.05 by default, or `--aa-probs`), oxidised methionine (+15.994915) appended with `--ox-m`,
/// trypsin terminus scoring (mixture weight = w_K + w_R) and the isotope-error range.
fn build_null(aa_probs: Option<&Path>, ox_m: bool, ti: (i32, i32)) -> Result<NullModel, String> {
    let probs: HashMap<u8, f64> = match aa_probs {
        Some(p) => load_aa_probs(p)?,
        None => RESIDUE_ORDER.iter().map(|&r| (r, 0.05)).collect(),
    };
    let prob_of = |r: u8| probs.get(&r).copied().unwrap_or(0.05);
    let extra: &[(u8, f64)] = if ox_m { &[(b'M', 15.994915)] } else { &[] };
    Ok(NullModel::from_probs(
        prob_of,
        extra,
        Cleavage::trypsin(),
        ti,
    ))
}

/// Load a residue→probability TSV (`R<TAB>0.0567`), one residue per line, `#` comments allowed.
fn load_aa_probs(path: &Path) -> Result<HashMap<u8, f64>, String> {
    let file = File::open(path).map_err(|e| format!("opening {}: {e}", path.display()))?;
    let mut out = HashMap::new();
    for (n, line) in BufReader::new(file).lines().enumerate() {
        let line = line.map_err(io_err)?;
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let mut f = t.split_whitespace();
        let (Some(res), Some(prob)) = (f.next(), f.next()) else {
            return Err(format!(
                "{}:{}: expected `residue<TAB>prob`",
                path.display(),
                n + 1
            ));
        };
        let r = res.as_bytes()[0].to_ascii_uppercase();
        let p: f64 = prob
            .parse()
            .map_err(|_| format!("{}:{}: bad probability `{prob}`", path.display(), n + 1))?;
        out.insert(r, p);
    }
    if out.is_empty() {
        return Err(format!("no probabilities in {}", path.display()));
    }
    Ok(out)
}

/// N-terminal (neighbouring) cleavage credit for trypsin, from the peptide string's context:
/// credit if the flanking N residue is K/R or a protein terminus (`-`); a **bare** peptide (no
/// `X.….Y` context) is assumed fully tryptic. The C-terminal term (credit iff the last residue is
/// K or R; ending the protein does not substitute) is applied by the scorer.
fn n_term_credit(pep: &str) -> bool {
    let b = pep.as_bytes();
    let has_ctx = b.len() >= 4 && b[1] == b'.' && b[b.len() - 2] == b'.';
    if has_ctx {
        matches!(b[0], b'K' | b'R' | b'-')
    } else {
        true
    }
}

pub fn io_err(e: io::Error) -> String {
    format!("I/O error: {e}")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn repo(rel: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../..")
            .join(rel)
    }

    /// The per-group tail-pruned driver must reproduce the unpruned one **bit for bit**: same
    /// outcome for every PSM, identical integer RawScore/DeNovoScore, and a SpecEValue equal as an
    /// `f64` bit pattern (not within an epsilon — the prune is exact, so anything less tests
    /// nothing). Skipped when the reference spectra/model/golden are absent.
    #[test]
    fn pruned_matches_unpruned_bitwise() {
        let mgf = repo("validation/data/spectra/F13.mgf");
        let param = repo("validation/data/models/HCD_HighRes_Tryp.param");
        let list = repo("validation/golden/iprg2013_F13.tsv");
        if !mgf.exists() || !param.exists() || !list.exists() {
            eprintln!("skip: F13 spectra / HighRes model / golden PSM list absent");
            return;
        }
        let (model, _) = crate::model::load(Some(param.as_path())).expect("model");
        let spectra = index_spectra(&mgf).expect("spectra");
        let null = build_null(None, true, (0, 1)).expect("null model");

        // MS-GF+'s own F13 output: ScanNum, Charge, Peptide. Several PSMs share a scan, which is
        // the case the group minimum exists for. Two deliberately broken rows exercise the skip
        // paths (unparseable peptide, and a scan that is not in the MGF).
        let text = std::fs::read_to_string(&list).expect("golden PSM list");
        let mut psms: Vec<Psm> = text
            .lines()
            .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
            .filter_map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                Some(Psm {
                    scan: f.get(2)?.to_string(),
                    peptide: f.get(9)?.to_string(),
                    charge: f.get(8).and_then(|c| c.parse().ok()),
                })
            })
            .take(250) // enough distinct groups and thresholds; keeps a debug `cargo test` brisk
            .collect();
        assert!(psms.len() > 100, "expected a populated golden PSM list");
        let borrowed_scan = psms[0].scan.clone();
        psms.push(Psm {
            scan: borrowed_scan,
            peptide: "PEPTIDEJ".into(), // J is not a standard residue
            charge: Some(2),
        });
        psms.push(Psm {
            scan: "no-such-scan".into(),
            peptide: "PEPTIDEK".into(),
            charge: Some(2),
        });

        let pruned = score_all(&model, &spectra, &psms, &null, true);
        let full = score_all(&model, &spectra, &psms, &null, false);

        assert_eq!(pruned.len(), psms.len());
        assert_eq!(full.len(), psms.len());
        let mut rows = 0usize;
        for (i, (p, f)) in pruned.iter().zip(&full).enumerate() {
            match (*p, *f) {
                (Outcome::Skip(a), Outcome::Skip(b)) => assert_eq!(a, b, "PSM {i} skip reason"),
                (
                    Outcome::Row {
                        charge: ca,
                        raw: ra,
                        denovo: da,
                        spec: sa,
                    },
                    Outcome::Row {
                        charge: cb,
                        raw: rb,
                        denovo: db,
                        spec: sb,
                    },
                ) => {
                    assert_eq!((ca, ra, da), (cb, rb, db), "PSM {i} integer scores");
                    assert_eq!(
                        sa.to_bits(),
                        sb.to_bits(),
                        "PSM {i}: SpecEValue {sa:e} vs {sb:e} differ in bits"
                    );
                    rows += 1;
                }
                (a, b) => panic!("PSM {i}: outcome kind differs: {a:?} vs {b:?}"),
            }
        }
        assert!(rows > 100, "expected scored rows, got {rows}");
        // The two injected rows must have skipped, for the reasons the old driver gave.
        assert!(matches!(
            pruned[psms.len() - 2],
            Outcome::Skip(Skip::BadPeptide)
        ));
        assert!(matches!(
            pruned[psms.len() - 1],
            Outcome::Skip(Skip::NoSpectrum)
        ));
    }

    /// `score_all` on a multi-thread pool must equal the single-thread run bit for bit, PSM by
    /// PSM (rows and skip reasons, in input order). Synthetic spectra and the bundled model, so it
    /// needs no external data.
    #[test]
    fn thread_count_does_not_change_outcomes() {
        let (model, _) = crate::model::load(None).expect("bundled model");
        let null = build_null(None, true, (-1, 2)).expect("null model");
        let peptides = [
            "SAMPLERK",
            "PEPTIDEK",
            "LGEHNIDVLEGNEQFINAAK",
            "VGAHAGEYGAEALER",
            "MSFVTTR",
            "AEFAEVSK",
            "YLYEIAR",
            "DGNASGTTLLEALDCILPPTRPTDKPLR",
            "HM+15.995VLAGR",
            "QNCELFEQLGEYK",
        ];
        let mut seed = 12345u64;
        let mut rnd = move || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut spectra = HashMap::new();
        let mut psms = Vec::new();
        for (k, pep) in peptides.iter().enumerate() {
            let plain: String = pep.chars().filter(|c| c.is_ascii_uppercase()).collect();
            let neutral = msgf_chem::peptide_neutral_mass(&plain).unwrap();
            let charge = 2 + (k % 2) as i32;
            // b/y-like peaks of the peptide plus noise, so some RawScores are high.
            let mut peaks: Vec<(f64, f64)> = Vec::new();
            let mut acc = 0.0;
            for c in plain.bytes() {
                acc += msgf_chem::residue_mass(c).unwrap();
                peaks.push((acc + msgf_chem::mass::PROTON, 1e5 * (0.5 + rnd())));
                peaks.push((neutral - acc + msgf_chem::mass::PROTON, 1e5 * (0.5 + rnd())));
            }
            for _ in 0..80 {
                peaks.push((100.0 + 1500.0 * rnd(), 1e4 * rnd()));
            }
            peaks.retain(|p| p.0 > 50.0 && p.0 < neutral);
            peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            let scan = format!("{}", 100 + k);
            spectra.insert(
                scan.clone(),
                RawSpectrum {
                    charge: Some(charge),
                    precursor_mz: (neutral + charge as f64 * msgf_chem::mass::PROTON)
                        / charge as f64,
                    peaks,
                },
            );
            // Each spectrum: its own peptide, two others, and an explicit second charge group.
            for (j, q) in [k, (k + 3) % peptides.len(), (k + 7) % peptides.len()]
                .iter()
                .enumerate()
            {
                psms.push(Psm {
                    scan: scan.clone(),
                    peptide: peptides[*q].to_string(),
                    charge: if j == 2 { Some(charge + 1) } else { None },
                });
            }
        }
        psms.push(Psm {
            scan: "100".into(),
            peptide: "PEPTIDEJ".into(),
            charge: None,
        });
        psms.push(Psm {
            scan: "no-such-scan".into(),
            peptide: "PEPTIDEK".into(),
            charge: Some(2),
        });

        let run = |n: usize| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .unwrap()
                .install(|| score_all(&model, &spectra, &psms, &null, true))
        };
        let key = |o: &Outcome| match *o {
            Outcome::Row {
                charge,
                raw,
                denovo,
                spec,
            } => format!("{charge} {raw} {denovo} {:016x}", spec.to_bits()),
            Outcome::Skip(s) => format!("{s:?}"),
        };
        let one: Vec<String> = run(1).iter().map(key).collect();
        assert_eq!(one.len(), psms.len());
        assert!(
            one.iter()
                .filter(|r| !r.starts_with(char::is_alphabetic))
                .count()
                >= 20
        );
        for n in [2, 3, 8] {
            let many: Vec<String> = run(n).iter().map(key).collect();
            assert_eq!(one, many, "{n} threads differ from 1 thread");
        }
    }
}
