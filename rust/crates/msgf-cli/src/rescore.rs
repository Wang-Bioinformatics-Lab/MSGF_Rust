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
//! [`ScoredSpectrum`]), and the **minimum** of those RawScores becomes the pruning threshold for
//! the one generating function they share ([`compute_tail_into`]). The tail is bit-identical to
//! the full DP at and above that threshold, and no PSM in the group is ever queried below it.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use msgf_chem::{mass, scaling};
use msgf_genfunc::graph::{build_reverse_graph, standard_aa_nominal, Aa, PeptideCleavage};
use msgf_genfunc::{compute, compute_tail_into, merge_group, Cleavage, DpScratch, GenFunc};
use msgf_io::MgfReader;
use msgf_scorer::preprocess::preprocess;
use msgf_scorer::scored_spectrum::ScoredSpectrum;
use msgf_scorer::ScoringModel;

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
}

impl Config {
    pub fn parse(args: &[String]) -> Result<Config, String> {
        let (mut spectra, mut param, mut psms, mut out, mut aa_probs) =
            (None, None, None, None, None);
        let mut ti = (0, 1);
        let (mut ox_m, mut db_size) = (false, None);
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
        })
    }
}

// ---- the rescoring driver --------------------------------------------------------------------

/// Raw spectrum data keyed by scan, indexed once from the MGF.
struct RawSpectrum {
    charge: Option<i32>,
    precursor_mz: f64,
    peaks: Vec<(f32, f32)>,
}

/// One PSM to rescore.
struct Psm {
    scan: String,
    peptide: String,
    charge: Option<i32>,
}

/// The candidate-independent half of a `(scan, charge)`: the scored spectrum (all a RawScore
/// needs) plus the isotope-error sink masses the generating function will be built over.
struct PreparedSpec<'m> {
    scored: ScoredSpectrum<'m>,
    sinks: Vec<i32>,
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
    let (model, model_source) = crate::model::load(cfg.param.as_deref())?;
    crate::model::announce(&model_source, &model);
    let spectra = index_spectra(&cfg.spectra)?;
    let psms = read_psms(&cfg.psms)?;
    let (aa, prob_cleavage) = build_alphabet(cfg.aa_probs.as_deref(), cfg.ox_m)?;

    // Open the output *before* scoring. Grouping by `(scan, charge)` means rows can only be emitted
    // once every group is done, but an unwritable `--out` must still fail in the first second rather
    // than after a full multi-minute generating-function run.
    let mut writer: Box<dyn Write> = match &cfg.out {
        Some(p) => Box::new(BufWriter::new(
            File::create(p).map_err(|e| format!("creating {}: {e}", p.display()))?,
        )),
        None => Box::new(BufWriter::new(io::stdout())),
    };

    let outcomes = score_all(
        &model,
        &spectra,
        &psms,
        &aa,
        prob_cleavage,
        cfg.ti,
        /* pruned = */ true,
    );

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
/// 1. builds the [`ScoredSpectrum`] once and takes the RawScore of every PSM in the group, then
/// 2. builds the group's single generating function pruned to the **minimum** of those RawScores —
///    the lowest score any of them will ever query — and reads each PSM's tail off it.
///
/// Grouping (rather than two passes over the whole list) is what keeps the memory profile honest:
/// exactly one `ScoredSpectrum` and one `GenFunc` are live at a time, where the previous
/// PSM-ordered driver kept one of each for **every** distinct `(scan, charge)` alive in a cache
/// until the run ended. What it costs is one `Vec<usize>` of PSM indices per group plus a 24-byte
/// [`Outcome`] per PSM, held so rows can be emitted in input order.
///
/// `pruned = false` runs the unpruned [`compute`] instead — the pre-pruning path, kept callable so
/// `pruned_matches_unpruned_bitwise` can assert the two agree to the last bit.
fn score_all(
    model: &ScoringModel,
    spectra: &HashMap<String, RawSpectrum>,
    psms: &[Psm],
    aa: &[Aa],
    prob_cleavage: f64,
    ti: (i32, i32),
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

    let cleave = Cleavage {
        credit: 2,
        penalty: -11,
        prob_cleavage_sites: prob_cleavage,
    };
    let mut scratch = DpScratch::default();
    let mut raws: Vec<i32> = Vec::new();
    for &((scan, charge), ref idxs) in &keyed {
        let raw_spectrum = &spectra[scan];
        let Some(prep) = prepare_spec(model, raw_spectrum, charge, ti) else {
            continue; // outcomes already carry Skip::NoGenFunc
        };

        // Pass 1: RawScores. `i32::MIN` marks an unparseable peptide — no real RawScore can reach
        // it, so it neither lowers the threshold nor becomes a row.
        raws.clear();
        let mut threshold = i32::MAX;
        for &i in idxs {
            match raw_score_of(&prep.scored, &psms[i].peptide) {
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
        let Some(gf) = build_gf(&prep, aa, cleave, &mut scratch, pruned.then_some(threshold))
        else {
            continue;
        };
        let denovo = gf.max_score();
        for (&i, &raw) in idxs.iter().zip(&raws) {
            outcomes[i] = if raw == i32::MIN {
                Outcome::Skip(Skip::BadPeptide)
            } else {
                Outcome::Row {
                    charge,
                    raw,
                    denovo,
                    spec: gf.spectral_probability(raw),
                }
            };
        }
    }
    outcomes
}

/// MS-GF+ RawScore = node+edge match score (`DBScanScorer.getScore`) + terminal cleavage.
/// `scored.raw_score` is the node+edge part; the peptide/neighboring cleavage the graph scores at
/// the termini is added so the SpecEValue tail is looked up at the same score MS-GF+ reports.
/// `None` if the peptide does not parse.
fn raw_score_of(scored: &ScoredSpectrum, peptide: &str) -> Option<i32> {
    let residues = msgf_chem::peptide::parse(peptide)?;
    let nominal = msgf_chem::peptide::nominal_prefix_masses(&residues);
    let accurate = msgf_chem::peptide::accurate_prefix_masses(&residues);
    let num_mods = msgf_chem::peptide::num_mods(&residues) as i32;
    Some(scored.raw_score(&nominal, &accurate, num_mods) + cleavage_score(peptide, &residues))
}

/// Build the scored spectrum and isotope-error sink set for one `(scan, charge)`. `None` if the
/// precursor is implausible or the sink range is empty.
fn prepare_spec<'m>(
    model: &'m ScoringModel,
    raw: &RawSpectrum,
    charge: i32,
    ti: (i32, i32),
) -> Option<PreparedSpec<'m>> {
    // Neutral precursor mass, then the candidate peptide's nominal mass (precursor − water).
    let parent_mass = raw.precursor_mz as f32 * charge as f32 - charge as f32 * mass::PROTON as f32;
    let pep_nominal = scaling::nominal_bin(parent_mass - mass::WATER as f32);
    if !(50..=10_000).contains(&pep_nominal) {
        return None;
    }
    // Isotope-error sink range (MS-GF+ -ti LO,HI): an isotope error of +k means the measured mass
    // is ~k Da high, so the true peptide mass is k bins lower.
    let sinks: Vec<i32> = (pep_nominal - ti.1..=pep_nominal - ti.0)
        .filter(|&p| p > 0)
        .collect();
    if sinks.is_empty() {
        return None;
    }

    let peaks = preprocess(model, charge, parent_mass, &raw.peaks);
    let scored = ScoredSpectrum::from_ranked_peaks(model, charge, parent_mass, peaks);
    Some(PreparedSpec { scored, sinks })
}

/// Build the merged generating function for one prepared spectrum. With `threshold = Some(t)` the
/// DP discards every score cell that provably cannot reach `t`; the resulting tail is bit-identical
/// to the unpruned one for every score `>= t`, and `max_score()` (DeNovoScore) is exact regardless.
/// `None` means the sinks are unreachable — the same condition the unpruned path returns `None` on,
/// because the cut is clamped to the DeNovoScore and so never empties a reachable graph.
fn build_gf(
    prep: &PreparedSpec,
    aa: &[Aa],
    cleave: Cleavage,
    scratch: &mut DpScratch,
    threshold: Option<i32>,
) -> Option<GenFunc> {
    // GeneratingFunctionGroup: one graph per candidate peptide mass (isotope range), then merged.
    // Tables and edges are candidate-independent, so build them once for the largest candidate and
    // only recompute node scores per candidate.
    let max_p = *prep.sinks.iter().max().unwrap(); // sinks is non-empty (prepare_spec checked)
    let tables = prep.scored.tables(max_p);
    let (mut graph, _) = build_reverse_graph(
        &prep.scored,
        &tables,
        max_p,
        &[max_p],
        aa,
        PeptideCleavage::TRYPSIN,
    );
    let mut gfs: Vec<GenFunc> = Vec::new();
    for &p in &prep.sinks {
        graph.recompute_node_scores(&tables, p, &[p]);
        let gf = match threshold {
            Some(t) => compute_tail_into(scratch, &graph, &[p as usize], Some(cleave), t),
            None => compute(&graph, &[p as usize], Some(cleave)),
        };
        if let Some(gf) = gf {
            gfs.push(gf);
        }
    }
    merge_group(&gfs)
}

// ---- input parsing ---------------------------------------------------------------------------

/// Index every MGF spectrum by its `SCANS=` value.
fn index_spectra(path: &Path) -> Result<HashMap<String, RawSpectrum>, String> {
    let file = File::open(path).map_err(|e| format!("opening {}: {e}", path.display()))?;
    let mut out = HashMap::new();
    for s in MgfReader::new(BufReader::new(file)) {
        let s = s.map_err(|e| format!("reading {}: {e}", path.display()))?;
        let (Some(scan), Some(mz)) = (s.scan.clone(), s.precursor_mz) else {
            continue; // need a scan id and a precursor to score
        };
        out.insert(
            scan,
            RawSpectrum {
                charge: s.charge,
                precursor_mz: mz,
                peaks: s
                    .peaks
                    .iter()
                    .map(|p| (p.mz as f32, p.intensity as f32))
                    .collect(),
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

/// Build the graph amino-acid alphabet and the K+R cleavage probability from a residue→probability
/// map (uniform 0.05 by default). With `--ox-m`, appends oxidized methionine (+15.994915).
fn build_alphabet(aa_probs: Option<&Path>, ox_m: bool) -> Result<(Vec<Aa>, f64), String> {
    let probs: HashMap<u8, f64> = match aa_probs {
        Some(p) => load_aa_probs(p)?,
        None => standard_aa_nominal()
            .iter()
            .map(|(r, _)| (*r, 0.05))
            .collect(),
    };
    let prob_of = |r: u8| probs.get(&r).copied().unwrap_or(0.05);
    let mut aa: Vec<Aa> = standard_aa_nominal()
        .into_iter()
        .map(|(residue, nominal)| Aa {
            residue,
            nominal,
            accurate_mass: msgf_chem::residue_mass(residue).expect("standard residue") as f32,
            prob: prob_of(residue),
        })
        .collect();
    if ox_m {
        let m_ox = msgf_chem::residue_mass(b'M').unwrap() + 15.994915;
        aa.push(Aa {
            residue: b'M',
            nominal: scaling::nominal_bin(m_ox as f32),
            accurate_mass: m_ox as f32,
            prob: prob_of(b'M'),
        });
    }
    Ok((aa, prob_of(b'K') + prob_of(b'R')))
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

/// Terminal cleavage contribution to the MS-GF+ RawScore (trypsin), using the same credit/penalty
/// (+2 / −11) the generating function applies at the peptide/neighboring termini.
///
/// - **C-terminal (peptide) cleavage:** credit if the last residue is K or R; otherwise penalty.
///   Sitting at the protein C-terminus does *not* earn credit.
/// - **N-terminal (neighboring) cleavage:** credit if the flanking N residue is K/R or a protein
///   terminus (`-`); otherwise penalty. A **bare** peptide (no `X.…​.Y` context) is assumed fully
///   tryptic (credit) — supply flanking context to score semi-tryptic termini correctly.
fn cleavage_score(pep: &str, residues: &[msgf_chem::peptide::Residue]) -> i32 {
    const CREDIT: i32 = 2;
    const PENALTY: i32 = -11;
    let b = pep.as_bytes();
    let has_ctx = b.len() >= 4 && b[1] == b'.' && b[b.len() - 2] == b'.';
    let (flank_n, flank_c) = if has_ctx {
        (Some(b[0]), Some(b[b.len() - 1]))
    } else {
        (None, None)
    };
    let is_kr = |c: u8| c == b'K' || c == b'R';

    let nterm = match flank_n {
        Some(c) if is_kr(c) || c == b'-' => CREDIT,
        Some(_) => PENALTY,
        None => CREDIT, // bare peptide: assume tryptic N-terminus
    };
    // Ending the protein does not substitute for a cleavage residue at the C-terminus (verified
    // against MS-GF+ on F13), so `flank_c` only matters for peptides that are not protein-terminal.
    let _ = flank_c;
    let last = residues.last().map(|r| r.aa).unwrap_or(0);
    let cterm = if is_kr(last) { CREDIT } else { PENALTY };
    nterm + cterm
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
        let (aa, prob_cleavage) = build_alphabet(None, true).expect("alphabet");

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

        let pruned = score_all(&model, &spectra, &psms, &aa, prob_cleavage, (0, 1), true);
        let full = score_all(&model, &spectra, &psms, &aa, prob_cleavage, (0, 1), false);

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
}
