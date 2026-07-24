//! `msgf rescore` — recompute MS-GF+ **RawScore**, **DeNovoScore** and **SpecEValue** for a list
//! of peptide-spectrum matches.
//!
//! This is *not* a database search: the candidate peptides come from the input PSM list. For a
//! search over a FASTA, see [`crate::search`].
//!
//! The generating function depends only on (spectrum, precursor mass, isotope range, amino-acid
//! alphabet) — not on any one peptide — so it is built **once per `(scan, charge)`** and cached;
//! every PSM against that spectrum is then a cheap RawScore + tail lookup.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use msgf_chem::{mass, scaling};
use msgf_genfunc::graph::{build_reverse_graph, standard_aa_nominal, Aa, PeptideCleavage};
use msgf_genfunc::{compute, merge_group, Cleavage, GenFunc};
use msgf_io::MgfReader;
use msgf_scorer::preprocess::preprocess;
use msgf_scorer::scored_spectrum::ScoredSpectrum;
use msgf_scorer::ScoringModel;

pub const USAGE: &str = "\
msgf rescore — recompute MS-GF+ scores for a PSM list

USAGE:
    msgf rescore --spectra <FILE.mgf> --param <MODEL.param> --psms <PSMS.tsv> [OPTIONS]

Recompute MS-GF+ RawScore, DeNovoScore and SpecEValue for each input PSM.

REQUIRED:
    -s, --spectra <FILE>   MS/MS spectra, MGF format (must carry SCANS=, CHARGE=, PEPMASS=)
    -p, --param   <FILE>   MS-GF+ scoring model (.param, e.g. HCD_HighRes_Tryp.param)
    -i, --psms    <FILE>   PSMs to rescore, TSV: columns `scan`, `peptide`, optional `charge`

OPTIONS:
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
    specific MS-GF+ *search* bit-for-bit, pass that search's --aa-probs (DB composition) and the
    same variable mods (--ox-m for oxidation on M); RawScore/DeNovoScore then match exactly and
    SpecEValue to f64 accumulation noise.
";

// ---- configuration / argument parsing --------------------------------------------------------

pub struct Config {
    spectra: PathBuf,
    param: PathBuf,
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
            param: param.ok_or("missing --param")?,
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

/// Everything a `(scan, charge)` spectrum needs to score any peptide against it: the scored
/// spectrum (for RawScore) and the prebuilt generating function (for DeNovoScore + SpecEValue).
struct Prepared<'m> {
    scored: ScoredSpectrum<'m>,
    gf: GenFunc,
}

pub fn run(cfg: &Config) -> Result<(), String> {
    let model = msgf_scorer::read_param_file(&cfg.param)
        .map_err(|e| format!("reading model {}: {e:?}", cfg.param.display()))?;
    let spectra = index_spectra(&cfg.spectra)?;
    let psms = read_psms(&cfg.psms)?;
    let (aa, prob_cleavage) = build_alphabet(cfg.aa_probs.as_deref(), cfg.ox_m)?;

    let mut writer: Box<dyn Write> = match &cfg.out {
        Some(p) => Box::new(BufWriter::new(
            File::create(p).map_err(|e| format!("creating {}: {e}", p.display()))?,
        )),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let mut header = String::from("scan\tpeptide\tcharge\traw_score\tdenovo_score\tspec_evalue");
    if cfg.db_size.is_some() {
        header.push_str("\tevalue");
    }
    writeln!(writer, "{header}").map_err(io_err)?;

    // GenFunc + ScoredSpectrum cache, one entry per (scan, charge). `None` = tried and unscorable.
    let mut cache: HashMap<(String, i32), Option<Prepared>> = HashMap::new();
    let (mut scored_n, mut skipped_n) = (0usize, 0usize);

    for psm in &psms {
        let raw = match spectra.get(&psm.scan) {
            Some(r) => r,
            None => {
                eprintln!(
                    "skip scan {} ({}): not in spectra file",
                    psm.scan, psm.peptide
                );
                skipped_n += 1;
                continue;
            }
        };
        let charge = match psm.charge.or(raw.charge) {
            Some(c) if c > 0 => c,
            _ => {
                eprintln!("skip scan {} ({}): no charge", psm.scan, psm.peptide);
                skipped_n += 1;
                continue;
            }
        };

        let prepared = cache
            .entry((psm.scan.clone(), charge))
            .or_insert_with(|| prepare(&model, raw, charge, &aa, prob_cleavage, cfg.ti));
        let Some(prepared) = prepared.as_ref() else {
            eprintln!(
                "skip scan {} ({}): could not build generating function",
                psm.scan, psm.peptide
            );
            skipped_n += 1;
            continue;
        };

        let Some(residues) = msgf_chem::peptide::parse(&psm.peptide) else {
            eprintln!(
                "skip scan {} ({}): unparseable peptide",
                psm.scan, psm.peptide
            );
            skipped_n += 1;
            continue;
        };
        let nominal = msgf_chem::peptide::nominal_prefix_masses(&residues);
        let accurate = msgf_chem::peptide::accurate_prefix_masses(&residues);
        let num_mods = msgf_chem::peptide::num_mods(&residues) as i32;

        // MS-GF+ RawScore = node+edge match score (DBScanScorer.getScore) + terminal cleavage.
        // `scored.raw_score` is the node+edge part; add the peptide/neighboring cleavage the graph
        // scores at the termini so the SpecEValue tail is looked up at the same score MS-GF+ reports.
        let raw_score = prepared.scored.raw_score(&nominal, &accurate, num_mods)
            + cleavage_score(&psm.peptide, &residues);
        let denovo = prepared.gf.max_score();
        let spec = prepared.gf.spectral_probability(raw_score);

        write!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{:.6e}",
            psm.scan, psm.peptide, charge, raw_score, denovo, spec
        )
        .map_err(io_err)?;
        if let Some(n) = cfg.db_size {
            write!(writer, "\t{:.6e}", spec * n).map_err(io_err)?;
        }
        writeln!(writer).map_err(io_err)?;
        scored_n += 1;
    }
    writer.flush().map_err(io_err)?;
    eprintln!("rescored {scored_n} PSM(s); skipped {skipped_n}");
    Ok(())
}

/// Build the scored spectrum and generating function for one `(scan, charge)`. `None` if the
/// precursor is implausible or the sinks are unreachable.
fn prepare<'m>(
    model: &'m ScoringModel,
    raw: &RawSpectrum,
    charge: i32,
    aa: &[Aa],
    prob_cleavage: f64,
    ti: (i32, i32),
) -> Option<Prepared<'m>> {
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
    let cleave = Cleavage {
        credit: 2,
        penalty: -11,
        prob_cleavage_sites: prob_cleavage,
    };
    // GeneratingFunctionGroup: one graph per candidate peptide mass (isotope range), then merged.
    // Tables and edges are candidate-independent, so build them once for the largest candidate and
    // only recompute node scores per candidate.
    let max_p = *sinks.iter().max().unwrap(); // sinks is non-empty (checked above)
    let tables = scored.tables(max_p);
    let (mut graph, _) = build_reverse_graph(
        &scored,
        &tables,
        max_p,
        &[max_p],
        aa,
        PeptideCleavage::TRYPSIN,
    );
    let mut gfs: Vec<GenFunc> = Vec::new();
    for &p in &sinks {
        graph.recompute_node_scores(&tables, p, &[p]);
        if let Some(gf) = compute(&graph, &[p as usize], Some(cleave)) {
            gfs.push(gf);
        }
    }
    let gf = merge_group(&gfs)?;
    Some(Prepared { scored, gf })
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
