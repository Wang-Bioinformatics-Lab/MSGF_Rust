//! `msgf-train` — count a fragment-scoring model (`.param`) from annotated spectra.

use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Instant;

use msgf_train::corpus::{self, CorpusFilter, CorpusStats};
use msgf_train::{counts, TrainConfig};

const USAGE: &str = "\
msgf-train — train a fragment-scoring model (.param) from annotated spectra

USAGE:
    msgf-train --corpus <FILE.mgf>... --out <MODEL.param> [OPTIONS]

The corpus is one or more MGF files whose spectra carry a peptide annotation in `SEQ=`
(MassIVE-KB peptide-library MGFs are exactly this). Training is a counting pass: same
corpus + same options => byte-identical model.

REQUIRED:
    -c, --corpus <FILE>     annotated MGF (repeatable)
    -o, --out    <FILE>     output .param

IDENTITY:
        --activation <S>    default HCD
        --instrument <S>    default HighRes
        --enzyme     <S>    default Tryp ('none' to disable the tryptic C-term filter)

COUNTING:
        --limit <N>         use at most N PSMs (after filtering)
        --max-rank <N>      rank ceiling (default 150)
        --segments <N>      mass segments per precursor (default 2)
        --min-psms <N>      target PSMs per mass partition (default 400)
        --max-parts <N>     cap on mass partitions per charge (default 30)
        --ion-threshold <F> keep ion types seen at >= F of sites (default 0.15)
        --max-ions <N>      cap on scored ion types per partition (default 6)
        --frag-charge <N>   highest fragment charge considered (default 2)
        --mme <F>           fragment matching tolerance in Da (default 0.5)
        --smoothing <F>     add-lambda on rank bins (default 0.005)
        --rank-smoothing <F> rank-pooling width as a fraction of rank (default 0.1)
        --no-precursor-defaults  do not fall back to chemistry-derived precursor offsets

OUTPUT:
        --report <FILE>     write a JSON training report
    -h, --help
";

struct Args {
    corpus: Vec<PathBuf>,
    out: Option<PathBuf>,
    report: Option<PathBuf>,
    limit: Option<usize>,
    cfg: TrainConfig,
}

fn parse_args(argv: &[String]) -> Result<Args, String> {
    let mut a = Args {
        corpus: Vec::new(),
        out: None,
        report: None,
        limit: None,
        cfg: TrainConfig::high_res_hcd_tryptic(),
    };
    let mut i = 0;
    while i < argv.len() {
        let k = argv[i].as_str();
        let mut val = || -> Result<String, String> {
            i += 1;
            argv.get(i)
                .cloned()
                .ok_or_else(|| format!("{k} needs a value"))
        };
        match k {
            "-c" | "--corpus" => a.corpus.push(PathBuf::from(val()?)),
            "-o" | "--out" => a.out = Some(PathBuf::from(val()?)),
            "--report" => a.report = Some(PathBuf::from(val()?)),
            "--limit" => a.limit = Some(val()?.parse().map_err(|_| "bad --limit")?),
            "--activation" => a.cfg.activation = val()?,
            "--instrument" => a.cfg.instrument = val()?,
            "--enzyme" => {
                let v = val()?;
                a.cfg.enzyme = if v == "none" { None } else { Some(v) };
            }
            "--max-rank" => a.cfg.max_rank = val()?.parse().map_err(|_| "bad --max-rank")?,
            "--segments" => a.cfg.num_segments = val()?.parse().map_err(|_| "bad --segments")?,
            "--min-psms" => {
                a.cfg.min_psms_per_partition = val()?.parse().map_err(|_| "bad --min-psms")?
            }
            "--max-parts" => {
                a.cfg.max_partitions_per_charge = val()?.parse().map_err(|_| "bad --max-parts")?
            }
            "--ion-threshold" => {
                a.cfg.ion_freq_threshold = val()?.parse().map_err(|_| "bad --ion-threshold")?
            }
            "--max-ions" => {
                a.cfg.max_ions_per_partition = val()?.parse().map_err(|_| "bad --max-ions")?
            }
            "--frag-charge" => {
                a.cfg.max_fragment_charge = val()?.parse().map_err(|_| "bad --frag-charge")?
            }
            "--mme" => {
                a.cfg.mme = msgf_chem::Tolerance::da(val()?.parse().map_err(|_| "bad --mme")?)
            }
            "--smoothing" => a.cfg.smoothing = val()?.parse().map_err(|_| "bad --smoothing")?,
            "--no-precursor-defaults" => {
                a.cfg.precursor_defaults = false;
                i -= 1; // no value
            }
            "--rank-smoothing" => {
                a.cfg.rank_smoothing = val()?.parse().map_err(|_| "bad --rank-smoothing")?
            }
            "-h" | "--help" => {
                print!("{USAGE}");
                std::process::exit(0);
            }
            other => return Err(format!("unknown flag {other}")),
        }
        i += 1;
    }
    if a.corpus.is_empty() {
        return Err("no --corpus given".into());
    }
    if a.out.is_none() {
        return Err("no --out given".into());
    }
    Ok(a)
}

fn main() -> ExitCode {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let args = match parse_args(&argv) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("msgf-train: {e}\n\n{USAGE}");
            return ExitCode::FAILURE;
        }
    };

    let filter = CorpusFilter {
        charge_min: args.cfg.charge_min,
        charge_max: args.cfg.charge_max,
        require_tryptic_cterm: args.cfg.enzyme.as_deref() == Some("Tryp"),
        ..CorpusFilter::default()
    };

    let t0 = Instant::now();
    let mut psms = Vec::new();
    let mut stats = CorpusStats::default();
    for p in &args.corpus {
        if let Err(e) = corpus::read_annotated_mgf(p, &filter, &mut psms, &mut stats) {
            eprintln!("msgf-train: {}: {e}", p.display());
            return ExitCode::FAILURE;
        }
        eprintln!(
            "  read {:<44} kept {:>7} / {:>7}",
            p.file_name().unwrap_or_default().to_string_lossy(),
            stats.kept,
            stats.read
        );
    }
    if let Some(n) = args.limit {
        psms.truncate(n);
    }
    let t_read = t0.elapsed();
    if psms.is_empty() {
        eprintln!("msgf-train: corpus is empty after filtering");
        return ExitCode::FAILURE;
    }
    eprintln!(
        "corpus: {} PSMs kept of {} read in {:.1}s  (rejected: charge {}, length {}, C-term {}, \
         residue {}, peaks {}, mass {}, unannotated {})",
        psms.len(),
        stats.read,
        t_read.as_secs_f64(),
        stats.charge,
        stats.length,
        stats.cterm,
        stats.bad_residue,
        stats.peaks,
        stats.mass_mismatch,
        stats.no_annotation
    );

    let t1 = Instant::now();
    let (model, report, scheme, spectra) = counts::train(&psms, &args.cfg);
    let t_train = t1.elapsed();

    let out = args.out.unwrap();
    if let Err(e) = msgf_scorer::write_param_file(&out, &model) {
        eprintln!("msgf-train: writing {}: {e}", out.display());
        return ExitCode::FAILURE;
    }
    let size = std::fs::metadata(&out).map(|m| m.len()).unwrap_or(0);

    let scored: usize = report.iter().filter(|r| !r.ions.is_empty()).count();
    eprintln!(
        "trained {} partitions ({} with scored ions) from {} spectra in {:.1}s -> {} ({} bytes)",
        scheme.partitions.len(),
        scored,
        spectra,
        t_train.as_secs_f64(),
        out.display(),
        size
    );
    for (c, b) in &scheme.boundaries {
        eprintln!("  charge {c}: {} mass partitions", b.len());
    }

    if let Some(rp) = args.report {
        let mut ion_hist: Vec<(String, usize)> = Vec::new();
        for r in &report {
            for (label, _, _) in &r.ions {
                match ion_hist.iter_mut().find(|(l, _)| l == label) {
                    Some(e) => e.1 += 1,
                    None => ion_hist.push((label.clone(), 1)),
                }
            }
        }
        ion_hist.sort_by(|a, b| b.1.cmp(&a.1));
        let mut s = String::from("{\n");
        s.push_str(&format!("  \"psms\": {},\n", psms.len()));
        s.push_str(&format!("  \"spectra_counted\": {spectra},\n"));
        s.push_str(&format!(
            "  \"read_seconds\": {:.2},\n",
            t_read.as_secs_f64()
        ));
        s.push_str(&format!(
            "  \"train_seconds\": {:.2},\n",
            t_train.as_secs_f64()
        ));
        s.push_str(&format!("  \"param_bytes\": {size},\n"));
        s.push_str(&format!("  \"partitions\": {},\n", scheme.partitions.len()));
        s.push_str(&format!("  \"max_rank\": {},\n", model.max_rank));
        s.push_str(&format!(
            "  \"precursor_offsets\": {},\n",
            model.precursor_off.len()
        ));
        s.push_str("  \"ion_types\": {");
        for (i, (l, n)) in ion_hist.iter().enumerate() {
            s.push_str(&format!("{}\"{l}\": {n}", if i == 0 { "" } else { ", " }));
        }
        s.push_str("},\n  \"partition_detail\": [\n");
        for (i, r) in report.iter().enumerate() {
            s.push_str(&format!(
                "    {{\"charge\": {}, \"parent_mass\": {:.4}, \"seg\": {}, \"spectra\": {}, \
                 \"main_ion\": \"{}\", \"ions\": [{}]}}{}\n",
                r.charge,
                r.parent_mass,
                r.seg,
                r.spectra,
                r.main_ion.clone().unwrap_or_default(),
                r.ions
                    .iter()
                    .map(|(l, n, f)| format!("[\"{l}\", \"{n}\", {f:.4}]"))
                    .collect::<Vec<_>>()
                    .join(", "),
                if i + 1 == report.len() { "" } else { "," }
            ));
        }
        s.push_str("  ]\n}\n");
        if let Err(e) = std::fs::write(&rp, s) {
            eprintln!("msgf-train: writing {}: {e}", rp.display());
            return ExitCode::FAILURE;
        }
        eprintln!("report -> {}", rp.display());
    }

    ExitCode::SUCCESS
}
