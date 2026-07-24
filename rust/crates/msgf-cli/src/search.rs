//! `msgf search` — a full database search: FASTA → digest → candidates → RawScore → SpecEValue →
//! target-decoy q-values.
//!
//! The engine lives in `msgf-search`; this module is argument parsing plus reporting.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use msgf_chem::{Tolerance, Unit};
use msgf_db::decoy::normalize_prefix;
use msgf_db::enzyme::{DigestParams, Enzyme};
use msgf_db::fasta::{DecoyStrategy, ProteinDb};
use msgf_search::index::PeptideIndex;
use msgf_search::mods::{ModSet, ModSpec};
use msgf_search::{assign_q_values, report, SearchEngine, SearchParams};

pub const USAGE: &str = "\
msgf search — search MS/MS spectra against a protein database

USAGE:
    msgf search --spectra <FILE.mgf> --param <MODEL.param> --fasta <DB.fasta> [OPTIONS]

REQUIRED:
    -s, --spectra <FILE>   MS/MS spectra, MGF format (SCANS=, CHARGE=, PEPMASS=)
    -p, --param   <FILE>   MS-GF+ scoring model (.param, e.g. HCD_HighRes_Tryp.param)
    -d, --fasta   <FILE>   Protein database, FASTA

OUTPUT:
    -o, --out     <FILE>   Results TSV (default: stdout), MS-GF+ column set
        --unroll           One row per protein occurrence (MS-GF+ `-unroll 1`)

SEARCH SPACE:
    -e, --enzyme  <SPEC>   Enzyme: a number 0-10, a name (Tryp, LysC, ...), or a full
                           `Name,CleaveAt,Terminus,Description` definition   [default: 1 = Tryp]
        --ntt     <N>      Min. enzymatic termini: 2 fully, 1 semi, 0 non-enzymatic  [default: 2]
    -c, --missed-cleavages <N>   `-1`/`unlimited` for no limit. MS-GF+ does not limit missed
                           cleavages by default, so the default here matches it; peptide length
                           still bounds the search. Pass `-c 2` for a conventional search (much
                           smaller index).                                   [default: unlimited]
        --min-len <N>      Minimum peptide length                            [default: 6]
        --max-len <N>      Maximum peptide length                            [default: 40]
    -t, --precursor-tol <TOL>   e.g. `10ppm` or `0.5Da`                      [default: 10ppm]
        --ti      <LO,HI>  Isotope-error range, like MS-GF+ -ti              [default: 0,1]
    -n, --num-matches <N>  Matches reported per spectrum                     [default: 1]
        --charges <LO,HI>  Charges tried when the spectrum declares none     [default: 2,3]

MODIFICATIONS:
        --mods    <FILE>   MS-GF+ `Mods.txt`-format configuration file
        --fixed-mod <SPEC> Fixed mod, repeatable. `C+57.021464` or a full
                           `Composition,Residues,fix,Position,Name` spec
        --var-mod <SPEC>   Variable mod, repeatable (same two forms)
        --num-mods <N>     Max variable mods per peptide                     [default: 2]

TARGET-DECOY:
        --tda              Generate reversed decoys before searching (MS-GF+ -tda 1). Not needed
                           if the FASTA is already concatenated (e.g. *.revCat.fasta) — decoys are
                           detected by accession prefix and never regenerated.
        --decoy-prefix <P> Decoy accession prefix                            [default: XXX]

OTHER:
        --db-size <N>      Override the EValue multiplier (default: candidate-index size)
        --threads <N>      Worker threads (default: all cores)
    -h, --help             Print this help

NOTES:
    Amino-acid background frequencies for the generating function are computed from the database
    being searched, so SpecEValue reflects that database rather than a uniform alphabet.
    Q-values come from SpecEValue and require decoys; without them the QValue column is 0 and is
    not an FDR estimate.
";

pub struct Config {
    spectra: PathBuf,
    param: PathBuf,
    fasta: PathBuf,
    out: Option<PathBuf>,
    unroll: bool,
    digest: DigestParams,
    mods: ModSet,
    params: SearchParams,
    tda: bool,
    decoy_prefix: String,
    threads: Option<usize>,
}

impl Config {
    pub fn parse(args: &[String]) -> Result<Config, String> {
        let (mut spectra, mut param, mut fasta, mut out) = (None, None, None, None);
        let mut digest = DigestParams::default();
        let mut params = SearchParams::default();
        let mut mods = ModSet::default();
        let (mut mods_file, mut num_mods): (Option<PathBuf>, Option<usize>) = (None, None);
        let (mut unroll, mut tda, mut threads) = (false, false, None);
        let mut decoy_prefix = String::from("XXX");
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
                "-d" | "--fasta" | "--db" => fasta = Some(PathBuf::from(want("--fasta")?)),
                "-o" | "--out" => out = Some(PathBuf::from(want("--out")?)),
                "--unroll" => unroll = true,
                "-e" | "--enzyme" => digest.enzyme = Enzyme::parse(&want("--enzyme")?)?,
                "--ntt" => {
                    digest.min_termini = want("--ntt")?
                        .parse()
                        .map_err(|_| "--ntt must be 0, 1 or 2")?;
                    if digest.min_termini > 2 {
                        return Err("--ntt must be 0, 1 or 2".into());
                    }
                }
                "-c" | "--missed-cleavages" => {
                    let v = want("--missed-cleavages")?;
                    digest.max_missed_cleavages = if v.trim() == "-1"
                        || v.trim().eq_ignore_ascii_case("unlimited")
                    {
                        msgf_db::enzyme::UNLIMITED_MISSED_CLEAVAGES
                    } else {
                        v.parse().map_err(|_| {
                            "--missed-cleavages must be a non-negative integer, -1, or `unlimited`"
                        })?
                    }
                }
                "--min-len" => {
                    digest.min_len = want("--min-len")?
                        .parse()
                        .map_err(|_| "--min-len must be a positive integer")?
                }
                "--max-len" => {
                    digest.max_len = want("--max-len")?
                        .parse()
                        .map_err(|_| "--max-len must be a positive integer")?
                }
                "-t" | "--precursor-tol" => {
                    params.precursor_tol = parse_tolerance(&want("--precursor-tol")?)?
                }
                "--ti" => params.isotope_errors = parse_pair(&want("--ti")?, "--ti")?,
                "-n" | "--num-matches" => {
                    params.num_matches = want("--num-matches")?
                        .parse()
                        .map_err(|_| "--num-matches must be a positive integer")?
                }
                "--charges" => params.charge_range = parse_pair(&want("--charges")?, "--charges")?,
                "--mods" => mods_file = Some(PathBuf::from(want("--mods")?)),
                "--fixed-mod" => mods
                    .mods
                    .push(ModSpec::parse_short(&want("--fixed-mod")?, true)?),
                "--var-mod" => mods
                    .mods
                    .push(ModSpec::parse_short(&want("--var-mod")?, false)?),
                "--num-mods" => {
                    num_mods = Some(
                        want("--num-mods")?
                            .parse()
                            .map_err(|_| "--num-mods must be a non-negative integer")?,
                    )
                }
                "--tda" => tda = true,
                "--decoy-prefix" => decoy_prefix = want("--decoy-prefix")?,
                "--db-size" => {
                    params.db_size = Some(
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
                "-h" | "--help" => {
                    print!("{USAGE}");
                    std::process::exit(0);
                }
                other => return Err(format!("unexpected argument `{other}`")),
            }
        }

        // A --mods file supplies both the mod list and NumMods; inline flags append to it.
        if let Some(path) = mods_file {
            let from_file = ModSet::read_file(&path)
                .map_err(|e| format!("reading mods {}: {e}", path.display()))?;
            let inline = std::mem::take(&mut mods.mods);
            mods = from_file;
            mods.mods.extend(inline);
        }
        if let Some(n) = num_mods {
            mods.max_var_mods = n;
        }
        if digest.min_len == 0 || digest.min_len > digest.max_len {
            return Err("--min-len must be >= 1 and <= --max-len".into());
        }
        if params.isotope_errors.0 > params.isotope_errors.1 {
            return Err("--ti LO must be <= HI".into());
        }
        if params.charge_range.0 < 1 || params.charge_range.0 > params.charge_range.1 {
            return Err("--charges LO must be >= 1 and <= HI".into());
        }
        if params.num_matches == 0 {
            return Err("--num-matches must be >= 1".into());
        }

        Ok(Config {
            spectra: spectra.ok_or("missing --spectra")?,
            param: param.ok_or("missing --param")?,
            fasta: fasta.ok_or("missing --fasta")?,
            out,
            unroll,
            digest,
            mods,
            params,
            tda,
            decoy_prefix,
            threads,
        })
    }
}

fn parse_tolerance(s: &str) -> Result<Tolerance, String> {
    let t = s.trim();
    let lower = t.to_ascii_lowercase();
    let (num, unit) = if let Some(v) = lower.strip_suffix("ppm") {
        (v, Unit::Ppm)
    } else if let Some(v) = lower.strip_suffix("da") {
        (v, Unit::Da)
    } else {
        return Err(format!("tolerance `{s}` must end in `ppm` or `Da`"));
    };
    let value: f64 = num
        .trim()
        .parse()
        .map_err(|_| format!("tolerance `{s}` has no numeric part"))?;
    if value <= 0.0 {
        return Err(format!("tolerance `{s}` must be positive"));
    }
    Ok(Tolerance { value, unit })
}

fn parse_pair(s: &str, flag: &str) -> Result<(i32, i32), String> {
    let (lo, hi) = s
        .split_once(',')
        .ok_or_else(|| format!("{flag} must be LO,HI"))?;
    Ok((
        lo.trim()
            .parse()
            .map_err(|_| format!("{flag} LO must be an integer"))?,
        hi.trim()
            .parse()
            .map_err(|_| format!("{flag} HI must be an integer"))?,
    ))
}

pub fn run(cfg: &Config) -> Result<(), String> {
    if let Some(n) = cfg.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .map_err(|e| format!("configuring {n} threads: {e}"))?;
    }

    let model = msgf_scorer::read_param_file(&cfg.param)
        .map_err(|e| format!("reading model {}: {e:?}", cfg.param.display()))?;

    let prefix = normalize_prefix(&cfg.decoy_prefix);
    let mut db = ProteinDb::read(&cfg.fasta, &prefix)
        .map_err(|e| format!("reading {}: {e}", cfg.fasta.display()))?;
    if db.proteins.is_empty() {
        return Err(format!("no proteins in {}", cfg.fasta.display()));
    }
    let found_decoys = db.n_decoys();
    let made = if cfg.tda {
        db.add_decoys(DecoyStrategy::Reverse, &prefix)
    } else {
        0
    };
    if cfg.tda && made == 0 && found_decoys > 0 {
        eprintln!("note: --tda ignored, {found_decoys} decoys already present in the database");
    }
    eprintln!(
        "database: {} proteins ({} decoy){}",
        db.proteins.len(),
        db.n_decoys(),
        if made > 0 { ", decoys generated" } else { "" }
    );
    if db.n_decoys() > 0 {
        let accessions: Vec<String> = db.proteins.iter().map(|p| p.name.clone()).collect();
        if let Err(e) = msgf_db::decoy::validate_concatenated(&accessions, &prefix) {
            eprintln!("warning: {e}");
        }
    } else {
        eprintln!(
            "warning: no decoys in the database — q-values will be 0 and are NOT an FDR estimate \
             (pass --tda, or search a concatenated *.revCat.fasta)"
        );
    }

    let index = PeptideIndex::build(&db, &cfg.digest, &cfg.mods);
    eprintln!(
        "index: {} peptides -> {} candidates (with modifications)",
        index.n_peptides,
        index.len()
    );
    if index.truncated_peptides > 0 {
        eprintln!(
            "warning: {} peptide(s) had more modified forms than the per-peptide cap ({}); \
             those extra forms were NOT searched — lower --num-mods to search exhaustively",
            index.truncated_peptides,
            msgf_search::index::MAX_VARIANTS_PER_PEPTIDE
        );
    }
    if index.is_empty() {
        return Err("no candidate peptides — check --enzyme, --min-len/--max-len".into());
    }

    let engine = SearchEngine::new(
        &model,
        &db,
        &index,
        &cfg.mods,
        &cfg.digest,
        cfg.params.clone(),
    );
    for w in engine.warnings() {
        eprintln!("warning: {w}");
    }

    let spectra = msgf_io::read_mgf_file(&cfg.spectra)
        .map_err(|e| format!("reading {}: {e}", cfg.spectra.display()))?;
    if spectra.is_empty() {
        return Err(format!("no spectra in {}", cfg.spectra.display()));
    }
    eprintln!("searching {} spectra...", spectra.len());

    let mut psms = engine.run(&spectra);
    assign_q_values(&mut psms);

    let spec_file = file_name(&cfg.spectra);
    let mut w: Box<dyn Write> = match &cfg.out {
        Some(p) => Box::new(BufWriter::new(
            File::create(p).map_err(|e| format!("creating {}: {e}", p.display()))?,
        )),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };
    let write = if cfg.unroll {
        report::write_tsv_unrolled
    } else {
        report::write_tsv
    };
    write(&mut w, &spec_file, &psms).map_err(|e| format!("writing results: {e}"))?;
    w.flush().map_err(|e| format!("writing results: {e}"))?;

    eprintln!("{}", report::summary(&psms, msgf_search::has_decoys(&psms)));
    Ok(())
}

fn file_name(p: &Path) -> String {
    p.file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| p.display().to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tolerances_parse() {
        let t = parse_tolerance("10ppm").unwrap();
        assert_eq!(t.unit, Unit::Ppm);
        assert!((t.value - 10.0).abs() < 1e-12);
        let t = parse_tolerance("0.5Da").unwrap();
        assert_eq!(t.unit, Unit::Da);
        assert!(parse_tolerance("10").is_err());
        assert!(parse_tolerance("-1ppm").is_err());
    }

    #[test]
    fn pairs_parse() {
        assert_eq!(parse_pair("0,1", "--ti").unwrap(), (0, 1));
        assert_eq!(parse_pair("-1, 2", "--ti").unwrap(), (-1, 2));
        assert!(parse_pair("0", "--ti").is_err());
    }

    fn args(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn minimal_config_parses() {
        let c = Config::parse(&args(&["-s", "a.mgf", "-p", "m.param", "-d", "db.fasta"])).unwrap();
        assert_eq!(c.digest.enzyme.name, "Tryp");
        assert_eq!(c.digest.min_termini, 2);
        assert_eq!(c.params.num_matches, 1);
        assert_eq!(c.mods.max_var_mods, 2);
    }

    #[test]
    fn missing_required_args_are_rejected() {
        assert!(Config::parse(&args(&["-s", "a.mgf"])).is_err());
        assert!(Config::parse(&args(&["--bogus"])).is_err());
    }

    #[test]
    fn inline_mods_and_enzyme_parse() {
        let c = Config::parse(&args(&[
            "-s",
            "a.mgf",
            "-p",
            "m.param",
            "-d",
            "db.fasta",
            "--fixed-mod",
            "C+57.021464",
            "--var-mod",
            "M+15.994915",
            "--num-mods",
            "3",
            "-e",
            "LysC",
            "--ntt",
            "1",
        ]))
        .unwrap();
        assert_eq!(c.mods.mods.len(), 2);
        assert!(c.mods.mods[0].is_fixed);
        assert!(!c.mods.mods[1].is_fixed);
        assert_eq!(c.mods.max_var_mods, 3);
        assert_eq!(c.digest.enzyme.name, "LysC");
        assert_eq!(c.digest.min_termini, 1);
    }

    #[test]
    fn nonsensical_ranges_are_rejected() {
        let base = ["-s", "a.mgf", "-p", "m.param", "-d", "db.fasta"];
        let with = |extra: &[&str]| {
            let mut v: Vec<&str> = base.to_vec();
            v.extend_from_slice(extra);
            Config::parse(&args(&v))
        };
        assert!(with(&["--min-len", "10", "--max-len", "5"]).is_err());
        assert!(with(&["--ti", "2,0"]).is_err());
        assert!(with(&["--charges", "0,3"]).is_err());
        assert!(with(&["--num-matches", "0"]).is_err());
        assert!(with(&["--ntt", "3"]).is_err());
    }
}
