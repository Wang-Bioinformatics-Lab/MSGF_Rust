//! `msgf decoy` — build a target-decoy FASTA, byte-compatible with MS-GF+'s `-tda 1`.

use std::path::PathBuf;

use msgf_db::decoy::{normalize_prefix, write_decoy_database, DecoyOptions, LineSep, Output};

pub const USAGE: &str = "\
msgf decoy — build a target-decoy FASTA (MS-GF+ -tda 1 compatible)

USAGE:
    msgf decoy --fasta <IN.fasta> --out <OUT.revCat.fasta> [OPTIONS]

REQUIRED:
    -d, --fasta   <FILE>   Input protein FASTA (targets)
    -o, --out     <FILE>   Output FASTA

OPTIONS:
        --prefix  <P>      Decoy accession prefix; a trailing `_` is normalised, so `XXX` and
                           `XXX_` are the same                               [default: XXX]
        --decoy-only       Write only the decoys (default: targets then decoys, concatenated)
        --crlf             Use CRLF line endings. MS-GF+ writes the JVM's line separator, so a
                           Windows-generated reference database is CRLF and a Linux one is LF.
    -h, --help             Print this help

Each protein is reversed whole-sequence and emitted as one unwrapped line under
`>` + prefix + the original header. Output is byte-identical to MS-GF+'s own `.revCat.fasta`.
";

pub struct Config {
    fasta: PathBuf,
    out: PathBuf,
    opts: DecoyOptions,
}

impl Config {
    pub fn parse(args: &[String]) -> Result<Config, String> {
        let (mut fasta, mut out) = (None, None);
        let mut opts = DecoyOptions::default();
        let mut it = args.iter();
        while let Some(a) = it.next() {
            let mut want = |name: &str| -> Result<String, String> {
                it.next()
                    .cloned()
                    .ok_or_else(|| format!("`{name}` needs a value"))
            };
            match a.as_str() {
                "-d" | "--fasta" | "--db" => fasta = Some(PathBuf::from(want("--fasta")?)),
                "-o" | "--out" => out = Some(PathBuf::from(want("--out")?)),
                "--prefix" | "--decoy-prefix" => opts.prefix = want("--prefix")?,
                "--decoy-only" => opts.output = Output::DecoyOnly,
                "--crlf" => opts.line_sep = LineSep::Crlf,
                "-h" | "--help" => {
                    print!("{USAGE}");
                    std::process::exit(0);
                }
                other => return Err(format!("unexpected argument `{other}`")),
            }
        }
        Ok(Config {
            fasta: fasta.ok_or("missing --fasta")?,
            out: out.ok_or("missing --out")?,
            opts,
        })
    }
}

pub fn run(cfg: &Config) -> Result<(), String> {
    if cfg.fasta == cfg.out {
        return Err("--fasta and --out must differ (refusing to overwrite the input)".into());
    }
    let n = write_decoy_database(&cfg.fasta, &cfg.out, &cfg.opts)
        .map_err(|e| format!("writing {}: {e}", cfg.out.display()))?;
    eprintln!(
        "wrote {n} decoys with prefix `{}` to {}",
        normalize_prefix(&cfg.opts.prefix),
        cfg.out.display()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn parses_and_normalizes() {
        let c = Config::parse(&args(&[
            "-d", "a.fasta", "-o", "b.fasta", "--prefix", "REV_",
        ]))
        .unwrap();
        assert_eq!(normalize_prefix(&c.opts.prefix), "REV_");
        assert_eq!(c.opts.output, Output::Concatenated);
        assert_eq!(c.opts.line_sep, LineSep::Lf);
    }

    #[test]
    fn requires_both_paths() {
        assert!(Config::parse(&args(&["-d", "a.fasta"])).is_err());
        assert!(Config::parse(&args(&["-o", "b.fasta"])).is_err());
    }

    #[test]
    fn refuses_to_overwrite_the_input() {
        let c = Config::parse(&args(&["-d", "a.fasta", "-o", "a.fasta"])).unwrap();
        assert!(run(&c).is_err());
    }
}
