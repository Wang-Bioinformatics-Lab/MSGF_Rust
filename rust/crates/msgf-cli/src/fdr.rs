//! `msgf fdr` — append MS-GF+-compatible `QValue` / `PepQValue` columns to an existing PSM table.
//!
//! Useful without a search engine: rescore a target+decoy PSM list, then run this over the result.
//! Input is a TSV with a header; the peptide, protein and SpecEValue columns are located by name
//! (MS-GF+'s own names are the defaults) and every original column is preserved.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

use msgf_db::decoy::normalize_prefix;
use msgf_fdr::{is_decoy_match, peptide_key, PsmRecord, TargetDecoyAnalysis};

pub const USAGE: &str = "\
msgf fdr — add target-decoy QValue / PepQValue columns to a PSM table

USAGE:
    msgf fdr --psms <IN.tsv> [OPTIONS]

REQUIRED:
    -i, --psms    <FILE>   PSM table (TSV with a header row)

OPTIONS:
    -o, --out     <FILE>   Output TSV (default: stdout)
        --decoy-prefix <P> Decoy accession prefix                            [default: XXX]
        --score-col <NAME> Score column, smaller-is-better       [default: SpecEValue]
        --peptide-col <NAME>                                     [default: Peptide]
        --protein-col <NAME>                                     [default: Protein]
    -h, --help             Print this help

Rows sharing a (peptide, score) are one match with several protein occurrences (MS-GF+ `-unroll 1`
output): they are rolled up before the sweep, since FDR counts matches, not occurrences. A match is
a decoy only when EVERY one of its protein occurrences carries the decoy prefix.
";

pub struct Config {
    psms: PathBuf,
    out: Option<PathBuf>,
    decoy_prefix: String,
    score_col: String,
    peptide_col: String,
    protein_col: String,
}

impl Config {
    pub fn parse(args: &[String]) -> Result<Config, String> {
        let (mut psms, mut out) = (None, None);
        let mut decoy_prefix = String::from("XXX");
        let mut score_col = String::from("SpecEValue");
        let mut peptide_col = String::from("Peptide");
        let mut protein_col = String::from("Protein");
        let mut it = args.iter();
        while let Some(a) = it.next() {
            let mut want = |name: &str| -> Result<String, String> {
                it.next()
                    .cloned()
                    .ok_or_else(|| format!("`{name}` needs a value"))
            };
            match a.as_str() {
                "-i" | "--psms" => psms = Some(PathBuf::from(want("--psms")?)),
                "-o" | "--out" => out = Some(PathBuf::from(want("--out")?)),
                "--decoy-prefix" => decoy_prefix = want("--decoy-prefix")?,
                "--score-col" => score_col = want("--score-col")?,
                "--peptide-col" => peptide_col = want("--peptide-col")?,
                "--protein-col" => protein_col = want("--protein-col")?,
                "-h" | "--help" => {
                    print!("{USAGE}");
                    std::process::exit(0);
                }
                other => return Err(format!("unexpected argument `{other}`")),
            }
        }
        Ok(Config {
            psms: psms.ok_or("missing --psms")?,
            out,
            decoy_prefix,
            score_col,
            peptide_col,
            protein_col,
        })
    }
}

pub fn run(cfg: &Config) -> Result<(), String> {
    let file = File::open(&cfg.psms).map_err(|e| format!("opening {}: {e}", cfg.psms.display()))?;
    let mut lines = BufReader::new(file).lines();
    let header = lines
        .next()
        .transpose()
        .map_err(|e| format!("reading {}: {e}", cfg.psms.display()))?
        .ok_or_else(|| format!("{} is empty", cfg.psms.display()))?;

    let cols: Vec<String> = header
        .split('\t')
        .map(|c| c.trim().trim_start_matches('#').to_string())
        .collect();
    let find = |name: &str| -> Result<usize, String> {
        cols.iter()
            .position(|c| c.eq_ignore_ascii_case(name))
            .ok_or_else(|| {
                format!(
                    "no `{name}` column in {} (found: {})",
                    cfg.psms.display(),
                    cols.join(", ")
                )
            })
    };
    let (si, pi, ri) = (
        find(&cfg.score_col)?,
        find(&cfg.peptide_col)?,
        find(&cfg.protein_col)?,
    );

    let rows: Vec<String> = lines
        .collect::<Result<_, _>>()
        .map_err(|e| format!("reading {}: {e}", cfg.psms.display()))?;

    // Roll protein occurrences up into matches, keyed by everything except the protein column.
    let prefix = normalize_prefix(&cfg.decoy_prefix);
    let mut order: Vec<String> = Vec::new();
    let mut by_match: HashMap<String, (f32, String, Vec<String>)> = HashMap::new();
    let mut row_key: Vec<Option<String>> = Vec::with_capacity(rows.len());

    for row in &rows {
        if row.trim().is_empty() {
            row_key.push(None);
            continue;
        }
        let f: Vec<&str> = row.split('\t').collect();
        let (Some(score), Some(peptide), Some(protein)) = (f.get(si), f.get(pi), f.get(ri)) else {
            row_key.push(None);
            continue;
        };
        let Ok(score) = score.trim().parse::<f32>() else {
            row_key.push(None);
            continue;
        };
        let key: String = f
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != ri)
            .map(|(_, v)| *v)
            .collect::<Vec<_>>()
            .join("\t");
        by_match
            .entry(key.clone())
            .and_modify(|m| m.2.push(protein.trim().to_string()))
            .or_insert_with(|| {
                order.push(key.clone());
                (
                    score,
                    peptide_key(peptide.trim()),
                    vec![protein.trim().to_string()],
                )
            });
        row_key.push(Some(key));
    }
    if by_match.is_empty() {
        return Err(format!("no usable rows in {}", cfg.psms.display()));
    }

    let records: Vec<PsmRecord> = order
        .iter()
        .map(|k| {
            let (score, peptide, proteins) = &by_match[k];
            PsmRecord {
                score: *score,
                peptide: peptide.clone(),
                is_decoy: is_decoy_match(proteins.iter().map(String::as_str), &prefix),
            }
        })
        .collect();
    let n_decoy = records.iter().filter(|r| r.is_decoy).count();
    if n_decoy == 0 {
        eprintln!(
            "warning: no decoy matches found with prefix `{prefix}` — q-values will be 0 and are \
             NOT an FDR estimate"
        );
    }
    let tda = TargetDecoyAnalysis::new(&records, 1.0);

    let mut w: Box<dyn Write> = match &cfg.out {
        Some(p) => Box::new(BufWriter::new(
            File::create(p).map_err(|e| format!("creating {}: {e}", p.display()))?,
        )),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };
    let io = |e: std::io::Error| format!("writing results: {e}");
    writeln!(w, "{header}\tQValue\tPepQValue").map_err(io)?;
    for (row, key) in rows.iter().zip(&row_key) {
        match key {
            Some(k) => {
                let (score, peptide, _) = &by_match[k];
                writeln!(
                    w,
                    "{row}\t{}\t{}",
                    tda.psm_q_value(*score),
                    tda.pep_q_value(peptide, *score)
                )
                .map_err(io)?;
            }
            None if row.trim().is_empty() => {}
            None => writeln!(w, "{row}\tNA\tNA").map_err(io)?,
        }
    }
    w.flush().map_err(io)?;
    eprintln!(
        "{} match(es) from {} row(s); {n_decoy} decoy",
        records.len(),
        rows.len()
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
    fn defaults_match_msgfplus_column_names() {
        let c = Config::parse(&args(&["-i", "a.tsv"])).unwrap();
        assert_eq!(c.score_col, "SpecEValue");
        assert_eq!(c.peptide_col, "Peptide");
        assert_eq!(c.protein_col, "Protein");
        assert_eq!(normalize_prefix(&c.decoy_prefix), "XXX_");
    }

    #[test]
    fn requires_psms() {
        assert!(Config::parse(&args(&[])).is_err());
    }

    #[test]
    fn rolls_up_and_annotates() {
        let dir = std::env::temp_dir().join("msgf-cli-tests");
        std::fs::create_dir_all(&dir).unwrap();
        let inp = dir.join("fdr_in.tsv");
        let out = dir.join("fdr_out.tsv");
        // Two protein occurrences of one match, plus a decoy match.
        std::fs::write(
            &inp,
            "#SpecFile\tPeptide\tProtein\tSpecEValue\n\
             a.mgf\tK.SAMPLER.A\tP1\t1e-10\n\
             a.mgf\tK.SAMPLER.A\tP2\t1e-10\n\
             a.mgf\tK.DECOYK.A\tXXX_P9\t1e-5\n",
        )
        .unwrap();
        let cfg = Config::parse(&args(&[
            "-i",
            inp.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ]))
        .unwrap();
        run(&cfg).unwrap();

        let text = std::fs::read_to_string(&out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert!(lines[0].ends_with("QValue\tPepQValue"));
        assert_eq!(lines.len(), 4); // header + 3 rows, all preserved
                                    // The target match beats the only decoy, so it sits at q = 0.
        assert!(lines[1].ends_with("\t0\t0"), "{}", lines[1]);
        assert!(lines[2].ends_with("\t0\t0"), "{}", lines[2]);
    }
}
