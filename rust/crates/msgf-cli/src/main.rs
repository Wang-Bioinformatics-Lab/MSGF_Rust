//! `msgf` — command-line front-end for MSGF_Rust.
//!
//! Four subcommands:
//!
//! | Command | What it does |
//! |---|---|
//! | `search` | full database search: FASTA → candidates → RawScore/SpecEValue → q-values |
//! | `rescore` | recompute MS-GF+ scores for a supplied PSM list (no candidate generation) |
//! | `decoy` | build a target-decoy FASTA, byte-compatible with MS-GF+ `-tda 1` |
//! | `fdr` | append `QValue`/`PepQValue` to an existing PSM table |
//!
//! The library behind them is the `msgf` crate (or the individual `msgf-*` crates).

mod decoy;
mod fdr;
mod model;
mod rescore;
mod search;

use std::process::ExitCode;

const USAGE: &str = "\
msgf — MSGF_Rust command line

USAGE:
    msgf <COMMAND> [OPTIONS]

COMMANDS:
    search    Search MS/MS spectra against a protein database (FASTA), with target-decoy FDR
    rescore   Recompute RawScore / DeNovoScore / SpecEValue for a supplied PSM list
    decoy     Build a target-decoy FASTA (MS-GF+ -tda 1 compatible)
    fdr       Add QValue / PepQValue columns to an existing PSM table

Run `msgf <COMMAND> --help` for that command's options.

Scoring uses the bundled HCD/HighRes/Tryptic model (trained from the CC0 MassIVE-KB corpus)
unless `--param <MODEL.param>` names another one. To reproduce or diff against MS-GF+'s own
output, pass MS-GF+'s .param — the bundled model is a different scoring function and gives
different peptides and SpecEValues by construction.

EXAMPLES:
    # Search a concatenated target-decoy database (no --param: uses the bundled model)
    msgf search -s run.mgf -d human.revCat.fasta \\
        --fixed-mod C+57.021464 --var-mod M+15.994915 -t 10ppm -o psms.tsv

    # Or let it build the decoys first
    msgf decoy -d human.fasta -o human.revCat.fasta
    msgf search -s run.mgf -d human.revCat.fasta -o psms.tsv

    # Rescore an existing PSM list, then add q-values
    msgf rescore -s run.mgf -i psms.tsv -o rescored.tsv
    msgf fdr -i rescored.tsv -o rescored.q.tsv
";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let rest = args.get(1..).unwrap_or(&[]);
    match args.first().map(String::as_str) {
        Some("search") => dispatch(
            "search",
            search::Config::parse(rest),
            search::run,
            search::USAGE,
        ),
        Some("rescore") => dispatch(
            "rescore",
            rescore::Config::parse(rest),
            rescore::run,
            rescore::USAGE,
        ),
        Some("decoy") => dispatch(
            "decoy",
            decoy::Config::parse(rest),
            decoy::run,
            decoy::USAGE,
        ),
        Some("fdr") => dispatch("fdr", fdr::Config::parse(rest), fdr::run, fdr::USAGE),
        Some("-h") | Some("--help") | Some("help") => {
            print!("{USAGE}");
            ExitCode::SUCCESS
        }
        Some("-V") | Some("--version") => {
            println!("msgf {}", env!("CARGO_PKG_VERSION"));
            ExitCode::SUCCESS
        }
        Some(other) => {
            eprintln!("msgf: unknown subcommand `{other}`\n\n{USAGE}");
            ExitCode::FAILURE
        }
        None => {
            print!("{USAGE}");
            ExitCode::FAILURE
        }
    }
}

/// Parse-then-run, reporting configuration errors with that subcommand's usage and runtime errors
/// on their own (a failing search should not bury its message under 60 lines of help text).
fn dispatch<C>(
    name: &str,
    cfg: Result<C, String>,
    run: fn(&C) -> Result<(), String>,
    usage: &str,
) -> ExitCode {
    match cfg {
        Ok(cfg) => match run(&cfg) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("msgf {name}: {e}");
                ExitCode::FAILURE
            }
        },
        Err(e) => {
            eprintln!("msgf {name}: {e}\n\n{usage}");
            ExitCode::FAILURE
        }
    }
}
