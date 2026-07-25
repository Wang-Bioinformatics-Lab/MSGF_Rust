//! End-to-end search validation.
//!
//! Both tests read the gitignored `validation/data/` (UC-licensed, see `CLAUDE.md`) and **skip
//! gracefully** when it is absent, so a fresh clone still passes `cargo test`.
//!
//! - [`search_recovers_known_peptides_with_exact_scores`] is a fast check that plants peptides
//!   MS-GF+ reported into a synthetic database and requires the search to rediscover them with
//!   MS-GF+'s exact scores.
//! - [`f13_search_matches_msgfplus`] is the real oracle comparison against MS-GF+'s own F13 output.
//!   It builds a ~48M-candidate index over the human database (~2 GB, ~10 s), so it is `#[ignore]`d
//!   and run explicitly:
//!
//!   ```text
//!   cargo test -p msgf-search --release -- --ignored --nocapture
//!   ```

use msgf_db::enzyme::{DigestParams, Enzyme};
use msgf_db::fasta::{ProteinDb, DEFAULT_DECOY_PREFIX};
use msgf_search::index::PeptideIndex;
use msgf_search::mods::{ModSet, ModSpec};
use msgf_search::{assign_q_values, SearchEngine, SearchParams};
use std::collections::HashMap;
use std::path::PathBuf;

fn repo(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

/// The F13 search configuration MS-GF+ used: `-inst 1 -m 3 -e 1 -t 10ppm -tda 1` with
/// `iprg-2013_Mods.txt` (NumMods=2, oxidation on M, no fixed mods) and MS-GF+'s defaults for
/// everything else — which notably means **unlimited missed cleavages**.
fn f13_config() -> (DigestParams, ModSet, SearchParams) {
    let digest = DigestParams {
        enzyme: Enzyme::builtin(1).unwrap(),
        ..Default::default()
    };
    let mods = ModSet {
        mods: vec![ModSpec::parse("O1,M,opt,any,Oxidation").unwrap()],
        max_var_mods: 2,
    };
    let params = SearchParams {
        precursor_tol: msgf_chem::Tolerance::ppm(10.0),
        isotope_errors: (0, 1),
        num_matches: 1,
        ..Default::default()
    };
    (digest, mods, params)
}

/// A fast, self-contained end-to-end check that still has real teeth.
///
/// The contaminant database shares nothing with the human F13 spectra, so searching it finds
/// nothing and proves nothing. Instead this builds a **synthetic protein** out of peptides MS-GF+
/// itself reported for these scans (chosen so both termini stay tryptic, keeping the cleavage
/// score identical to the original protein context), then requires the search to rediscover them
/// with MS-GF+'s exact RawScore and DeNovoScore.
#[test]
fn search_recovers_known_peptides_with_exact_scores() {
    let mgf = repo("validation/data/spectra/F13.mgf");
    let param = repo("validation/data/models/HCD_HighRes_Tryp.param");
    let golden = repo("validation/golden/iprg2013_F13.tsv");
    if !mgf.exists() || !param.exists() || !golden.exists() {
        eprintln!("skip: validation/data absent");
        return;
    }

    // Best MS-GF+ PSM per scan, keeping only fully-tryptic unmodified peptides so a synthetic
    // protein reproduces their flanking context exactly.
    let text = std::fs::read_to_string(&golden).unwrap();
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let col = |n: &str| header.iter().position(|c| *c == n).expect(n);
    let (c_scan, c_pep, c_raw, c_dn, c_se) = (
        col("ScanNum"),
        col("Peptide"),
        col("MSGFScore"),
        col("DeNovoScore"),
        col("SpecEValue"),
    );
    let mut best: HashMap<String, (f64, String, i32, i32)> = HashMap::new();
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() <= c_se {
            continue;
        }
        let e: f64 = f[c_se].parse().unwrap();
        let v = (
            e,
            f[c_pep].to_string(),
            f[c_raw].parse().unwrap(),
            f[c_dn].parse().unwrap(),
        );
        best.entry(f[c_scan].to_string())
            .and_modify(|cur| {
                if e < cur.0 {
                    *cur = v.clone()
                }
            })
            .or_insert(v);
    }

    let strip = |p: &str| -> String {
        let b = p.as_bytes();
        if b.len() >= 4 && b[1] == b'.' && b[b.len() - 2] == b'.' {
            p[2..p.len() - 2].to_string()
        } else {
            p.to_string()
        }
    };
    // Fully tryptic (N-flank K/R, last residue K/R), no modifications, so concatenating them into
    // one protein preserves every peptide's cleavage context.
    let mut chosen: Vec<(String, String, i32, i32)> = best
        .iter()
        .filter(|(_, (_, pep, _, _))| {
            let core = strip(pep);
            let n_ok = matches!(pep.as_bytes()[0], b'K' | b'R');
            let c_ok = matches!(core.as_bytes()[core.len() - 1], b'K' | b'R');
            n_ok && c_ok && !core.contains('+') && !core.contains('-') && core.len() >= 7
        })
        .map(|(scan, (_, pep, raw, dn))| (scan.clone(), strip(pep), *raw, *dn))
        .collect();
    chosen.sort();
    chosen.truncate(40);
    assert!(
        chosen.len() >= 10,
        "not enough usable golden peptides: {}",
        chosen.len()
    );

    // One protein: every peptide ends in K/R, so each is flanked by the previous peptide's K/R.
    let protein: String = chosen.iter().map(|c| c.1.as_str()).collect();
    let dir = std::env::temp_dir().join("msgf-search-tests");
    std::fs::create_dir_all(&dir).unwrap();
    let fasta = dir.join("synthetic.fasta");
    std::fs::write(&fasta, format!(">SYN synthetic\n{protein}\n")).unwrap();

    let model = msgf_scorer::read_param_file(&param).unwrap();
    let db = ProteinDb::read(&fasta, DEFAULT_DECOY_PREFIX).unwrap();
    let (digest, mods, mut params) = f13_config();
    // Concatenating peptides creates junction sequences that are themselves valid tryptic
    // candidates, and one can outscore a planted peptide for its own scan. Report the top few so
    // the assertion is "the planted peptide is found and scored exactly", not "it happens to win
    // against artifacts the real database would never contain".
    params.num_matches = 5;
    let index = PeptideIndex::build(&db, &digest, &mods);
    assert!(!index.is_empty());

    let engine = SearchEngine::new(&model, &db, &index, &mods, &digest, params);
    let spectra = msgf_io::read_mgf_file(&mgf).unwrap();
    let wanted: HashMap<&str, &(String, String, i32, i32)> =
        chosen.iter().map(|c| (c.0.as_str(), c)).collect();
    let subset: Vec<msgf_io::Spectrum> = spectra
        .into_iter()
        .filter(|s| s.scan.as_deref().is_some_and(|n| wanted.contains_key(n)))
        .collect();
    assert_eq!(
        subset.len(),
        chosen.len(),
        "not every chosen scan was found in the MGF"
    );

    let mut psms = engine.run(&subset);
    assign_q_values(&mut psms);

    let (mut found, mut raw_ok, mut dn_ok) = (0, 0, 0);
    let mut bad = Vec::new();
    let mut seen: std::collections::HashSet<&str> = std::collections::HashSet::new();
    for p in &psms {
        let Some(want) = wanted.get(p.scan.as_str()) else {
            continue;
        };
        if strip(&p.peptide) != want.1 || !seen.insert(p.scan.as_str()) {
            continue;
        }
        found += 1;
        if p.raw_score == want.2 {
            raw_ok += 1;
        } else {
            bad.push(format!(
                "scan {} {}: RawScore java {} rust {}",
                p.scan, want.1, want.2, p.raw_score
            ));
        }
        if p.denovo_score == want.3 {
            dn_ok += 1;
        }
        assert!(!p.proteins.is_empty());
        assert!(p.q_value.is_finite() && (0.0..=1.0).contains(&p.q_value));
        assert!(!p.peptide_key.contains('.'));
    }
    eprintln!("synthetic-DB search: recovered {found}/{} planted peptides; RawScore {raw_ok} exact, DeNovoScore {dn_ok} exact",
              chosen.len());
    assert_eq!(
        found,
        chosen.len(),
        "the search must rediscover every planted peptide"
    );
    assert_eq!(raw_ok, found, "RawScore must match MS-GF+ exactly: {bad:?}");
    assert_eq!(dn_ok, found, "DeNovoScore must match MS-GF+ exactly");
}

/// The Phase-6 oracle comparison. MS-GF+'s own F13 output is the reference.
///
/// What is asserted, and why these are the right bars:
///
/// - **No scan where our best candidate scores lower than MS-GF+'s.** This is the candidate-
///   generation gate: a lower top score means we failed to generate a peptide MS-GF+ found.
/// - **On every scan where we agree on the peptide, RawScore and DeNovoScore match exactly** and
///   SpecEValue is within the `|log10(rust/java)| <= 0.05` contract from `CLAUDE.md`.
///
/// Not asserted: that we always pick the *same* peptide. ~7% of scans are exact score ties between
/// isobaric alternatives (`R.RLTALR.G` vs `R.RIVVSR.G`, both RawScore 16 and the same SpecEValue),
/// where the choice is an arbitrary tie-break. A further ~1% are scans where we score *higher* than
/// MS-GF+ — consistent with its two-stage `FastScorer` pre-filter dropping candidates before full
/// node+edge scoring, which we do not replicate. See `plans/PLAN2.md` §4 for why F13 cannot support an
/// ID-count-at-1%-FDR gate.
#[test]
#[ignore = "builds a ~48M-candidate index over the human database (~2 GB, ~10 s)"]
fn f13_search_matches_msgfplus() {
    let mgf = repo("validation/data/spectra/F13.mgf");
    let param = repo("validation/data/models/HCD_HighRes_Tryp.param");
    let fasta = repo("validation/data/fasta/iprg2013_human.revCat.fasta");
    let golden = repo("validation/golden/iprg2013_F13.tsv");
    if !mgf.exists() || !param.exists() || !fasta.exists() || !golden.exists() {
        eprintln!("skip: validation/data absent (run fetch_reference_data.sh --full)");
        return;
    }

    // MS-GF+'s best PSM per scan.
    struct Ref {
        spec_evalue: f64,
        peptide: String,
        raw: i32,
        denovo: i32,
    }
    let mut want: HashMap<String, Ref> = HashMap::new();
    let text = std::fs::read_to_string(&golden).unwrap();
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let col = |n: &str| header.iter().position(|c| *c == n).expect(n);
    let (c_scan, c_pep, c_raw, c_dn, c_se) = (
        col("ScanNum"),
        col("Peptide"),
        col("MSGFScore"),
        col("DeNovoScore"),
        col("SpecEValue"),
    );
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() <= c_se {
            continue;
        }
        let e: f64 = f[c_se].parse().unwrap();
        let entry = want.entry(f[c_scan].to_string());
        let r = Ref {
            spec_evalue: e,
            peptide: f[c_pep].to_string(),
            raw: f[c_raw].parse().unwrap(),
            denovo: f[c_dn].parse().unwrap(),
        };
        entry
            .and_modify(|cur| {
                if e < cur.spec_evalue {
                    *cur = Ref {
                        spec_evalue: e,
                        peptide: f[c_pep].to_string(),
                        raw: f[c_raw].parse().unwrap(),
                        denovo: f[c_dn].parse().unwrap(),
                    };
                }
            })
            .or_insert(r);
    }
    assert!(
        want.len() > 1000,
        "golden TSV looks wrong: {} scans",
        want.len()
    );

    let model = msgf_scorer::read_param_file(&param).unwrap();
    let db = ProteinDb::read(&fasta, DEFAULT_DECOY_PREFIX).unwrap();
    let (digest, mods, params) = f13_config();
    let index = PeptideIndex::build(&db, &digest, &mods);
    let engine = SearchEngine::new(&model, &db, &index, &mods, &digest, params);
    let spectra = msgf_io::read_mgf_file(&mgf).unwrap();
    let mut psms = engine.run(&spectra);
    assign_q_values(&mut psms);

    // Our best PSM per scan.
    let mut got: HashMap<&str, &msgf_search::Psm> = HashMap::new();
    for p in &psms {
        got.entry(p.scan.as_str())
            .and_modify(|cur| {
                if p.spec_evalue < cur.spec_evalue {
                    *cur = p;
                }
            })
            .or_insert(p);
    }

    let strip = |p: &str| -> String {
        let b = p.as_bytes();
        if b.len() >= 4 && b[1] == b'.' && b[b.len() - 2] == b'.' {
            p[2..p.len() - 2].to_string()
        } else {
            p.to_string()
        }
    };

    let (mut common, mut lower, mut agreed, mut raw_ok, mut dn_ok, mut se_ok) = (0, 0, 0, 0, 0, 0);
    let mut lower_examples = Vec::new();
    let mut score_bad = Vec::new();
    for (scan, r) in &want {
        let Some(p) = got.get(scan.as_str()) else {
            continue;
        };
        common += 1;
        if p.raw_score < r.raw {
            lower += 1;
            if lower_examples.len() < 5 {
                lower_examples.push(format!(
                    "scan {scan}: java {} raw {} | rust {} raw {}",
                    r.peptide, r.raw, p.peptide, p.raw_score
                ));
            }
        }
        if strip(&r.peptide) == strip(&p.peptide) {
            agreed += 1;
            if p.raw_score == r.raw {
                raw_ok += 1;
            } else if score_bad.len() < 5 {
                score_bad.push(format!(
                    "scan {scan} {}: RawScore java {} rust {}",
                    r.peptide, r.raw, p.raw_score
                ));
            }
            if p.denovo_score == r.denovo {
                dn_ok += 1;
            }
            if (p.spec_evalue / r.spec_evalue).log10().abs() <= 0.05 {
                se_ok += 1;
            }
        }
    }

    eprintln!(
        "F13 search vs MS-GF+: {common} common scans, {agreed} identical top peptides \
         ({:.1}%); on those RawScore {raw_ok}/{agreed} exact, DeNovoScore {dn_ok}/{agreed} exact, \
         SpecEValue {se_ok}/{agreed} in tolerance; {lower} scans scored lower than MS-GF+",
        100.0 * agreed as f64 / common as f64
    );

    assert!(common > 1000, "only {common} scans in common");
    assert_eq!(
        lower, 0,
        "our best candidate scored lower than MS-GF+'s on {lower} scan(s) — a candidate-generation \
         gap. Examples: {lower_examples:?}"
    );
    assert_eq!(
        raw_ok, agreed,
        "RawScore must be exact where the peptide agrees. Mismatches: {score_bad:?}"
    );
    assert_eq!(
        dn_ok, agreed,
        "DeNovoScore must be exact where the peptide agrees"
    );
    assert_eq!(
        se_ok, agreed,
        "SpecEValue must be within |log10 ratio| <= 0.05"
    );
    // The corpus identifies essentially nothing (PLAN2 §4); MS-GF+ itself finds exactly one PSM
    // at 1% FDR. Reproducing that count is the only ID-level check F13 can support.
    assert_eq!(
        msgf_search::n_targets_below(&psms, 0.01),
        1,
        "MS-GF+ reports exactly one target PSM at 1% FDR on F13"
    );
}
