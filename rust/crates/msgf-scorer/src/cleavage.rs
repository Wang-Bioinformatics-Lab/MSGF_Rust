//! The enzymatic-terminus rule shared by RawScore and the null distribution.

/// Enzymatic-terminus scoring rule.
#[derive(Clone, Debug)]
pub struct Cleavage {
    /// Residues after which the enzyme cleaves (C-terminal to them).
    pub sites: Vec<u8>,
    pub credit: i32,
    pub penalty: i32,
    /// False for an unspecific enzyme and for N-terminal-cleaving enzymes (not supported: the
    /// graph is built in the C-terminal direction), which removes every terminus term.
    pub enabled: bool,
}

impl Cleavage {
    /// Trypsin: K/R, credit +2, penalty -11.
    pub fn trypsin() -> Cleavage {
        Cleavage {
            sites: b"KR".to_vec(),
            credit: 2,
            penalty: -11,
            enabled: true,
        }
    }
    /// No terminus scoring.
    pub fn disabled() -> Cleavage {
        Cleavage {
            sites: Vec::new(),
            credit: 0,
            penalty: 0,
            enabled: false,
        }
    }
    /// The terminus term for a peptide whose C-terminal residue is `letter`.
    pub fn term(&self, letter: u8) -> i32 {
        if !self.enabled {
            0
        } else if self.sites.contains(&letter) {
            self.credit
        } else {
            self.penalty
        }
    }
}
