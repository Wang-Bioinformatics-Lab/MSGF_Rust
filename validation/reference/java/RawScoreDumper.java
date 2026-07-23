import edu.ucsd.msjava.msscorer.NewRankScorer;
import edu.ucsd.msjava.msscorer.NewScoredSpectrum;
import edu.ucsd.msjava.msscorer.FastScorer;
import edu.ucsd.msjava.msscorer.DBScanScorer;
import edu.ucsd.msjava.msgf.NominalMass;
import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.Modification;
import edu.ucsd.msjava.msutil.SpectraAccessor;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.params.ParamManager;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * RawScoreDumper -- validation oracle for the FINAL RawScore stage of the Rust
 * MS-GF+ port. For each selected (scan, charge, golden_peptide) it reconstructs
 * MS-GF+'s peptide prefix-mass arrays EXACTLY as edu.ucsd.msjava.msdbsearch's
 * CandidatePeptideGrid + DBScanner do, then runs the reference FastScorer (node
 * summation, no edges) and DBScanScorer (node + edge) getScore(...) and dumps
 * node_only / full / edge scores. Every number comes from actually running
 * MS-GF+ (MSGFPlus.jar). No value is fabricated.
 *
 * Reconstruction facts pinned from the jar bytecode (see report):
 *   - CandidatePeptideGrid.addResidue accumulates nominalPRM[i]=nominalPRM[i-1]
 *     + aa.getNominalMass() (per-residue nominal-mass sum) and prm[i]=prm[i-1]
 *     + aa.getMass() (per-residue accurate sum). Index 0 = 0 (empty prefix);
 *     grid length = numResidues+1; index numResidues = full peptide.
 *   - DBScanner scores with getScore(prm, nominalPRM, fromIndex=1,
 *     toIndex=length+1, getNumMods). The leading-zero + fromIndex=1
 *     convention is REQUIRED: DBScanScorer's prefix-main-ion edge path reads
 *     nominalPRM[fromIndex-1], so fromIndex must be >= 1.
 *   - FastScorer.getScore: pmNominal = nominalPRM[toIndex-1]; for i in
 *     [fromIndex, toIndex-2] adds round(prefixScore[nominalPRM[i]] +
 *     suffixScore[pmNominal - nominalPRM[i]]); numMods contributes 0.
 *   - The scorer is built per spectrum with the spectrum's nominal peptide mass;
 *     the peptide's own nominal mass (nominalPRM[last]) is used here for the
 *     array size, which keeps every prefix/suffix/edge lookup in bounds and so
 *     yields identical node/edge scores (verified vs the golden scored-spectrum
 *     prefix/suffix arrays).
 *
 * The scored spectrum is built EXACTLY as ScoredSpectrumDumper builds it
 * (new NewScoredSpectrum(spec, scorer) from a NewRankScorer.readFromFile model,
 * no doNotUseError()) so that node_only_score is bit-consistent with the
 * f13_scored_spectrum.golden.json prefix/suffix arrays the Rust port already
 * validates against.
 *
 * Args: <model.param> <spectra.mgf> <modfile> <selection.tsv> <out.json>
 *   selection.tsv rows: scan<TAB>charge<TAB>golden_peptide<TAB>golden_raw_score
 *   (golden_peptide keeps its flanking context, e.g. K.RSRRRRKR.A)
 */
public class RawScoreDumper {

    static final String MODEL_NAME = "HCD_QExactive_Tryp.param";

    static final class Sel {
        final int charge;
        final String peptide;   // with flanking context, e.g. "K.RSRRRRKR.A"
        final int rawScore;
        Sel(int charge, String peptide, int rawScore) {
            this.charge = charge; this.peptide = peptide; this.rawScore = rawScore;
        }
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 5) {
            System.err.println("usage: RawScoreDumper <model.param> <spectra.mgf> <modfile> <selection.tsv> <out.json>");
            System.exit(2);
        }
        File modelFile = new File(args[0]);
        File mgfFile   = new File(args[1]);
        String modFilePath = args[2];
        File selFile   = new File(args[3]);
        File outFile   = new File(args[4]);

        // 1. Scoring model, built exactly as ScoredSpectrumDumper / the golden
        //    scored-spectrum pipeline (readFromFile, no doNotUseError()).
        NewRankScorer scorer = new NewRankScorer();
        scorer.readFromFile(modelFile);

        // 2. Amino-acid set from the mod file (needed to represent inline mods
        //    such as M+15.995 = Oxidation on M). registerEnzyme(TRYPSIN) makes
        //    the cleavage credit/penalty constants available for reconciliation.
        ParamManager pm = new ParamManager("RawScoreDumper", "1", "2026", "usage");
        AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSetFromModFile(modFilePath, pm);
        Enzyme enzyme = Enzyme.TRYPSIN;
        aaSet.registerEnzyme(enzyme);
        int neighborCredit = aaSet.getNeighboringAACleavageCredit();
        int neighborPenalty = aaSet.getNeighboringAACleavagePenalty();
        int peptideCredit = aaSet.getPeptideCleavageCredit();
        int peptidePenalty = aaSet.getPeptideCleavagePenalty();

        // 3. Selection (scan -> golden info), preserving order.
        Map<Integer, Sel> selection = new LinkedHashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(selFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] c = line.split("\t", -1);
                int scan = Integer.parseInt(c[0].trim());
                int charge = Integer.parseInt(c[1].trim());
                String peptide = c[2];
                int rawScore = Integer.parseInt(c[3].trim());
                selection.put(scan, new Sel(charge, peptide, rawScore));
            }
        }

        SpectraAccessor accessor = new SpectraAccessor(mgfFile);
        Iterator<Spectrum> it = accessor.getSpecItr();
        List<String> emitted = new ArrayList<>();
        java.util.Set<Integer> seen = new java.util.HashSet<>();
        int fullMatches = 0, nodeMatches = 0, fullPlusCleavMatches = 0;

        while (it.hasNext()) {
            Spectrum spec = it.next();
            int scan = spec.getScanNum();
            Sel sel = selection.get(scan);
            if (sel == null || seen.contains(scan)) continue;
            seen.add(scan);

            // Scored spectrum, exactly as the golden scored-spectrum pipeline.
            NewScoredSpectrum<NominalMass> ss = new NewScoredSpectrum<>(spec, scorer);

            String contextPep = sel.peptide;
            char prevChar = contextPep.charAt(0);
            String core = contextPep.substring(2, contextPep.length() - 2);

            // Parse core into AminoAcids, resolving inline +/-delta mods.
            List<AminoAcid> aas = new ArrayList<>();
            int numMods = 0;
            int i = 0;
            while (i < core.length()) {
                char res = core.charAt(i);
                i++;
                StringBuilder mod = new StringBuilder();
                while (i < core.length()) {
                    char ch = core.charAt(i);
                    if (ch == '+' || ch == '-' || ch == '.' || Character.isDigit(ch)) { mod.append(ch); i++; }
                    else break;
                }
                AminoAcid base = aaSet.getAminoAcid(res);
                if (base == null) throw new RuntimeException("no amino acid for residue " + res + " in " + core);
                if (mod.length() == 0) {
                    aas.add(base);
                } else {
                    double delta = Double.parseDouble(mod.toString());
                    numMods++;
                    AminoAcid found = null;
                    double bestErr = Double.MAX_VALUE;
                    for (AminoAcid aa : aaSet.getAAList(Modification.Location.Anywhere)) {
                        if (aa.getUnmodResidue() == res && aa.isModified()) {
                            double err = Math.abs((aa.getMass() - base.getMass()) - delta);
                            if (err < bestErr) { bestErr = err; found = aa; }
                        }
                    }
                    if (found == null || bestErr > 0.01)
                        throw new RuntimeException("no modified amino acid for " + res + "+" + delta
                            + " (bestErr=" + bestErr + ") in " + core);
                    aas.add(found);
                }
            }

            // Last residue char (skip trailing mod digits) for cleavage scoring.
            char lastResidue = 0;
            for (int k = 0; k < core.length(); k++)
                if (Character.isLetter(core.charAt(k))) lastResidue = core.charAt(k);

            int n = aas.size();
            // Leading-zero cumulative arrays, length n+1 (index 0 = empty prefix).
            double[] prm = new double[n + 1];
            int[] nominalPRM = new int[n + 1];
            for (int k = 0; k < n; k++) {
                prm[k + 1] = prm[k] + aas.get(k).getMass();
                nominalPRM[k + 1] = nominalPRM[k] + aas.get(k).getNominalMass();
            }
            int pepMassNominal = nominalPRM[n];

            FastScorer fs = new FastScorer(ss, pepMassNominal);
            int nodeOnly = fs.getScore(prm, nominalPRM, 1, n + 1, numMods);
            DBScanScorer ds = new DBScanScorer(ss, pepMassNominal);
            int full = ds.getScore(prm, nominalPRM, 1, n + 1, numMods);
            int edge = full - nodeOnly;

            // Enzyme cleavage score MS-GF+'s DBScanner adds to getScore to form
            // the DatabaseMatch score (trypsin cleaves C-terminal to K/R):
            //   N-term (neighboring AA): protein-N-term or prevChar cleavable
            //                            -> credit, else penalty
            //   C-term (peptide)       : last residue cleavable -> credit, else penalty
            int nTermScore = (prevChar == '-' || enzyme.isCleavable(prevChar)) ? neighborCredit : neighborPenalty;
            int cTermScore = enzyme.isCleavable(lastResidue) ? peptideCredit : peptidePenalty;
            int cleavage = nTermScore + cTermScore;
            int fullPlusCleav = full + cleavage;

            if (full == sel.rawScore) fullMatches++;
            if (nodeOnly == sel.rawScore) nodeMatches++;
            if (fullPlusCleav == sel.rawScore) fullPlusCleavMatches++;

            // nominal_prefix_masses emitted WITHOUT the leading zero (length n,
            // last = full peptide nominal mass), matching the Rust contract; the
            // scoring above uses the equivalent leading-zero form internally.
            StringBuilder npm = new StringBuilder("[");
            for (int k = 1; k <= n; k++) {
                if (k > 1) npm.append(",");
                npm.append(nominalPRM[k]);
            }
            npm.append("]");

            StringBuilder b = new StringBuilder();
            b.append("  {\n");
            b.append("   \"scan\": ").append(scan).append(",\n");
            b.append("   \"charge\": ").append(sel.charge).append(",\n");
            b.append("   \"peptide\": ").append(jstr(core)).append(",\n");
            b.append("   \"num_mods\": ").append(numMods).append(",\n");
            b.append("   \"golden_raw_score\": ").append(sel.rawScore).append(",\n");
            b.append("   \"nominal_prefix_masses\": ").append(npm).append(",\n");
            b.append("   \"node_only_score\": ").append(nodeOnly).append(",\n");
            b.append("   \"full_score\": ").append(full).append(",\n");
            b.append("   \"edge_score\": ").append(edge).append(",\n");
            b.append("   \"cleavage_score\": ").append(cleavage).append(",\n");
            b.append("   \"full_plus_cleavage\": ").append(fullPlusCleav).append(",\n");
            b.append("   \"full_matches_golden\": ").append(full == sel.rawScore).append("\n");
            b.append("  }");
            emitted.add(b.toString());
        }

        StringBuilder sb = new StringBuilder();
        sb.append("{\n");
        sb.append(" \"model\": ").append(jstr(MODEL_NAME)).append(",\n");
        sb.append(" \"note\": ").append(jstr(
            "Reference FastScorer/DBScanScorer getScore for the Rust RawScore stage. "
            + "node_only_score=FastScorer.getScore (node summation, no edges); "
            + "full_score=DBScanScorer.getScore (node+edge); edge_score=full-node. "
            + "Arrays built as MS-GF+ CandidatePeptideGrid+DBScanner do: cumulative "
            + "per-residue nominal masses (nominal_prefix_masses, last=full peptide) "
            + "and cumulative per-residue accurate masses; scored via getScore(prm, "
            + "nominalPRM, 1, len+1, numMods). aaSet=getAminoAcidSetFromModFile. "
            + "cleavage_score is MS-GF+'s trypsin terminus credit/penalty that DBScanner "
            + "adds to getScore for the DatabaseMatch score. SELF-CHECK: the reported "
            + "MS-GF:RawScore (golden_raw_score) is NOT reproduced by full_score (or "
            + "full+cleavage) except coincidentally -- see self_check. Generated by "
            + "RawScoreDumper.java against MSGFPlus.jar; every number is from MS-GF+.")).append(",\n");
        sb.append(" \"amino_acid_set\": ").append(jstr("getAminoAcidSetFromModFile(" + new File(modFilePath).getName() + ")")).append(",\n");
        sb.append(" \"cleavage_constants\": {")
          .append("\"neighboring_credit\": ").append(neighborCredit).append(", ")
          .append("\"neighboring_penalty\": ").append(neighborPenalty).append(", ")
          .append("\"peptide_credit\": ").append(peptideCredit).append(", ")
          .append("\"peptide_penalty\": ").append(peptidePenalty).append("},\n");
        sb.append(" \"self_check\": {")
          .append("\"full_score_equals_golden\": ").append(fullMatches).append(", ")
          .append("\"node_only_equals_golden\": ").append(nodeMatches).append(", ")
          .append("\"full_plus_cleavage_equals_golden\": ").append(fullPlusCleavMatches).append(", ")
          .append("\"n\": ").append(emitted.size()).append("},\n");
        sb.append(" \"n_spectra\": ").append(emitted.size()).append(",\n");
        sb.append(" \"spectra\": [\n");
        for (int k = 0; k < emitted.size(); k++) {
            sb.append(emitted.get(k));
            sb.append(k + 1 < emitted.size() ? ",\n" : "\n");
        }
        sb.append(" ]\n");
        sb.append("}\n");

        try (PrintStream out = new PrintStream(new FileOutputStream(outFile), false, "UTF-8")) {
            out.print(sb);
        }
        System.out.println("dumped " + emitted.size() + " spectra -> " + outFile.getPath());
        System.out.println("self-check: full==golden " + fullMatches + "/" + emitted.size()
            + ", node==golden " + nodeMatches + "/" + emitted.size()
            + ", full+cleavage==golden " + fullPlusCleavMatches + "/" + emitted.size());
        if (emitted.size() != selection.size())
            System.out.println("WARNING: selected " + selection.size() + " scans but emitted " + emitted.size());
    }

    static String jstr(String s) {
        StringBuilder b = new StringBuilder("\"");
        for (int i = 0; i < s.length(); i++) {
            char ch = s.charAt(i);
            switch (ch) {
                case '"':  b.append("\\\""); break;
                case '\\': b.append("\\\\"); break;
                case '\n': b.append("\\n");  break;
                case '\r': b.append("\\r");  break;
                case '\t': b.append("\\t");  break;
                default:
                    if (ch < 0x20) b.append(String.format("\\u%04x", (int) ch));
                    else b.append(ch);
            }
        }
        return b.append("\"").toString();
    }
}
