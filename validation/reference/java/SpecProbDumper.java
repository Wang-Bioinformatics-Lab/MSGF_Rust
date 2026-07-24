import edu.ucsd.msjava.msscorer.NewRankScorer;
import edu.ucsd.msjava.msscorer.NewScorerFactory;
import edu.ucsd.msjava.msscorer.NewScoredSpectrum;
import edu.ucsd.msjava.msscorer.DBScanScorer;
import edu.ucsd.msjava.msgf.DeNovoGraph;
import edu.ucsd.msjava.msgf.FlexAminoAcidGraph;
import edu.ucsd.msjava.msgf.GeneratingFunction;
import edu.ucsd.msjava.msgf.GeneratingFunctionGroup;
import edu.ucsd.msjava.msgf.NominalMass;
import edu.ucsd.msjava.msgf.ScoreDist;
import edu.ucsd.msjava.msgf.Tolerance;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.AminoAcid;
import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Composition;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.InstrumentType;
import edu.ucsd.msjava.msutil.Modification;
import edu.ucsd.msjava.msutil.Protocol;
import edu.ucsd.msjava.msutil.SpectraAccessor;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.params.ParamManager;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * SpecProbDumper -- validation oracle for MS-GF+'s GENERATING FUNCTION stage
 * (MS-GF:DeNovoScore and the spectral probability / MS-GF:SpecEValue) for the
 * Rust MS-GF+ port.
 *
 * For each frozen (scan, charge, peptide) PSM it MIRRORS EXACTLY the construction
 * that edu.ucsd.msjava.msdbsearch.DBScanner.computeSpecEValue uses to compute the
 * reported MS-GF:SpecEValue / MS-GF:DeNovoScore in the real search, then dumps the
 * generating-function max score (DeNovoScore), the spectral probability, and (for
 * a few spectra) the full ScoreDist so the Rust DP can be validated bin-by-bin.
 * Every number comes from actually running MSGFPlus.jar. No value is fabricated.
 *
 * ---------------------------------------------------------------------------
 * GF CONSTRUCTION -- pinned from DBScanner.java / ScoredSpectraMap (bytecode),
 * reproducing the search `-inst 1 -m 3 -e 1 -t 10ppm -tda 1` (default -ti 0,1):
 *
 *   scorer  = NewScorerFactory.get(HCD, HIGH_RESOLUTION_LTQ, TRYPSIN, STANDARD); // -inst 1; edge scoring ON
 *   spec.setCharge(charge);
 *   NewScoredSpectrum<NominalMass> ss = scorer.getScoredSpectrum(spec);
 *   float  peptideMass         = ss.getPrecursorPeak().getMass() - Composition.H2O;
 *   int    nominalPeptideMass  = NominalMass.toNominalMass(peptideMass);
 *   // DBScanScorer built with the precursor-derived nominal mass (same as ScoredSpectraMap):
 *   int    dbNominal           = nominalPeptideMass + round(tolDaLeft-0.4999) - minIso;
 *   DBScanScorer scoredSpec    = new DBScanScorer(ss, dbNominal);   // node+edge scorer
 *
 *   // GeneratingFunctionGroup over the isotope+tolerance mass-index range:
 *   minPeptideMassIndex = (nominalPeptideMass - maxIso) - round(tolDaRight-0.4999);
 *   maxPeptideMassIndex = (nominalPeptideMass - minIso) + round(tolDaLeft -0.4999);
 *   for (idx = min..max):
 *       graph = new FlexAminoAcidGraph(aaSet, idx, TRYPSIN, scoredSpec, false, false);
 *       gfi   = new GeneratingFunction(graph).doNotBacktrack().doNotCalcNumber();
 *       gf.registerGF(graph.getPMNode(), gfi);
 *   gf.computeGeneratingFunction();
 *   deNovoScore = gf.getMaxScore() - 1;                 // == TSV MS-GF:DeNovoScore
 *   specProb    = gf.getSpectralProbability(rawScore);  // == TSV MS-GF:SpecEValue
 *
 * where rawScore is the DBScanner match score for this peptide:
 *   rawScore = cleavageScore + scoredSpec.getScore(prm, nominalPRM, 1, len+1, numMods)
 * with cleavageScore the trypsin terminus credit/penalty (identical to RawScoreDumper).
 *
 * NOTE 1 (edge scoring ON): FlexAminoAcidGraph reads BOTH per-node scores
 *   (scoredSpec.getNodeScore) AND per-edge error scores (scoredSpec.getEdgeScore),
 *   so the graph and its ScoreDist are the WITH-EDGES construction.
 * NOTE 2 (getMaxScore()-1): ScoreDist's maxScore is exclusive; the top achievable
 *   score is getMaxScore()-1, which is exactly what DBScanner stores as DeNovoScore.
 * NOTE 3 (setUpScoreThreshold): DBScanner additionally calls setUpScoreThreshold(minScore)
 *   as a pruning optimization. It does NOT change getMaxScore() nor
 *   getSpectralProbability(score>=minScore); we OMIT it so getScoreDist() returns the
 *   FULL distribution for the Rust DP to validate. (SpecEValue is reproduced identically;
 *   verified by the self_check counters below.)
 * NOTE 4 (no doNotUseError): the DEFAULT high-res search keeps edge scoring ON. doNotUseError()
 *   is only called under -turnOffEdgeScoring; it zeroes errorScalingFactor and makes
 *   supportEdgeScores() false (node-only FastScorer). We do NOT call it. The residual gap the
 *   earlier rawscore golden saw (reconstructed full_score != reported MSGFScore) came instead from
 *   sizing the scorer with the peptide's own nominal mass; the real search sizes it with the
 *   PRECURSOR-derived nominal mass (maxNominalPeptideMass), which this driver reproduces.
 *
 * Args: <model.param> <spectra.mgf> <modfile> <selection.tsv> <msgf.tsv> <out.json>
 *   selection.tsv rows: scan<TAB>charge<TAB>context_peptide<TAB>golden_raw_score
 *                       (context_peptide keeps flanking, e.g. K.RSRRRRKR.A)
 *   msgf.tsv          : the MS-GF+ search TSV (DeNovoScore / MSGFScore / SpecEValue columns)
 *   model.param       : used only for the MODEL_NAME label; the scorer is built from the
 *                       equivalent bundled NewScorerFactory resource exactly as the search does.
 */
public class SpecProbDumper {

    // The iprg2013_F13 search used `-inst 1` = HighRes (Orbitrap/FTICR/Lumos), which loads the
    // HCD_HighRes_Tryp.param scoring model -- NOT QExactive (`-inst 3`, HCD_QExactive_Tryp.param).
    // Using the HighRes model is REQUIRED to reproduce the search's DeNovoScore / MSGFScore / SpecEValue.
    static final String MODEL_NAME = "HCD_HighRes_Tryp.param";

    // Reproduce search `-t 10ppm` and default `-ti 0,1`.
    static final int MIN_ISOTOPE_ERROR = 0;
    static final int MAX_ISOTOPE_ERROR = 1;
    static final Tolerance PRECURSOR_TOL = new Tolerance(10f, true); // 10 ppm

    static final int NUM_SCORE_DIST_SAMPLES = 3;

    static final class Sel {
        final int charge;
        final String contextPep; // with flanking, e.g. "K.RSRRRRKR.A"
        final int rawScore;
        Sel(int charge, String contextPep, int rawScore) {
            this.charge = charge; this.contextPep = contextPep; this.rawScore = rawScore;
        }
    }

    static final class TsvRow {
        final int deNovoScore;
        final int msgfScore;
        final double specEValue;
        TsvRow(int deNovoScore, int msgfScore, double specEValue) {
            this.deNovoScore = deNovoScore; this.msgfScore = msgfScore; this.specEValue = specEValue;
        }
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 6) {
            System.err.println("usage: SpecProbDumper <model.param> <spectra.mgf> <modfile> <selection.tsv> <msgf.tsv> <out.json>");
            System.exit(2);
        }
        File mgfFile = new File(args[1]);
        String modFilePath = args[2];
        File selFile = new File(args[3]);
        File tsvFile = new File(args[4]);
        File outFile = new File(args[5]);

        // 1. Scorer built EXACTLY as the default high-res search (ScoredSpectraMap): resource scorer.
        //    doNotUseError() is ONLY called when -turnOffEdgeScoring is set (it zeroes errorScalingFactor
        //    and makes supportEdgeScores() false -> node-only FastScorer). The default search does NOT set
        //    it, so we do NOT call it and edge scoring stays ON (DBScanScorer).
        NewRankScorer scorer = NewScorerFactory.get(
                ActivationMethod.HCD, InstrumentType.HIGH_RESOLUTION_LTQ, Enzyme.TRYPSIN, Protocol.STANDARD);
        if (!scorer.supportEdgeScores()) {
            System.err.println("WARNING: scorer does not support edge scores (edge scoring would be OFF)");
        }

        // 2. Amino-acid set from the mod file; register trypsin for cleavage constants.
        ParamManager pm = new ParamManager("SpecProbDumper", "1", "2026", "usage");
        AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSetFromModFile(modFilePath, pm);
        Enzyme enzyme = Enzyme.TRYPSIN;
        aaSet.registerEnzyme(enzyme);
        int neighborCredit = aaSet.getNeighboringAACleavageCredit();
        int neighborPenalty = aaSet.getNeighboringAACleavagePenalty();
        int peptideCredit = aaSet.getPeptideCleavageCredit();
        int peptidePenalty = aaSet.getPeptideCleavagePenalty();

        // 3. Selection (scan -> golden info).
        Map<Integer, Sel> selection = new LinkedHashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(selFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] c = line.split("\t", -1);
                selection.put(Integer.parseInt(c[0].trim()),
                        new Sel(Integer.parseInt(c[1].trim()), c[2], Integer.parseInt(c[3].trim())));
            }
        }

        // 4. TSV lookup: (scan, core-peptide) -> DeNovoScore / MSGFScore / SpecEValue.
        Map<String, TsvRow> tsvMap = loadTsv(tsvFile);

        SpectraAccessor accessor = new SpectraAccessor(mgfFile);
        Iterator<Spectrum> it = accessor.getSpecItr();
        List<String> emitted = new ArrayList<>();
        java.util.Set<Integer> seen = new java.util.HashSet<>();

        int nRawEqTsv = 0, nDenovoReconcile = 0, nSpecProbReconcile = 0, nTsvFound = 0;
        int sampleCount = 0;

        while (it.hasNext()) {
            Spectrum spec = it.next();
            int scan = spec.getScanNum();
            Sel sel = selection.get(scan);
            if (sel == null || seen.contains(scan)) continue;
            seen.add(scan);

            spec.setCharge(sel.charge);

            // Scored spectrum EXACTLY as the search.
            NewScoredSpectrum<NominalMass> ss = scorer.getScoredSpectrum(spec);

            // ---- peptide -> prefix mass arrays (leading zero; length n+1) ----
            String contextPep = sel.contextPep;
            char prevChar = contextPep.charAt(0);
            String core = contextPep.substring(2, contextPep.length() - 2);
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
                        throw new RuntimeException("no modified amino acid for " + res + "+" + delta + " in " + core);
                    aas.add(found);
                }
            }
            char lastResidue = 0;
            for (int k = 0; k < core.length(); k++)
                if (Character.isLetter(core.charAt(k))) lastResidue = core.charAt(k);

            int n = aas.size();
            double[] prm = new double[n + 1];
            int[] nominalPRM = new int[n + 1];
            for (int k = 0; k < n; k++) {
                prm[k + 1] = prm[k] + aas.get(k).getMass();
                nominalPRM[k + 1] = nominalPRM[k] + aas.get(k).getNominalMass();
            }
            int pepMassNominal = nominalPRM[n];

            // ---- precursor-derived nominal peptide mass + DBScanScorer (as ScoredSpectraMap) ----
            float peptideMass = spec.getPrecursorMass() - (float) Composition.H2O;
            int nominalPeptideMass = NominalMass.toNominalMass(peptideMass);
            float tolDaLeft = PRECURSOR_TOL.getToleranceAsDa(peptideMass);
            float tolDaRight = PRECURSOR_TOL.getToleranceAsDa(peptideMass);
            int dbNominal = nominalPeptideMass + Math.round(tolDaLeft - 0.4999f) - MIN_ISOTOPE_ERROR;
            DBScanScorer scoredSpec = new DBScanScorer(ss, dbNominal);

            // ---- raw score = cleavage + node+edge (== DBScanner match score) ----
            int nodeEdge = scoredSpec.getScore(prm, nominalPRM, 1, n + 1, numMods);
            int nTermScore = (prevChar == '-' || enzyme.isCleavable(prevChar)) ? neighborCredit : neighborPenalty;
            int cTermScore = enzyme.isCleavable(lastResidue) ? peptideCredit : peptidePenalty;
            int cleavage = nTermScore + cTermScore;
            int rawScore = nodeEdge + cleavage;

            // ---- generating-function GROUP over the isotope+tolerance mass-index range ----
            int minNominalPeptideMass = nominalPeptideMass - MAX_ISOTOPE_ERROR;
            int maxNominalPeptideMass = nominalPeptideMass - MIN_ISOTOPE_ERROR;
            int maxPeptideMassIndex = maxNominalPeptideMass + Math.round(tolDaLeft - 0.4999f);
            int minPeptideMassIndex = minNominalPeptideMass - Math.round(tolDaRight - 0.4999f);

            GeneratingFunctionGroup<NominalMass> gf = new GeneratingFunctionGroup<>();
            List<Integer> massIndices = new ArrayList<>();
            for (int idx = minPeptideMassIndex; idx <= maxPeptideMassIndex; idx++) {
                if (idx <= 0) continue;
                DeNovoGraph<NominalMass> graph = new FlexAminoAcidGraph(
                        aaSet, idx, enzyme, scoredSpec, false, false);
                GeneratingFunction<NominalMass> gfi = new GeneratingFunction<NominalMass>(graph)
                        .doNotBacktrack()
                        .doNotCalcNumber();
                gf.registerGF(graph.getPMNode(), gfi);
                massIndices.add(idx);
            }
            boolean gfComputed = gf.computeGeneratingFunction();

            int gfMaxScore = gfComputed ? gf.getMaxScore() : Integer.MIN_VALUE;
            int deNovoScore = gfComputed ? gfMaxScore - 1 : Integer.MIN_VALUE;
            double specProb = gfComputed ? gf.getSpectralProbability(rawScore) : 1.0;
            boolean inRange = pepMassNominal >= minPeptideMassIndex && pepMassNominal <= maxPeptideMassIndex;

            // ---- TSV oracle values ----
            String key = scan + "" + core;
            TsvRow tsv = tsvMap.get(key);
            boolean haveTsv = tsv != null;
            if (haveTsv) nTsvFound++;
            boolean rawEqTsv = haveTsv && rawScore == tsv.msgfScore;
            if (rawEqTsv) nRawEqTsv++;
            boolean denovoReconciles = haveTsv && deNovoScore == tsv.deNovoScore;
            if (denovoReconciles) nDenovoReconcile++;
            // spec prob reproduces SpecEValue if float-cast matches (MS-GF+ prints (float)specProb).
            boolean specProbReconciles = haveTsv
                    && Float.floatToIntBits((float) specProb) == Float.floatToIntBits((float) tsv.specEValue);
            if (specProbReconciles) nSpecProbReconcile++;
            double specProbRatio = (haveTsv && tsv.specEValue != 0) ? specProb / tsv.specEValue : 0;

            boolean emitDist = sampleCount < NUM_SCORE_DIST_SAMPLES;
            sampleCount++;

            // ---- emit ----
            StringBuilder b = new StringBuilder();
            b.append("  {\n");
            b.append("   \"scan\": ").append(scan).append(",\n");
            b.append("   \"charge\": ").append(sel.charge).append(",\n");
            b.append("   \"peptide\": ").append(jstr(core)).append(",\n");
            b.append("   \"num_mods\": ").append(numMods).append(",\n");
            b.append("   \"peptide_mass_nominal\": ").append(pepMassNominal).append(",\n");
            b.append("   \"precursor_nominal_mass\": ").append(nominalPeptideMass).append(",\n");
            b.append("   \"mass_index_range\": [").append(minPeptideMassIndex).append(",").append(maxPeptideMassIndex).append("],\n");
            b.append("   \"peptide_in_mass_range\": ").append(inRange).append(",\n");
            b.append("   \"node_edge_score\": ").append(nodeEdge).append(",\n");
            b.append("   \"cleavage_score\": ").append(cleavage).append(",\n");
            b.append("   \"raw_score\": ").append(rawScore).append(",\n");
            b.append("   \"gf_max_score\": ").append(gfMaxScore).append(",\n");
            b.append("   \"denovo_score\": ").append(deNovoScore).append(",\n");
            b.append("   \"spec_prob\": ").append(jd(specProb)).append(",\n");
            b.append("   \"tsv_denovo_score\": ").append(haveTsv ? Integer.toString(tsv.deNovoScore) : "null").append(",\n");
            b.append("   \"tsv_msgf_score\": ").append(haveTsv ? Integer.toString(tsv.msgfScore) : "null").append(",\n");
            b.append("   \"tsv_spec_evalue\": ").append(haveTsv ? jd(tsv.specEValue) : "null").append(",\n");
            b.append("   \"raw_score_equals_tsv_msgf\": ").append(rawEqTsv).append(",\n");
            b.append("   \"denovo_reconciles\": ").append(denovoReconciles).append(",\n");
            b.append("   \"spec_prob_reconciles\": ").append(specProbReconciles).append(",\n");
            b.append("   \"spec_prob_over_tsv_ratio\": ").append(jd(specProbRatio)).append(emitDist ? ",\n" : "\n");

            if (emitDist && gfComputed) {
                ScoreDist dist = gf.getScoreDist();
                int minS = dist.getMinScore();
                int maxS = dist.getMaxScore();
                StringBuilder probs = new StringBuilder("[");
                for (int s = minS; s < maxS; s++) {
                    if (s > minS) probs.append(",");
                    probs.append(jd(dist.getProbability(s)));
                }
                probs.append("]");
                b.append("   \"score_dist_sample\": {\"min_score\": ").append(minS)
                 .append(", \"max_score\": ").append(maxS)
                 .append(", \"probs\": ").append(probs).append("}\n");
            } else if (emitDist) {
                b.append("   \"score_dist_sample\": null\n");
            }
            b.append("  }");
            emitted.add(b.toString());
        }

        // ---- top-level JSON ----
        StringBuilder sb = new StringBuilder();
        sb.append("{\n");
        sb.append(" \"model\": ").append(jstr(MODEL_NAME)).append(",\n");
        sb.append(" \"note\": ").append(jstr(
            "MS-GF+ generating-function oracle (MS-GF:DeNovoScore + spectral probability / "
            + "MS-GF:SpecEValue) for the Rust port. GF built EXACTLY as DBScanner.computeSpecEValue: "
            + "scorer=NewScorerFactory.get(HCD,HIGH_RESOLUTION_LTQ,TRYPSIN,STANDARD) [-inst 1 HighRes] (edge scoring ON); "
            + "scoredSpec=DBScanScorer(ss, precursorNominalMass) (node+edge scoring ON); "
            + "graph=FlexAminoAcidGraph(aaSet, massIndex, TRYPSIN, scoredSpec, false, false) over the "
            + "isotope(-ti 0,1)+10ppm mass-index range, registered in a GeneratingFunctionGroup; "
            + "denovo_score=gf.getMaxScore()-1; spec_prob=gf.getSpectralProbability(raw_score); "
            + "raw_score=cleavage+DBScanScorer.getScore(prm,nominalPRM,1,len+1,numMods)=DBScanner match score. "
            + "setUpScoreThreshold is omitted (pruning only) so score_dist_sample is the FULL ScoreDist. "
            + "Generated by SpecProbDumper.java against MSGFPlus.jar; every number is from MS-GF+.")).append(",\n");
        sb.append(" \"gf_construction\": {")
          .append("\"scorer\": ").append(jstr("NewScorerFactory.get(HCD,HIGH_RESOLUTION_LTQ,TRYPSIN,STANDARD)")).append(", ")
          .append("\"graph_class\": ").append(jstr("FlexAminoAcidGraph")).append(", ")
          .append("\"scored_spectrum\": ").append(jstr("DBScanScorer (node+edge)")).append(", ")
          .append("\"edge_scoring\": true, ")
          .append("\"aa_set\": ").append(jstr("getAminoAcidSetFromModFile(iprg-2013_Mods.txt)")).append(", ")
          .append("\"enzyme\": ").append(jstr("Trypsin")).append(", ")
          .append("\"isotope_error\": [").append(MIN_ISOTOPE_ERROR).append(",").append(MAX_ISOTOPE_ERROR).append("], ")
          .append("\"precursor_tolerance\": ").append(jstr("10ppm")).append("},\n");
        sb.append(" \"cleavage_constants\": {")
          .append("\"neighboring_credit\": ").append(neighborCredit).append(", ")
          .append("\"neighboring_penalty\": ").append(neighborPenalty).append(", ")
          .append("\"peptide_credit\": ").append(peptideCredit).append(", ")
          .append("\"peptide_penalty\": ").append(peptidePenalty).append("},\n");
        sb.append(" \"self_check\": {")
          .append("\"n\": ").append(emitted.size()).append(", ")
          .append("\"tsv_found\": ").append(nTsvFound).append(", ")
          .append("\"raw_score_equals_tsv_msgf\": ").append(nRawEqTsv).append(", ")
          .append("\"denovo_reconciles\": ").append(nDenovoReconcile).append(", ")
          .append("\"spec_prob_reconciles\": ").append(nSpecProbReconcile).append("},\n");
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
        System.out.println("self-check: tsvFound " + nTsvFound + "/" + emitted.size()
            + ", raw==tsvMSGF " + nRawEqTsv + "/" + emitted.size()
            + ", denovo reconciles " + nDenovoReconcile + "/" + emitted.size()
            + ", specProb reconciles " + nSpecProbReconcile + "/" + emitted.size());
        if (emitted.size() != selection.size())
            System.out.println("WARNING: selected " + selection.size() + " scans but emitted " + emitted.size());
    }

    // Parse the MS-GF+ TSV into (scan, core-peptide) -> row. Handles flanking (X.PEP.Y)
    // and inline mod deltas containing '.' (e.g. M+15.995).
    static Map<String, TsvRow> loadTsv(File tsvFile) throws Exception {
        Map<String, TsvRow> map = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(tsvFile))) {
            String header = br.readLine();
            if (header == null) return map;
            String[] cols = header.split("\t", -1);
            int iScan = -1, iPep = -1, iDeNovo = -1, iMsgf = -1, iSpecE = -1;
            for (int j = 0; j < cols.length; j++) {
                String c = cols[j].trim();
                if (c.equals("ScanNum")) iScan = j;
                else if (c.equals("Peptide")) iPep = j;
                else if (c.equals("DeNovoScore")) iDeNovo = j;
                else if (c.equals("MSGFScore")) iMsgf = j;
                else if (c.equals("SpecEValue")) iSpecE = j;
            }
            if (iScan < 0 || iPep < 0 || iDeNovo < 0 || iMsgf < 0 || iSpecE < 0)
                throw new RuntimeException("TSV missing required columns");
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] c = line.split("\t", -1);
                int scan = Integer.parseInt(c[iScan].trim());
                String pep = c[iPep];
                String core = stripFlanking(pep);
                int deNovo = Integer.parseInt(c[iDeNovo].trim());
                int msgf = Integer.parseInt(c[iMsgf].trim());
                double specE = Double.parseDouble(c[iSpecE].trim());
                String key = scan + "" + core;
                // Keep the first (best) occurrence for a (scan, peptide) pair.
                if (!map.containsKey(key)) map.put(key, new TsvRow(deNovo, msgf, specE));
            }
        }
        return map;
    }

    // Strip single-char flanking context "X.PEPTIDE.Y" -> "PEPTIDE" (keeping inline mods).
    static String stripFlanking(String pep) {
        if (pep.length() >= 4 && pep.charAt(1) == '.' && pep.charAt(pep.length() - 2) == '.')
            return pep.substring(2, pep.length() - 2);
        return pep;
    }

    // Round-trip-safe double -> JSON number. Fails loudly on non-finite values.
    static String jd(double x) {
        if (!Double.isFinite(x)) throw new RuntimeException("non-finite double: " + x);
        return Double.toString(x);
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
