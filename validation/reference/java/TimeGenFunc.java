import edu.ucsd.msjava.msscorer.*;
import edu.ucsd.msjava.msgf.*;
import edu.ucsd.msjava.msutil.*;
import edu.ucsd.msjava.params.ParamManager;
import java.io.*;
import java.util.*;

public class TimeGenFunc {
    static final int MIN_ISO = 0, MAX_ISO = 1;
    static final Tolerance TOL = new Tolerance(10f, true); // 10 ppm

    public static void main(String[] a) throws Exception {
        String model = a[0], mgf = a[1], modFile = a[2], db = a[3];
        NewRankScorer scorer = NewScorerFactory.get(
                ActivationMethod.HCD, InstrumentType.HIGH_RESOLUTION_LTQ, Enzyme.TRYPSIN, Protocol.STANDARD);
        ParamManager pm = new ParamManager("t", "1", "2", "u");
        AminoAcidSet aaSet = AminoAcidSet.getAminoAcidSetFromModFile(modFile, pm);
        edu.ucsd.msjava.msdbsearch.DBScanner.setAminoAcidProbabilities(db, aaSet);
        aaSet.registerEnzyme(Enzyme.TRYPSIN);

        for (int pass = 0; pass < 5; pass++) {
            SpectraAccessor acc = new SpectraAccessor(new File(mgf));
            List<Spectrum> specs = new ArrayList<>();
            Iterator<Spectrum> it = acc.getSpecItr();
            while (it.hasNext()) { Spectrum s = it.next(); if (s.getCharge() > 0 && s.size() > 0) specs.add(s); }
            long t0 = System.nanoTime(); long chk = 0; int n = 0;
            for (Spectrum spec : specs) {
                float peptideMass = spec.getPrecursorMass() - (float) Composition.H2O;
                int nominalPeptideMass = NominalMass.toNominalMass(peptideMass);
                if (nominalPeptideMass < 200 || nominalPeptideMass > 6000) continue;
                NewScoredSpectrum<NominalMass> ss = scorer.getScoredSpectrum(spec);
                float tolDa = TOL.getToleranceAsDa(peptideMass);
                int dbNominal = nominalPeptideMass + Math.round(tolDa - 0.4999f) - MIN_ISO;
                DBScanScorer scoredSpec = new DBScanScorer(ss, dbNominal);
                int minIdx = (nominalPeptideMass - MAX_ISO) - Math.round(tolDa - 0.4999f);
                int maxIdx = (nominalPeptideMass - MIN_ISO) + Math.round(tolDa - 0.4999f);
                GeneratingFunctionGroup<NominalMass> gf = new GeneratingFunctionGroup<>();
                for (int idx = minIdx; idx <= maxIdx; idx++) {
                    if (idx <= 0) continue;
                    DeNovoGraph<NominalMass> graph = new FlexAminoAcidGraph(aaSet, idx, Enzyme.TRYPSIN, scoredSpec, false, false);
                    GeneratingFunction<NominalMass> gfi = new GeneratingFunction<NominalMass>(graph).doNotBacktrack().doNotCalcNumber();
                    gf.registerGF(graph.getPMNode(), gfi);
                }
                if (gf.computeGeneratingFunction()) { chk += gf.getMaxScore(); n++; }
            }
            long t1 = System.nanoTime(); double ms = (t1 - t0) / 1e6;
            System.out.printf("pass %d: %d spectra  %.1f ms  %.3f ms/spec  %.0f spec/s  (chk=%d)%n",
                    pass, n, ms, ms / n, n / (ms / 1000.0), chk);
        }
    }
}
