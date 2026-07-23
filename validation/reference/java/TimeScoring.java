import edu.ucsd.msjava.msscorer.*;
import edu.ucsd.msjava.msgf.NominalMass;
import edu.ucsd.msjava.msutil.*;
import java.io.*;
import java.util.*;

public class TimeScoring {
    public static void main(String[] a) throws Exception {
        NewRankScorer scorer = new NewRankScorer();
        scorer.readFromFile(new File(a[0]));
        for (int pass = 0; pass < 5; pass++) {
            // fresh spectra each pass (NewScoredSpectrum mutates the spectrum); read NOT timed
            SpectraAccessor acc = new SpectraAccessor(new File(a[1]));
            List<Spectrum> specs = new ArrayList<>();
            Iterator<Spectrum> it = acc.getSpecItr();
            while (it.hasNext()) { Spectrum s = it.next(); if (s.getCharge() > 0 && s.size() > 0) specs.add(s); }
            long t0 = System.nanoTime();
            long chk = 0;
            for (Spectrum s : specs) {
                NewScoredSpectrum<NominalMass> ss = new NewScoredSpectrum<>(s, scorer);
                int pepMass = NominalMass.toNominalMass(s.getPeptideMass());
                if (pepMass < 1) continue;
                FastScorer fs = new FastScorer(ss, pepMass);
                chk += pepMass;
            }
            long t1 = System.nanoTime();
            double ms = (t1 - t0) / 1e6;
            System.out.printf("pass %d: %d spectra  %.1f ms  %.3f ms/spec  %.0f spec/s  (chk=%d)%n",
                    pass, specs.size(), ms, ms / specs.size(), specs.size() / (ms / 1000.0), chk);
        }
    }
}
