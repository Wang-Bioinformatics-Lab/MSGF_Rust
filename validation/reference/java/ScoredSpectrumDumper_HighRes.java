import edu.ucsd.msjava.msscorer.NewRankScorer;
import edu.ucsd.msjava.msscorer.NewScoredSpectrum;
import edu.ucsd.msjava.msgf.NominalMass;
import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.SpectraAccessor;
import edu.ucsd.msjava.msutil.Spectrum;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.io.FileOutputStream;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * ScoredSpectrumDumper_HighRes — HighRes twin of ScoredSpectrumDumper.
 *
 * IDENTICAL to ScoredSpectrumDumper.java in every respect except the MODEL_NAME
 * label constant, which is set to "HCD_HighRes_Tryp.param" so the emitted JSON
 * "model" field honestly names the HighRes model the F13 search actually used
 * (-inst 1 = HighRes). The actual scoring model is loaded from args[0] in both
 * classes; this copy exists only to keep the emitted label correct.
 *
 * For a selected set of high-res MS/MS spectra it runs the reference MS-GF+
 * NewRankScorer + NewScoredSpectrum preprocessing/scoring and dumps, per
 * spectrum: the preprocessed (filtered/deconvolved/ranked) peak list and the
 * full per-nominal-mass prefix/suffix node-score vectors. Every number comes
 * from actually running MS-GF+.
 *
 * Args: <model.param> <spectra.mgf> <selection.tsv> <out.json>
 *   selection.tsv rows: scan<TAB>charge<TAB>golden_peptide<TAB>golden_raw_score
 *
 * Emits the JSON contract expected by the Rust diff harness with a small
 * structural indent; each peak row and each score vector is kept compact on a
 * single line so the artifact stays inspectable.
 */
public class ScoredSpectrumDumper_HighRes {

    static final String MODEL_NAME = "HCD_HighRes_Tryp.param";

    // Selection entry from the TSV.
    static final class Sel {
        final int charge;
        final String peptide;
        final int rawScore;
        Sel(int charge, String peptide, int rawScore) {
            this.charge = charge; this.peptide = peptide; this.rawScore = rawScore;
        }
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 4) {
            System.err.println("usage: ScoredSpectrumDumper_HighRes <model.param> <spectra.mgf> <selection.tsv> <out.json>");
            System.exit(2);
        }
        File modelFile = new File(args[0]);
        File mgfFile = new File(args[1]);
        File selFile = new File(args[2]);
        File outFile = new File(args[3]);

        // 1. Load the scoring model.
        NewRankScorer scorer = new NewRankScorer();
        scorer.readFromFile(modelFile);

        // 2. Load the selection (scan -> golden info), preserving insertion order.
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

        // Reflection handle for the private preprocessed peak list.
        Field specField = NewScoredSpectrum.class.getDeclaredField("spec");
        specField.setAccessible(true);

        // 3. Iterate spectra in file order; dump the selected ones.
        SpectraAccessor accessor = new SpectraAccessor(mgfFile);
        Iterator<Spectrum> it = accessor.getSpecItr();
        List<String> emitted = new ArrayList<>();
        java.util.Set<Integer> seen = new java.util.HashSet<>();

        while (it.hasNext()) {
            Spectrum spec = it.next();
            int scan = spec.getScanNum();
            Sel sel = selection.get(scan);
            if (sel == null || seen.contains(scan)) continue;
            seen.add(scan);

            // Build the scored spectrum (preprocesses the spectrum: precursor
            // filtering, deconvolution, rank assignment).
            NewScoredSpectrum<NominalMass> ss = new NewScoredSpectrum<>(spec, scorer);

            int pepMassNominal = NominalMass.toNominalMass(spec.getPeptideMass());

            // Preprocessed peak list via reflection (file order).
            Spectrum ps = (Spectrum) specField.get(ss);

            emitted.add(buildSpectrumJson(scan, spec, ss, sel, pepMassNominal, ps));
        }

        // 4. Assemble the top-level JSON with a small structural indent.
        StringBuilder sb = new StringBuilder();
        sb.append("{\n");
        sb.append(" \"model\": ").append(jstr(MODEL_NAME)).append(",\n");
        sb.append(" \"note\": ").append(jstr(
            "MS-GF+ NewScoredSpectrum dump for the Rust RawScore oracle: preprocessed peaks "
            + "(post precursor-filter/deconvolution/rank) and per-nominal-mass prefix/suffix "
            + "node scores (peptide-independent). prefix_score[nm]=getNodeScore(NominalMass(nm),true), "
            + "suffix_score[nm]=getNodeScore(...,false); index 0 = 0.0. Generated by "
            + "ScoredSpectrumDumper_HighRes.java against MSGFPlus.jar.")).append(",\n");
        sb.append(" \"n_spectra\": ").append(emitted.size()).append(",\n");
        sb.append(" \"spectra\": [\n");
        for (int i = 0; i < emitted.size(); i++) {
            sb.append(emitted.get(i));
            sb.append(i + 1 < emitted.size() ? ",\n" : "\n");
        }
        sb.append(" ]\n");
        sb.append("}\n");

        try (PrintStream out = new PrintStream(new FileOutputStream(outFile), false, "UTF-8")) {
            out.print(sb);
        }
        System.out.println("dumped " + emitted.size() + " spectra -> " + outFile.getPath());
        if (emitted.size() != selection.size()) {
            System.out.println("WARNING: selected " + selection.size()
                + " scans but only " + emitted.size() + " found in the MGF");
        }
    }

    static String buildSpectrumJson(int scan, Spectrum spec, NewScoredSpectrum<NominalMass> ss,
                                    Sel sel, int pepMassNominal, Spectrum ps) {
        StringBuilder b = new StringBuilder();
        b.append("  {\n");
        b.append("   \"scan\": ").append(scan).append(",\n");
        b.append("   \"charge\": ").append(spec.getCharge()).append(",\n");
        b.append("   \"precursor_mass\": ").append(jf(spec.getPrecursorMass())).append(",\n");
        b.append("   \"peptide_mass_nominal\": ").append(pepMassNominal).append(",\n");
        b.append("   \"main_ion_is_prefix\": ").append(ss.getMainIonDirection()).append(",\n");
        b.append("   \"prob_peak\": ").append(jf(ss.getProbPeak())).append(",\n");

        // peaks: [[mz, intensity, rank], ...] in file order, single line.
        b.append("   \"peaks\": [");
        for (int i = 0; i < ps.size(); i++) {
            Peak pk = ps.get(i);
            if (i > 0) b.append(",");
            b.append("[").append(jf(pk.getMz())).append(",")
             .append(jf(pk.getIntensity())).append(",")
             .append(pk.getRank()).append("]");
        }
        b.append("],\n");

        // prefix/suffix node scores, index 0..pepMassNominal-1; index 0 = 0.0.
        float[] prefix = new float[pepMassNominal];
        float[] suffix = new float[pepMassNominal];
        for (int nm = 1; nm < pepMassNominal; nm++) {
            prefix[nm] = ss.getNodeScore(new NominalMass(nm), true);
            suffix[nm] = ss.getNodeScore(new NominalMass(nm), false);
        }
        b.append("   \"prefix_score\": ").append(jfArray(prefix)).append(",\n");
        b.append("   \"suffix_score\": ").append(jfArray(suffix)).append(",\n");

        b.append("   \"golden_peptide\": ").append(jstr(sel.peptide)).append(",\n");
        b.append("   \"golden_raw_score\": ").append(sel.rawScore).append("\n");
        b.append("  }");
        return b.toString();
    }

    // Compact JSON array of floats on a single line.
    static String jfArray(float[] a) {
        StringBuilder b = new StringBuilder("[");
        for (int i = 0; i < a.length; i++) {
            if (i > 0) b.append(",");
            b.append(jf(a[i]));
        }
        return b.append("]").toString();
    }

    // Round-trip-safe float -> JSON number. Fails loudly on non-finite values so
    // we never silently emit invalid JSON or fabricate a value.
    static String jf(float x) {
        if (!Float.isFinite(x)) {
            throw new RuntimeException("non-finite float encountered: " + x);
        }
        return Float.toString(x);
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
