/*
 * DumpFdrMap.java -- the PLAN2.md TD-2 "Gate 2" oracle for target-decoy FDR.
 *
 * The F13 search golden pins the q-value columns end-to-end, but it is a weak discriminator: it
 * yields only two distinct q-values (PLAN2.md section 4), so it cannot tell apart implementations
 * that differ in tie handling, in the `targetIndex > 0` guard, or in the map-lookup rule. This
 * probe drives edu.ucsd.msjava.fdr.TargetDecoyAnalysis directly on small synthetic score lists
 * chosen to separate exactly those behaviours, and freezes both the FDR/q-value map and the
 * lookups so the Rust port can be compared pair-by-pair.
 *
 * Two MS-GF+ entry points are exercised, both as the real search uses them (SpecEValue, so
 * isGreaterBetter = false, pit = 1):
 *
 *   TargetDecoyAnalysis.getFDRMap(target, decoy, false, 1f)   -> the threshold -> q-value map
 *   new TargetDecoyAnalysis(target, decoy).getPSMQValue(s)    -> the per-PSM lookup
 *
 * Probe scores include each map key's immediate float neighbours (Math.nextDown / Math.nextUp),
 * which is what pins the lookup rule: a floor lookup and a strictly-greater lookup disagree
 * exactly at the keys themselves.
 *
 * Floats are written as Float.toString() decimal strings, which round-trip exactly (and which
 * Rust's `str::parse::<f32>()` recovers bit-for-bit); infinities appear as "Infinity"/"-Infinity".
 *
 * Output: golden/fdr/fdrmap_cases.golden.json      Usage: java DumpFdrMap <out.json>
 */

import edu.ucsd.msjava.fdr.PSMSet;
import edu.ucsd.msjava.fdr.ScoredString;
import edu.ucsd.msjava.fdr.TargetDecoyAnalysis;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;
import java.util.TreeSet;

public class DumpFdrMap {

    /**
     * A PSMSet whose populations are handed over directly, so a TargetDecoyAnalysis can be built
     * from bare score lists. Both populations are filled with the same scores: MS-GF+ runs the PSM
     * and peptide levels through the identical getFDRMap code, so one list exercises both.
     */
    static class ProbeSet extends PSMSet {
        ProbeSet(float[] scores) {
            psmList = new ArrayList<ScoredString>();
            peptideScoreTable = new HashMap<String, Float>();
            for (int i = 0; i < scores.length; i++) {
                psmList.add(new ScoredString("p" + i, scores[i]));
                peptideScoreTable.put("p" + i, scores[i]);
            }
        }

        @Override
        public boolean isGreaterBetter() {
            return false; // SpecEValue: smaller is better
        }

        @Override
        public void read() {
        }
    }

    static class Case {
        final String name;
        final String note;
        final float[] targets;
        final float[] decoys;

        Case(String name, String note, float[] targets, float[] decoys) {
            this.name = name;
            this.note = note;
            this.targets = targets;
            this.decoys = decoys;
        }
    }

    /** Pseudo-random SpecEValue-like scores in [1e-12, 1e-1], fixed seed for reproducibility. */
    static float[] randomScores(long seed, int n) {
        Random rnd = new Random(seed);
        float[] out = new float[n];
        for (int i = 0; i < n; i++)
            out[i] = (float) Math.pow(10, -1 - 11 * rnd.nextDouble());
        return out;
    }

    /** Random scores drawn from a *better* range, so targets and decoys only partly overlap. */
    static float[] randomScoresBiased(long seed, int n, double lo, double hi) {
        Random rnd = new Random(seed);
        float[] out = new float[n];
        for (int i = 0; i < n; i++)
            out[i] = (float) Math.pow(10, lo + (hi - lo) * rnd.nextDouble());
        return out;
    }

    static Case[] cases() {
        return new Case[]{
                // Deliberately synthetic. This case used to carry the literal head of the F13
                // golden, but those SpecEValues are frozen output of *running* MS-GF+ on UC test
                // data, which LICENSING.md excludes from the repository. The shape is what has
                // discriminating power -- one target better than every decoy, then interleaved --
                // and it is reproduced here with numbers of our own.
                new Case("head_one_target_then_interleaved",
                        "one target better than every decoy, then the two lists interleave",
                        new float[]{9e-10f, 5e-9f, 6.2e-9f, 7.8e-9f, 1.4e-8f},
                        new float[]{2.4e-9f, 2.6e-9f, 2.7e-9f, 8.3e-9f, 8.7e-9f}),

                new Case("tie_within_decoys",
                        "two decoys share a score: does the sweep step by one decoy or by the whole run?",
                        new float[]{1e-9f, 1e-8f, 1e-7f, 1e-6f, 1e-5f},
                        new float[]{1e-7f, 1e-7f}),

                new Case("tie_run_of_three",
                        "a longer run of equal decoy scores, with targets on both sides of it",
                        new float[]{1e-10f, 2e-10f, 3e-10f, 1e-6f, 2e-6f, 3e-6f},
                        new float[]{5e-9f, 5e-9f, 5e-9f, 1e-5f}),

                new Case("target_equals_decoy_key",
                        "targets sitting exactly on decoy keys: separates a floor lookup from a strictly-greater one",
                        new float[]{1e-9f, 1e-8f, 1e-7f, 1e-6f, 1e-5f, 1e-4f, 1e-3f},
                        new float[]{1e-7f, 1e-5f, 1e-3f}),

                new Case("guard_no_target_better",
                        "no target beats the first decoys: MS-GF+ writes no entry and keeps sweeping (targetIndex > 0 guard)",
                        new float[]{5e-9f, 6e-9f, 7e-9f},
                        new float[]{1e-9f, 1e-8f}),

                new Case("all_decoys_better",
                        "every decoy beats every target",
                        new float[]{1e-3f},
                        new float[]{1e-9f, 1e-8f}),

                new Case("empty_decoys",
                        "no decoys at all -- the sweep body never runs",
                        new float[]{1e-9f, 1e-8f, 1e-7f},
                        new float[]{}),

                new Case("empty_targets",
                        "no targets at all",
                        new float[]{},
                        new float[]{1e-9f, 1e-8f}),

                new Case("single_each",
                        "one target, one decoy",
                        new float[]{1e-9f},
                        new float[]{1e-8f}),

                new Case("early_break",
                        "FDR reaches 1 at an early threshold although later thresholds would be lower: the sweep stops",
                        new float[]{1e-9f, 1e-4f, 2e-4f, 3e-4f, 4e-4f, 5e-4f, 6e-4f},
                        new float[]{5e-9f, 6e-9f, 7e-9f, 1e-3f}),

                new Case("monotone_spike",
                        "raw FDR is non-monotone across thresholds; the running-minimum pass must flatten it",
                        new float[]{1e-9f, 2e-9f, 3e-9f, 4e-9f, 5e-9f, 6e-9f, 7e-9f, 8e-9f, 9e-9f, 1e-8f,
                                1.1e-8f, 1.2e-8f, 1.3e-8f, 1.4e-8f, 1.5e-8f, 1.6e-8f, 1.7e-8f, 1.8e-8f,
                                1.9e-8f, 2e-8f},
                        new float[]{5e-10f, 1.9e-8f}),

                new Case("random_overlapping", "50 vs 50 drawn from the same range",
                        randomScores(20260724L, 50), randomScores(981231L, 50)),

                new Case("random_separated", "targets biased better than decoys, the realistic shape",
                        randomScoresBiased(4242L, 60, -12, -6), randomScoresBiased(7777L, 40, -8, -1)),

                new Case("random_sparse", "few decoys against many targets",
                        randomScores(13L, 40), randomScoresBiased(99L, 5, -10, -2)),
        };
    }

    /** Every score worth asking about: the inputs, each map key, and each key's float neighbours. */
    static float[] probeScores(Case c, TreeMap<Float, Float> map) {
        TreeSet<Float> probes = new TreeSet<Float>();
        for (float t : c.targets) probes.add(t);
        for (float d : c.decoys) probes.add(d);
        for (Float k : map.keySet()) {
            if (k.isInfinite()) continue;
            probes.add(k);
            probes.add(Math.nextDown(k));
            probes.add(Math.nextUp(k));
        }
        probes.add(0f);
        probes.add(Float.MIN_VALUE);
        probes.add(1e-30f);
        probes.add(1f);
        probes.add(1e6f);
        float[] out = new float[probes.size()];
        int i = 0;
        for (Float p : probes) out[i++] = p;
        return out;
    }

    static String q(String s) {
        return "\"" + s + "\"";
    }

    static String f(float v) {
        return q(Float.toString(v));
    }

    static String floats(float[] v) {
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < v.length; i++) {
            if (i > 0) sb.append(", ");
            sb.append(f(v[i]));
        }
        return sb.append("]").toString();
    }

    public static void main(String[] argv) throws Exception {
        if (argv.length != 1) {
            System.err.println("usage: java DumpFdrMap <out.json>");
            System.exit(1);
        }
        PrintStream out = new PrintStream(argv[0], "UTF-8");
        out.println("{");
        out.println("  \"generator\": \"validation/reference/java/DumpFdrMap.java against MSGFPlus.jar\",");
        out.println("  \"source\": \"edu.ucsd.msjava.fdr.TargetDecoyAnalysis (getFDRMap + getPSMQValue)\",");
        out.println("  \"score_semantics\": \"SpecEValue, smaller is better (isGreaterBetter=false), pit=1\",");
        out.println("  \"encoding\": \"every float is a Float.toString() decimal string; infinities are Infinity/-Infinity\",");
        out.println("  \"cases\": [");

        Case[] cases = cases();
        for (int ci = 0; ci < cases.length; ci++) {
            Case c = cases[ci];

            ArrayList<Float> t = new ArrayList<Float>();
            for (float v : c.targets) t.add(v);
            ArrayList<Float> d = new ArrayList<Float>();
            for (float v : c.decoys) d.add(v);
            // getFDRMap sorts its arguments in place; hand it copies so the dumped inputs stay as given.
            TreeMap<Float, Float> map = TargetDecoyAnalysis.getFDRMap(
                    new ArrayList<Float>(t), new ArrayList<Float>(d), false, 1f);

            TargetDecoyAnalysis tda = new TargetDecoyAnalysis(new ProbeSet(c.targets), new ProbeSet(c.decoys));

            out.println("    {");
            out.println("      \"name\": " + q(c.name) + ",");
            out.println("      \"note\": " + q(c.note) + ",");
            out.println("      \"targets\": " + floats(c.targets) + ",");
            out.println("      \"decoys\": " + floats(c.decoys) + ",");

            out.print("      \"map\": [");
            boolean first = true;
            for (Map.Entry<Float, Float> e : map.entrySet()) {
                if (!first) out.print(", ");
                first = false;
                out.print("{\"key\": " + f(e.getKey()) + ", \"q\": " + f(e.getValue()) + "}");
            }
            out.println("],");

            float[] probes = probeScores(c, map);
            out.print("      \"lookups\": [");
            for (int i = 0; i < probes.length; i++) {
                if (i > 0) out.print(", ");
                String value;
                try {
                    value = f(tda.getPSMQValue(probes[i]));
                } catch (NullPointerException npe) {
                    // getPSMQValue does an unguarded higherEntry(...).getValue(); record where that
                    // has no answer rather than pretending it returns something.
                    value = "null";
                }
                out.print("{\"score\": " + f(probes[i]) + ", \"q\": " + value + "}");
            }
            out.println("]");
            out.println("    }" + (ci + 1 < cases.length ? "," : ""));
        }
        out.println("  ]");
        out.println("}");
        out.close();
        System.err.println("DumpFdrMap: wrote " + cases.length + " cases to " + argv[0]);
    }
}
