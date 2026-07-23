import edu.ucsd.msjava.msscorer.*;
import edu.ucsd.msjava.msutil.IonType;
import java.io.*;
import java.util.*;
import java.lang.reflect.*;

public class ScoreDumper {
    @SuppressWarnings("unchecked")
    public static void main(String[] a) throws Exception {
        NewRankScorer s = new NewRankScorer();
        s.readFromFile(new File(a[0]));
        PrintStream out = new PrintStream(new FileOutputStream(a[1]));
        Method gps = NewRankScorer.class.getDeclaredMethod("getPartitionSet");
        gps.setAccessible(true);
        TreeSet<Partition> parts = (TreeSet<Partition>) gps.invoke(s);
        int[] ranks = {1, 2, 3, 5, 10, 50, 100, 149, 150, 151};
        out.println("pi\tcharge\tseg\tparentMass\tion\tion_charge\trank\tscore");
        int pi = 0;
        for (Partition p : parts) {
            IonType[] ions = s.getIonTypes(p.getCharge(), p.getParentMass(), p.getSegNum());
            for (IonType ion : ions) {
                for (int r : ranks)
                    out.println(pi + "\t" + p.getCharge() + "\t" + p.getSegNum() + "\t" + p.getParentMass()
                            + "\t" + ion.getName() + "\t" + ion.getCharge() + "\t" + r + "\t" + s.getNodeScore(p, ion, r));
                out.println(pi + "\t" + p.getCharge() + "\t" + p.getSegNum() + "\t" + p.getParentMass()
                        + "\t" + ion.getName() + "\t" + ion.getCharge() + "\tMISSING\t" + s.getMissingIonScore(p, ion));
            }
            pi++;
        }
        out.close();
        System.out.println("dumped " + pi + " partitions -> " + a[1]);
    }
}
