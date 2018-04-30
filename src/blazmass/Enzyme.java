/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package blazmass;

/**
 *
 * @author rpark
 */
public class Enzyme {

    public static final int SIZE = 256;
    private static boolean[] enzymeArr = new boolean[SIZE];
    //private int miscleavage=2;    
    //private int maxInternalMiscleavage=-1;    

    public Enzyme() {
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();

        sb.append("Enzymes:");
        for (int i = 0; i < enzymeArr.length; ++i) {
            if (enzymeArr[i] == true) {
                sb.append((char) i);
            }
        }
        sb.append(" ");

        return sb.toString();
    }


    public static int checkCleavage(String seq, int start, int end) {

        int status = 0;

        if (start <= 0) {
            status++;
        } else {
            if (enzymeArr[seq.charAt(start - 1)]) {
                status++;
            }
        }

        if (enzymeArr[seq.charAt(end)]) {
            status++;
        }

        return status;
    }

    public void addEnzyme(char ch) {
        enzymeArr[ch] = true;
    }

    public static boolean[] getEnzymeArr() {
        return enzymeArr;
    }

    public static void setEnzymeArr(boolean[] enzymeArr) {
        Enzyme.enzymeArr = enzymeArr;
    }


    public static boolean isEnzyme(char ch) {
        return enzymeArr[ch];
    }

    /*
     public int getMaxInternalMiscleavage() {
     return maxInternalMiscleavage;
     }

     public void setMaxInternalMiscleavage(int maxInternalMiscleavage) {
     this.maxInternalMiscleavage = maxInternalMiscleavage;
     }

     public int getMiscleavage() {
     return miscleavage;
     }

     public void setMiscleavage(int miscleavage) {
     this.miscleavage = miscleavage;
     }
     */

    public static int checkCleavage(String seq, int start, int end, String enzymeNoCutAA) {
        int status = 0;

        if (start <= 0) {

            status++;

        } else {

            if (enzymeArr[seq.charAt(start - 1)]) {

                status++;

            }
        }

        if (enzymeArr[seq.charAt(end)]) {
            status++;
        }

        // new modification by Salva (July 2015)
        if (seq.length() > end + 1 && enzymeNoCutAA != null
                && enzymeNoCutAA.contains(seq.substring(end + 1, end + 2))) {

            status--;
        }

        // MODIFICATION BY SALVA
        if (status != 2 && end == seq.length() - 1) {
            status++;
        }

        return status;
    }

}
