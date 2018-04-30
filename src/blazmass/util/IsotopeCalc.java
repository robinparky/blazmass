/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package blazmass.util;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 *
 * @author rpark2
 */
public class IsotopeCalc {
    
    //private static final int NUMCHARS = 256;
    private static final int NUMCHARS = 128;
    private static final double [] avgAaMasses = new double[NUMCHARS];
    private static final double [] monoAaMasses = new double[NUMCHARS];
    private static final int [] monoIntAaMasses = new int[NUMCHARS];
    private static final double [][] isotopicDistribution = new double[20][20];
    public static final int[][] intensePeaks = new int[20][];
    public static final int[] mostIntensePeaks = new int[20];

  
    /**
     * get the isotopic peaks that has intensity greater than or equals to the threshold 
     * relative to the most intense peak.  e.g., as mass == 500, the isotopic  
     * distribution is M+0 = 1, M+1 = 0.27, M+2 = 0.0512 and M+3 = 0.007109.
     * if the threshold is 0.05, a list of 0, 1, 2 will be returned, 
     * and if the threshold is 0.1, then a list of 0, 1 will be returned
     */
    public static List<Integer> mass2IsotopicPeaks(double mass, double threshold) {
        ArrayList<Integer> intensPeaks = new ArrayList<Integer>();
        int massIndex = ((int)(mass/500 + 0.5)) - 1;
        massIndex = massIndex < 0? 0 : massIndex;
         //       System.out.println("\nmass: " + mass + "\tmassIndex: " + massIndex);
        for(int i = 0; i < 20; i++) {
            if(isotopicDistribution[massIndex][i] >= threshold) {
                intensPeaks.add(new Integer(i));
                //System.out.print(i + "\t" + "intens: " + isotopicDistribution[massIndex][i] + "\n");
            }
        }
        
        return intensPeaks;
    }
    public static int [] mass2IntensePeaks(double mass) {
 
        int massIndex = ((int)(mass/500 + 0.5)) - 1;
        massIndex = massIndex < 0? 0 : massIndex;
        return intensePeaks[massIndex];
    }
    public static void loadIntensePeaks() {
        for(int i = 0; i < 20; i++) {
            int numPeaks = 0;
            for(int j = 0; j < 20; j++) {
                
                if(isotopicDistribution[i][j] >= .25) {
                    numPeaks++;
                } else {
                }
            }
            intensePeaks[i]  = new int[numPeaks]; 
            for(int j = 0; j < 20; j++) {
                if(isotopicDistribution[i][j] >= .25) {
                    intensePeaks[i][--numPeaks] = j;
                }
            }
             
        }
    }
    public static void loadMostIntensePeaks() {
        for(int i = 0; i < 20; i++) {
            for(int j = 0; j < 20; j++) {
                if(isotopicDistribution[i][j] == 1.0f) {
                    mostIntensePeaks[i] = j;
                    break;
                }
            }
             
        }
    }
    public static int getMostIntensePeak(double mass) {

        int massIndex = ((int)(mass/500 + 0.5)) - 1;
        massIndex = massIndex < 0? 0 : massIndex;
        return mostIntensePeaks[massIndex];
    }
    /*
    public double getPrecursorMass(byte [] seq, int start, int end) {
        double mass = MassSpecConstants.MASSH2O + MassSpecConstants.MASSPROTON;
        for(int i = start; i <= end; i++) {
            mass += getPrecursorMass(seq[i]);
        }
        return mass;
    }
    public double getPrecursorMass(String sequence) {
        int len = sequence.length();
        double mass = MassSpecConstants.MASSH2O + MassSpecConstants.MASSPROTON;
        for(int i = 0; i < len; i++) {
            mass += getPrecursorMass(sequence.charAt(i));
        }
        return mass;
    }*/
   
    private double [] getDBinWidthMasses(double masses[], double dBinWidth) {
        int length = masses.length;
        double [] dBinWidthMasses = new double[length];
        for(int i = 0; i < length; i++) {
            if(masses[i] !=0) {
                dBinWidthMasses[i] = masses[i]*dBinWidth;
            }
        }
        return dBinWidthMasses;
    }
    
    /*
    public double[] getPrecMasses(double dBinWidth){
        return getDBinWidthMasses(getPrecMasses(), dBinWidth);
    } */
 
    public double [] getTheoreticMasses(byte [] seq) {
        int numIons = (seq.length - 1 ) * 2;
        double [] masses = new double [numIons];

        return masses;
    }
    
    /**
     *
     *@param isotope - isotopic type, either 'mono' or 'avg'
     *@param ms - modifications
     */
    /*
    private double [] getMasses(String isotope, Iterator<Modification> ms) {
        double [] masses = new double[NUMCHARS];
        String monoIsotope = MassSpecConstants.MONOISOTOPE;
        String avgIsotope = MassSpecConstants.AVGISOTOPE;

        double [] origMasses = null;
        if (monoIsotope.equals(isotope)) {
            origMasses = getMonoMasses(); 
        } else if(avgIsotope.equals(isotope)){
            origMasses = getAvgMasses();
        } else {
            throw new InvalidArgumentException("Unknown isotope type, should be mono or avg");
        }

        for(int i = 0; i < NUMCHARS; i++) {
            masses[i] = origMasses[i];
        }

        while (ms.hasNext()) {
            Modification m = ms.next();
            masses[m.getResidue()] += m.getMassShift();
        }
        return masses; 
    }*/

    public static double getAvgMass(byte aa) {
        return avgAaMasses[aa];
    }
    public static double getAvgMass(char aa) {
        return avgAaMasses[aa];
    }

    public static double getMonoMass(byte aa) {
        return monoAaMasses[aa];
    }
    public static int getIntMonoMass(byte aa) {
        return monoIntAaMasses[aa];
    }
    public static int getIntMonoMass(char aa) {
        return monoIntAaMasses[aa];
    }
    public static double getMonoMass(char aa) {
        return monoAaMasses[aa];
    }

    public static double [] getAvgMasses() {
        return avgAaMasses;
    }
    public static double [] getMonoMasses() {
        return monoAaMasses;
    }
    public static int [] getIntMonoMasses() {
        return monoIntAaMasses;
    }
    
    /*
    public static double getMonoMass(String sequence, int chargeState) {
        int len = sequence.length();
        double mass = MassSpecConstants.MASSH2O + MassSpecConstants.MASSPROTON * (chargeState-1);
        for(int i = 0; i < len; i++) {
            mass += getMonoMass(sequence.charAt(i));
            if(sequence.charAt(i) == 'C') {
                mass += 57.02146;
            }
                
        }

        return mass;
    }*/

    /*
    public static double getAvgMass(String sequence, int chargeState) {
        int len = sequence.length();
        double mass = MassSpecConstants.MASSH2O + MassSpecConstants.MASSPROTON * (chargeState-1);
        for(int i = 0; i < len; i++) {
            mass += getAvgMass(sequence.charAt(i));
            if(sequence.charAt(i) == 'C') {
                mass += 57.02146;
            }
        }
        return mass;
    }*/

      static {
        // from averagine
        isotopicDistribution[0][0] = 1.000000;
        isotopicDistribution[0][1] = 0.276950;
        isotopicDistribution[0][2] = 0.051217;
        isotopicDistribution[0][3] = 0.007109;
        isotopicDistribution[0][4] = 0.000000;
        isotopicDistribution[0][5] = 0.000000;
        isotopicDistribution[0][6] = 0.000000;
        isotopicDistribution[0][7] = 0.000000;
        isotopicDistribution[0][8] = 0.000000;
        isotopicDistribution[0][9] = 0.000000;
        isotopicDistribution[0][10] = 0.000000;
        isotopicDistribution[0][11] = 0.000000;
        isotopicDistribution[0][12] = 0.000000;
        isotopicDistribution[0][13] = 0.000000;
        isotopicDistribution[0][14] = 0.000000;
        isotopicDistribution[0][15] = 0.000000;
        isotopicDistribution[0][16] = 0.000000;
        isotopicDistribution[0][17] = 0.000000;
        isotopicDistribution[0][18] = 0.000000;
        isotopicDistribution[0][19] = 0.000000;
        isotopicDistribution[1][0] = 1.000000;
        isotopicDistribution[1][1] = 0.555776;
        isotopicDistribution[1][2] = 0.178136;
        isotopicDistribution[1][3] = 0.041790;
        isotopicDistribution[1][4] = 0.007894;
        isotopicDistribution[1][5] = 0.001263;
        isotopicDistribution[1][6] = 0.000000;
        isotopicDistribution[1][7] = 0.000000;
        isotopicDistribution[1][8] = 0.000000;
        isotopicDistribution[1][9] = 0.000000;
        isotopicDistribution[1][10] = 0.000000;
        isotopicDistribution[1][11] = 0.000000;
        isotopicDistribution[1][12] = 0.000000;
        isotopicDistribution[1][13] = 0.000000;
        isotopicDistribution[1][14] = 0.000000;
        isotopicDistribution[1][15] = 0.000000;
        isotopicDistribution[1][16] = 0.000000;
        isotopicDistribution[1][17] = 0.000000;
        isotopicDistribution[1][18] = 0.000000;
        isotopicDistribution[1][19] = 0.000000;
        isotopicDistribution[2][0] = 1.000000;
        isotopicDistribution[2][1] = 0.845212;
        isotopicDistribution[2][2] = 0.437964;
        isotopicDistribution[2][3] = 0.168570;
        isotopicDistribution[2][4] = 0.052286;
        isotopicDistribution[2][5] = 0.013663;
        isotopicDistribution[2][6] = 0.003097;
        isotopicDistribution[2][7] = 0.000000;
        isotopicDistribution[2][8] = 0.000000;
        isotopicDistribution[2][9] = 0.000000;
        isotopicDistribution[2][10] = 0.000000;
        isotopicDistribution[2][11] = 0.000000;
        isotopicDistribution[2][12] = 0.000000;
        isotopicDistribution[2][13] = 0.000000;
        isotopicDistribution[2][14] = 0.000000;
        isotopicDistribution[2][15] = 0.000000;
        isotopicDistribution[2][16] = 0.000000;
        isotopicDistribution[2][17] = 0.000000;
        isotopicDistribution[2][18] = 0.000000;
        isotopicDistribution[2][19] = 0.000000;
        isotopicDistribution[3][0] = 0.891137;
        isotopicDistribution[3][1] = 1.000000;
        isotopicDistribution[3][2] = 0.644526;
        isotopicDistribution[3][3] = 0.303221;
        isotopicDistribution[3][4] = 0.114261;
        isotopicDistribution[3][5] = 0.036227;
        isotopicDistribution[3][6] = 0.009967;
        isotopicDistribution[3][7] = 0.002431;
        isotopicDistribution[3][8] = 0.000000;
        isotopicDistribution[3][9] = 0.000000;
        isotopicDistribution[3][10] = 0.000000;
        isotopicDistribution[3][11] = 0.000000;
        isotopicDistribution[3][12] = 0.000000;
        isotopicDistribution[3][13] = 0.000000;
        isotopicDistribution[3][14] = 0.000000;
        isotopicDistribution[3][15] = 0.000000;
        isotopicDistribution[3][16] = 0.000000;
        isotopicDistribution[3][17] = 0.000000;
        isotopicDistribution[3][18] = 0.000000;
        isotopicDistribution[3][19] = 0.000000;
        isotopicDistribution[4][0] = 0.712904;
        isotopicDistribution[4][1] = 1.000000;
        isotopicDistribution[4][2] = 0.775825;
        isotopicDistribution[4][3] = 0.432162;
        isotopicDistribution[4][4] = 0.191232;
        isotopicDistribution[4][5] = 0.070896;
        isotopicDistribution[4][6] = 0.022756;
        isotopicDistribution[4][7] = 0.006467;
        isotopicDistribution[4][8] = 0.001654;
        isotopicDistribution[4][9] = 0.000000;
        isotopicDistribution[4][10] = 0.000000;
        isotopicDistribution[4][11] = 0.000000;
        isotopicDistribution[4][12] = 0.000000;
        isotopicDistribution[4][13] = 0.000000;
        isotopicDistribution[4][14] = 0.000000;
        isotopicDistribution[4][15] = 0.000000;
        isotopicDistribution[4][16] = 0.000000;
        isotopicDistribution[4][17] = 0.000000;
        isotopicDistribution[4][18] = 0.000000;
        isotopicDistribution[4][19] = 0.000000;
        isotopicDistribution[5][0] = 0.595410;
        isotopicDistribution[5][1] = 1.000000;
        isotopicDistribution[5][2] = 0.909613;
        isotopicDistribution[5][3] = 0.587263;
        isotopicDistribution[5][4] = 0.299194;
        isotopicDistribution[5][5] = 0.127213;
        isotopicDistribution[5][6] = 0.046724;
        isotopicDistribution[5][7] = 0.015174;
        isotopicDistribution[5][8] = 0.004431;
        isotopicDistribution[5][9] = 0.001178;
        isotopicDistribution[5][10] = 0.000000;
        isotopicDistribution[5][11] = 0.000000;
        isotopicDistribution[5][12] = 0.000000;
        isotopicDistribution[5][13] = 0.000000;
        isotopicDistribution[5][14] = 0.000000;
        isotopicDistribution[5][15] = 0.000000;
        isotopicDistribution[5][16] = 0.000000;
        isotopicDistribution[5][17] = 0.000000;
        isotopicDistribution[5][18] = 0.000000;
        isotopicDistribution[5][19] = 0.000000;
        isotopicDistribution[6][0] = 0.484866;
        isotopicDistribution[6][1] = 0.953177;
        isotopicDistribution[6][2] = 1.000000;
        isotopicDistribution[6][3] = 0.737825;
        isotopicDistribution[6][4] = 0.427036;
        isotopicDistribution[6][5] = 0.205462;
        isotopicDistribution[6][6] = 0.085169;
        isotopicDistribution[6][7] = 0.031161;
        isotopicDistribution[6][8] = 0.010239;
        isotopicDistribution[6][9] = 0.003062;
        isotopicDistribution[6][10] = 0.000000;
        isotopicDistribution[6][11] = 0.000000;
        isotopicDistribution[6][12] = 0.000000;
        isotopicDistribution[6][13] = 0.000000;
        isotopicDistribution[6][14] = 0.000000;
        isotopicDistribution[6][15] = 0.000000;
        isotopicDistribution[6][16] = 0.000000;
        isotopicDistribution[6][17] = 0.000000;
        isotopicDistribution[6][18] = 0.000000;
        isotopicDistribution[6][19] = 0.000000;
        isotopicDistribution[7][0] = 0.368757;
        isotopicDistribution[7][1] = 0.828938;
        isotopicDistribution[7][2] = 1.000000;
        isotopicDistribution[7][3] = 0.851429;
        isotopicDistribution[7][4] = 0.570087;
        isotopicDistribution[7][5] = 0.317911;
        isotopicDistribution[7][6] = 0.152962;
        isotopicDistribution[7][7] = 0.065031;
        isotopicDistribution[7][8] = 0.024851;
        isotopicDistribution[7][9] = 0.008647;
        isotopicDistribution[7][10] = 0.002767;
        isotopicDistribution[7][11] = 0.000000;
        isotopicDistribution[7][12] = 0.000000;
        isotopicDistribution[7][13] = 0.000000;
        isotopicDistribution[7][14] = 0.000000;
        isotopicDistribution[7][15] = 0.000000;
        isotopicDistribution[7][16] = 0.000000;
        isotopicDistribution[7][17] = 0.000000;
        isotopicDistribution[7][18] = 0.000000;
        isotopicDistribution[7][19] = 0.000000;
        isotopicDistribution[8][0] = 0.295369;
        isotopicDistribution[8][1] = 0.745770;
        isotopicDistribution[8][2] = 1.000000;
        isotopicDistribution[8][3] = 0.939923;
        isotopicDistribution[8][4] = 0.691490;
        isotopicDistribution[8][5] = 0.422289;
        isotopicDistribution[8][6] = 0.221980;
        isotopicDistribution[8][7] = 0.102927;
        isotopicDistribution[8][8] = 0.042845;
        isotopicDistribution[8][9] = 0.016224;
        isotopicDistribution[8][10] = 0.005647;
        isotopicDistribution[8][11] = 0.001821;
        isotopicDistribution[8][12] = 0.000000;
        isotopicDistribution[8][13] = 0.000000;
        isotopicDistribution[8][14] = 0.000000;
        isotopicDistribution[8][15] = 0.000000;
        isotopicDistribution[8][16] = 0.000000;
        isotopicDistribution[8][17] = 0.000000;
        isotopicDistribution[8][18] = 0.000000;
        isotopicDistribution[8][19] = 0.000000;
        isotopicDistribution[9][0] = 0.234948;
        isotopicDistribution[9][1] = 0.658283;
        isotopicDistribution[9][2] = 0.971762;
        isotopicDistribution[9][3] = 1.000000;
        isotopicDistribution[9][4] = 0.802247;
        isotopicDistribution[9][5] = 0.532683;
        isotopicDistribution[9][6] = 0.303778;
        isotopicDistribution[9][7] = 0.152559;
        isotopicDistribution[9][8] = 0.068695;
        isotopicDistribution[9][9] = 0.028112;
        isotopicDistribution[9][10] = 0.010566;
        isotopicDistribution[9][11] = 0.003678;
        isotopicDistribution[9][12] = 0.001194;
        isotopicDistribution[9][13] = 0.000000;
        isotopicDistribution[9][14] = 0.000000;
        isotopicDistribution[9][15] = 0.000000;
        isotopicDistribution[9][16] = 0.000000;
        isotopicDistribution[9][17] = 0.000000;
        isotopicDistribution[9][18] = 0.000000;
        isotopicDistribution[9][19] = 0.000000;
        isotopicDistribution[10][0] = 0.178508;
        isotopicDistribution[10][1] = 0.551598;
        isotopicDistribution[10][2] = 0.891811;
        isotopicDistribution[10][3] = 1.000000;
        isotopicDistribution[10][4] = 0.870814;
        isotopicDistribution[10][5] = 0.625779;
        isotopicDistribution[10][6] = 0.385337;
        isotopicDistribution[10][7] = 0.208576;
        isotopicDistribution[10][8] = 0.101081;
        isotopicDistribution[10][9] = 0.044468;
        isotopicDistribution[10][10] = 0.017951;
        isotopicDistribution[10][11] = 0.006706;
        isotopicDistribution[10][12] = 0.002335;
        isotopicDistribution[10][13] = 0.000000;
        isotopicDistribution[10][14] = 0.000000;
        isotopicDistribution[10][15] = 0.000000;
        isotopicDistribution[10][16] = 0.000000;
        isotopicDistribution[10][17] = 0.000000;
        isotopicDistribution[10][18] = 0.000000;
        isotopicDistribution[10][19] = 0.000000;
        isotopicDistribution[11][0] = 0.139841;
        isotopicDistribution[11][1] = 0.470844;
        isotopicDistribution[11][2] = 0.825470;
        isotopicDistribution[11][3] = 1.000000;
        isotopicDistribution[11][4] = 0.938111;
        isotopicDistribution[11][5] = 0.724608;
        isotopicDistribution[11][6] = 0.478745;
        isotopicDistribution[11][7] = 0.277647;
        isotopicDistribution[11][8] = 0.144001;
        isotopicDistribution[11][9] = 0.067736;
        isotopicDistribution[11][10] = 0.029214;
        isotopicDistribution[11][11] = 0.011655;
        isotopicDistribution[11][12] = 0.004331;
        isotopicDistribution[11][13] = 0.001508;
        isotopicDistribution[11][14] = 0.000000;
        isotopicDistribution[11][15] = 0.000000;
        isotopicDistribution[11][16] = 0.000000;
        isotopicDistribution[11][17] = 0.000000;
        isotopicDistribution[11][18] = 0.000000;
        isotopicDistribution[11][19] = 0.000000;
        isotopicDistribution[12][0] = 0.110880;
        isotopicDistribution[12][1] = 0.404250;
        isotopicDistribution[12][2] = 0.764124;
        isotopicDistribution[12][3] = 0.994691;
        isotopicDistribution[12][4] = 1.000000;
        isotopicDistribution[12][5] = 0.825970;
        isotopicDistribution[12][6] = 0.582522;
        isotopicDistribution[12][7] = 0.360096;
        isotopicDistribution[12][8] = 0.198832;
        isotopicDistribution[12][9] = 0.099472;
        isotopicDistribution[12][10] = 0.045591;
        isotopicDistribution[12][11] = 0.019314;
        isotopicDistribution[12][12] = 0.007618;
        isotopicDistribution[12][13] = 0.002814;
        isotopicDistribution[12][14] = 0.000000;
        isotopicDistribution[12][15] = 0.000000;
        isotopicDistribution[12][16] = 0.000000;
        isotopicDistribution[12][17] = 0.000000;
        isotopicDistribution[12][18] = 0.000000;
        isotopicDistribution[12][19] = 0.000000;
        isotopicDistribution[13][0] = 0.081380;
        isotopicDistribution[13][1] = 0.319627;
        isotopicDistribution[13][2] = 0.652303;
        isotopicDistribution[13][3] = 0.918383;
        isotopicDistribution[13][4] = 1.000000;
        isotopicDistribution[13][5] = 0.895645;
        isotopicDistribution[13][6] = 0.685621;
        isotopicDistribution[13][7] = 0.460421;
        isotopicDistribution[13][8] = 0.276378;
        isotopicDistribution[13][9] = 0.150407;
        isotopicDistribution[13][10] = 0.075030;
        isotopicDistribution[13][11] = 0.034611;
        isotopicDistribution[13][12] = 0.014871;
        isotopicDistribution[13][13] = 0.005987;
        isotopicDistribution[13][14] = 0.002270;
        isotopicDistribution[13][15] = 0.000000;
        isotopicDistribution[13][16] = 0.000000;
        isotopicDistribution[13][17] = 0.000000;
        isotopicDistribution[13][18] = 0.000000;
        isotopicDistribution[13][19] = 0.000000;
        isotopicDistribution[14][0] = 0.062581;
        isotopicDistribution[14][1] = 0.263702;
        isotopicDistribution[14][2] = 0.575328;
        isotopicDistribution[14][3] = 0.863482;
        isotopicDistribution[14][4] = 1.000000;
        isotopicDistribution[14][5] = 0.950817;
        isotopicDistribution[14][6] = 0.771506;
        isotopicDistribution[14][7] = 0.548469;
        isotopicDistribution[14][8] = 0.348161;
        isotopicDistribution[14][9] = 0.200189;
        isotopicDistribution[14][10] = 0.105433;
        isotopicDistribution[14][11] = 0.051316;
        isotopicDistribution[14][12] = 0.023251;
        isotopicDistribution[14][13] = 0.009867;
        isotopicDistribution[14][14] = 0.003941;
        isotopicDistribution[14][15] = 0.001489;
        isotopicDistribution[14][16] = 0.000000;
        isotopicDistribution[14][17] = 0.000000;
        isotopicDistribution[14][18] = 0.000000;
        isotopicDistribution[14][19] = 0.000000;
        isotopicDistribution[15][0] = 0.049015;
        isotopicDistribution[15][1] = 0.220215;
        isotopicDistribution[15][2] = 0.510673;
        isotopicDistribution[15][3] = 0.812617;
        isotopicDistribution[15][4] = 0.995728;
        isotopicDistribution[15][5] = 1.000000;
        isotopicDistribution[15][6] = 0.855816;
        isotopicDistribution[15][7] = 0.640920;
        isotopicDistribution[15][8] = 0.428150;
        isotopicDistribution[15][9] = 0.258844;
        isotopicDistribution[15][10] = 0.143228;
        isotopicDistribution[15][11] = 0.073195;
        isotopicDistribution[15][12] = 0.034802;
        isotopicDistribution[15][13] = 0.015490;
        isotopicDistribution[15][14] = 0.006487;
        isotopicDistribution[15][15] = 0.002568;
        isotopicDistribution[15][16] = 0.000000;
        isotopicDistribution[15][17] = 0.000000;
        isotopicDistribution[15][18] = 0.000000;
        isotopicDistribution[15][19] = 0.000000;
        isotopicDistribution[16][0] = 0.037099;
        isotopicDistribution[16][1] = 0.176952;
        isotopicDistribution[16][2] = 0.434584;
        isotopicDistribution[16][3] = 0.730909;
        isotopicDistribution[16][4] = 0.945007;
        isotopicDistribution[16][5] = 1.000000;
        isotopicDistribution[16][6] = 0.900675;
        isotopicDistribution[16][7] = 0.709153;
        isotopicDistribution[16][8] = 0.497628;
        isotopicDistribution[16][9] = 0.315790;
        isotopicDistribution[16][10] = 0.183301;
        isotopicDistribution[16][11] = 0.098210;
        isotopicDistribution[16][12] = 0.048934;
        isotopicDistribution[16][13] = 0.022815;
        isotopicDistribution[16][14] = 0.010006;
        isotopicDistribution[16][15] = 0.004146;
        isotopicDistribution[16][16] = 0.001630;
        isotopicDistribution[16][17] = 0.000000;
        isotopicDistribution[16][18] = 0.000000;
        isotopicDistribution[16][19] = 0.000000;
        isotopicDistribution[17][0] = 0.028489;
        isotopicDistribution[17][1] = 0.143772;
        isotopicDistribution[17][2] = 0.372798;
        isotopicDistribution[17][3] = 0.660813;
        isotopicDistribution[17][4] = 0.899122;
        isotopicDistribution[17][5] = 1.000000;
        isotopicDistribution[17][6] = 0.945614;
        isotopicDistribution[17][7] = 0.780955;
        isotopicDistribution[17][8] = 0.574359;
        isotopicDistribution[17][9] = 0.381740;
        isotopicDistribution[17][10] = 0.231932;
        isotopicDistribution[17][11] = 0.130003;
        isotopicDistribution[17][12] = 0.067735;
        isotopicDistribution[17][13] = 0.033010;
        isotopicDistribution[17][14] = 0.015127;
        isotopicDistribution[17][15] = 0.006548;
        isotopicDistribution[17][16] = 0.002688;
        isotopicDistribution[17][17] = 0.001050;
        isotopicDistribution[17][18] = 0.021473;
        isotopicDistribution[17][19] = 0.000000;
        isotopicDistribution[18][0] = 0.114423;
        isotopicDistribution[18][1] = 0.313596;
        isotopicDistribution[18][2] = 0.588027;
        isotopicDistribution[18][3] = 0.846974;
        isotopicDistribution[18][4] = 0.997820;
        isotopicDistribution[18][5] = 1.000000;
        isotopicDistribution[18][6] = 0.875691;
        isotopicDistribution[18][7] = 0.683169;
        isotopicDistribution[18][8] = 0.481829;
        isotopicDistribution[18][9] = 0.310749;
        isotopicDistribution[18][10] = 0.184951;
        isotopicDistribution[18][11] = 0.102348;
        isotopicDistribution[18][12] = 0.052989;
        isotopicDistribution[18][13] = 0.025802;
        isotopicDistribution[18][14] = 0.011870;
        isotopicDistribution[18][15] = 0.005179;
        isotopicDistribution[18][16] = 0.002150;
        isotopicDistribution[18][17] = 0.000000;
        isotopicDistribution[18][18] = 0.000000;
        isotopicDistribution[18][19] = 0.000000;
        isotopicDistribution[19][0] = 0.016073;
        isotopicDistribution[19][1] = 0.090246;
        isotopicDistribution[19][2] = 0.260110;
        isotopicDistribution[19][3] = 0.512071;
        isotopicDistribution[19][4] = 0.773268;
        isotopicDistribution[19][5] = 0.953917;
        isotopicDistribution[19][6] = 1.000000;
        isotopicDistribution[19][7] = 0.915161;
        isotopicDistribution[19][8] = 0.745551;
        isotopicDistribution[19][9] = 0.548712;
        isotopicDistribution[19][10] = 0.369064;
        isotopicDistribution[19][11] = 0.228958;
        isotopicDistribution[19][12] = 0.132003;
        isotopicDistribution[19][13] = 0.071173;
        isotopicDistribution[19][14] = 0.036079;
        isotopicDistribution[19][15] = 0.017273;
        isotopicDistribution[19][16] = 0.007841;
        isotopicDistribution[19][17] = 0.003386;
        isotopicDistribution[19][18] = 0.001395;
        isotopicDistribution[19][19] = 0.000000;


        loadIntensePeaks();
        loadMostIntensePeaks();
    }
    public static void main(String args []) throws Exception {
        String line = null;
        
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        for(int i = 0; i < 20; i++) {
            double mass = i*500;
            System.out.println(i + "\tmost intense peak for " + mass + ": " + getMostIntensePeak(mass));
            int [] peaks = mass2IntensePeaks(mass);
            for(int j : peaks) {
                System.out.print(j + "\t");
            }
            System.out.println();
        }
        
        while(true) {
            System.out.print("input the mass: ");
            line = br.readLine();
            if("end".equals(line) || "exit".equals(line) || "quit".equals(line)) {
                break;
            }
            double mass = Double.parseDouble(line);
            System.out.print("input the threshold: ");
            line = br.readLine();
            if("end".equals(line) || "exit".equals(line) || "quit".equals(line)) {
                break;
            }
            double threshold = Double.parseDouble(line);
            List<Integer> peaks = IsotopeCalc.mass2IsotopicPeaks(mass, threshold);        
            for(Integer i : peaks) {
                System.out.println(i);
            }
            
            System.out.println();
 
        }
        
    }
}



