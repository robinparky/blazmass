package blazmass;
/**
 *
 * @author rpark
 */

public class Constants {
    
    
    public static float PROTON = 1.007276466f;
    public static float H = 1.007825f;
    public static final float H2O = 15.9949145f + 2*H;
    public static final float H2O_PROTON = H2O + PROTON;
    public static final float H2O_PROTON_SCALED_DOWN = H2O_PROTON*1000;
    public static final int MIN_PEP_LENGTH=6;
    //public static final int MIN_FRAG_NUM=10;
//    public static final int MIN_PRECURSOR=600;
//    public static final int MAX_PRECURSOR=6000;
    
    public static final int PREFIX_RESIDUES=4;
    public static final int PREFIX_RESIDUES2=8;
    
    public static final float MADD_DIFF_C12C13 = 1.003354826f;
    public static final int MADD_DIFF_C12C13_PPM = 1003;
    
    //////////SCORING //////////////
    public static final int SCORE_BIN_SIZE=200;
    public static final double MAX_PRECURSOR_MASS = 6000.0; 
    //public static final int MAX_CHARGE_STATE = 6; 
    public static final int MAX_CHARGE_STATE_PRECURSOR = 2; 
    //public static final int MAX_CHARGE_STATE = 16; 

    
    /////////OUTPUT ///////////////////
    public static final int PEPTIDE_NUM_DISPLAY = 5; 
    public static final int NUM_REFERENCE = 50; 
   // public static final int CON_SIZE = 7;
    

    /*    
    public static final float log10 = 0.434294481f;
    //public static final float MASSH2 = MASSH * 2;
    // MASSH2O + MASSPROTON
    public static final float MASSH3O = MASSH2O + MASSPROTON;
    public static final float DBINWIDTH = 1.f/1.0005079f;
    public static final float MASSHDB = MASSH*DBINWIDTH;
    public static final float MASSPROTONDB = MASSPROTON*DBINWIDTH;
    public static final float MASSH3ODB = MASSH3O*DBINWIDTH;
    public static final String AVGISOTOPE = "avg";
    public static final String MONOISOTOPE = "mono";
    // Thermo source is about 1.002806, which does not look right
    public static final float MASSDIFFC12C13 = 1.003354826f;
    public static final float MASSDIFFN14N15 = 0.997034968f;
    public static final int NUMCHARS = 256;
 */

    public static final int MAX_INDEX_RESIDUE_LEN = 3; //how many max. residue ions on each side of sequence to store in index    
}
