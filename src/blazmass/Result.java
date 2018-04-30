
package blazmass;
public class Result implements Comparable<Result> {
	private double deltCN;
	private double pepMass;
	private double rawMass;
	private double zScore;
	private double xCorr;
	private int lengthSeq;
	private int matchedIons;
	private int totalIons;
	private boolean decoy;
	private int MAX_PEPTIDE_SIZE=512;
	private int duplicateCount;
	private int[] startingPosition ;	
	private int[] piDiffSearchSites = new int[MAX_PEPTIDE_SIZE];
	private String peptide;
	private  char[] prevNextAA;
	private  String[] reference;
	public Result() {
	}
	public Result(double dDeltCN, double dXCorr, double dZScore,
        	int duplicateCoun, int NUM_REFERENCE) {
		this.prevNextAA = new char[4];
		this.deltCN = dDeltCN;
		this.xCorr = dXCorr;
		this.zScore = dZScore;
		this.duplicateCount = duplicateCount;		
		this.startingPosition = new int[NUM_REFERENCE];
		for (int i = 0; i < NUM_REFERENCE; i++) {
			startingPosition[i] = 0;
		}
		this.reference=new String[NUM_REFERENCE];
		this.peptide="";
	}
    public int compareTo(Result r) {
        //if(r == null) return -1;
        
        double d = r.getxCorr() - this.getxCorr();
        if(d>0) return 1;
        else if(d<0) return -1;
        else return 0;
    }
    public int getMAX_PEPTIDE_SIZE() {
        return MAX_PEPTIDE_SIZE;
    }
    public void setMAX_PEPTIDE_SIZE(int MAX_PEPTIDE_SIZE) {
        this.MAX_PEPTIDE_SIZE = MAX_PEPTIDE_SIZE;
    }
    public double getDeltCN() {
        return deltCN;
    }
    public void setDeltCN(double deltCN) {
        this.deltCN = deltCN;
    }
    public int getDuplicateCount() {
        return duplicateCount;
    }
    public void setDuplicateCount(int duplicateCount) {
        this.duplicateCount = duplicateCount;
    }
    public boolean isDecoy() {
        return decoy;
    }
    public void setDecoy(boolean decoy) {
        this.decoy = decoy;
    }
    
    public int getLengthSeq() {
        return lengthSeq;
    }
    public void setLengthSeq(int lengthSeq) {
        this.lengthSeq = lengthSeq;
    }
    public int getMatchedIons() {
        return matchedIons;
    }
    public void setMatchedIons(int matchedIons) {
        this.matchedIons = matchedIons;
    }
    public double getPepMass() {
        return pepMass;
    }
    public void setPepMass(double pepMass) {
        this.pepMass = pepMass;
    }
    public String getPeptide() {
        return peptide;
    }
    public void setPeptide(String peptide) {
        this.peptide = peptide;
    }
    public int[] getPiDiffSearchSites() {
        return piDiffSearchSites;
    }
    public void setPiDiffSearchSites(int[] piDiffSearchSites) {
        this.piDiffSearchSites = piDiffSearchSites;
    }
    public char[] getPrevNextAA() {
        return prevNextAA;
    }
    public void setPrevNextAA(char[] prevNextAA) {
        this.prevNextAA = prevNextAA;
    }
    
    public void setPrevNextAA(int pos, char ch) {
        this.prevNextAA[pos] = ch;
    }    
    public double getRawMass() {
        return rawMass;
    }
    public void setRawMass(double rawMass) {
        this.rawMass = rawMass;
    }
    public String[] getReference() {
        return reference;
    }
    public void setReference(String[] reference) {
        this.reference = reference;
    }
    public int[] getStartingPosition() {
        return startingPosition;
    }
    public void setStartingPosition(int[] startingPosition) {
        this.startingPosition = startingPosition;
    }
    public int getTotalIons() {
        return totalIons;
    }
    public void setTotalIons(int totalIons) {
        this.totalIons = totalIons;
    }
    public double getxCorr() {
        return xCorr;
    }
    public void setxCorr(double xCorr) {
        this.xCorr = xCorr;
    }
    public double getzScore() {
        return zScore;
    }
    public void setzScore(double zScore) {
        this.zScore = zScore;
    }
    
    
}
