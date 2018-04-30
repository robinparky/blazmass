
package blazmass;

public class PDbEntry {
	          /* pDbEntry stores temporary name & sequence info. */
	          /* for each protein in the database to be analyzed */
	    int iLengthSeq;
	  //  int iMaxLenAllocated;
	   String szDescript;
	   String szSequence;
	   public PDbEntry(){
		   
	   }
	public int getiLengthSeq() {
		return iLengthSeq;
	}
	public void setiLengthSeq(int iLengthSeq) {
		this.iLengthSeq = iLengthSeq;
	}
	/*public int getiMaxLenAllocated() {
		return iMaxLenAllocated;
	}
	public void setiMaxLenAllocated(int iMaxLenAllocated) {
		this.iMaxLenAllocated = iMaxLenAllocated;
	}*/
	public String getSzDescript() {
		return szDescript;
	}
	public void setSzDescript(String szDescript) {
		this.szDescript = szDescript;
	}
	public String getSzSequence() {
		return szSequence;
	}
	public void setSzSequence(String szSequence) {
		this.szSequence = szSequence;
		setiLengthSeq(this.szSequence.length());
	}
	

}
