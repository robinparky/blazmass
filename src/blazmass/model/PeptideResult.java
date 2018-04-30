/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package blazmass.model;

import blazmass.dbindex.IndexedSequence;
import java.util.LinkedList;
import java.util.List;

/**
 *
 * @author rpark
 */
public class PeptideResult implements Comparable {
    
    private IndexedSequence indexedSeq;
    private int pepLength;
    private float peptideMass;
    private int ion;
    private boolean isDecoy;
    private int totalIon;
    private float xCorr=-100f;
    private float zScore;
    private int matchedIon;
    private boolean isModified;        
    
    public class MatchedIon {
        private String ionType;
        private float mz;

        public MatchedIon(String ionType, float mz) {
            this.ionType = ionType;
            this.mz = mz;
        }
        
        public String getIonType() {
            return ionType;
        }

        public void setIonType(String ionType) {
            this.ionType = ionType;
        }

        public float getMz() {
            return mz;
        }

        public void setMz(float mz) {
            this.mz = mz;
        }
        
        
    }
    
    private final List<MatchedIon> matchedIonList = new LinkedList<MatchedIon>();
    
    public int getPepLength() {
        return pepLength;
    }

    public void setPepLength(int pepLength) {
        this.pepLength = pepLength;
    }

    public float getPeptideMass() {
        return peptideMass;
    }

    public void setPeptideMass(float peptideMass) {
        this.peptideMass = peptideMass;
    }
    
    public void addPeptideMass(double peptideMass) {
        this.peptideMass += peptideMass;
    }

    public int getIon() {
        return ion;
    }

    public void setIon(int ion) {
        this.ion = ion;
    }

    public boolean isIsDecoy() {
        return isDecoy;
    }

    public void setIsDecoy(boolean isDecoy) {
        this.isDecoy = isDecoy;
    }

    public int getTotalIon() {
        return totalIon;
    }

    public void setTotalIon(int totalIon) {
        this.totalIon = totalIon;
    }

    public float getxCorr() {
        return xCorr;
    }

    public void setxCorr(float xCorr) {
        this.xCorr = xCorr;
    }
 
    public IndexedSequence getIndexedSeq() {
        return indexedSeq;
    }

    public void setIndexedSeq(IndexedSequence indexedSeq) {
        this.indexedSeq = indexedSeq;
    }

    public int compareTo(Object obj) {
        
        if(null == obj)
            return 1;
        
        PeptideResult pr = (PeptideResult) obj;

        float diff = pr.getxCorr()-xCorr;
        
        if( diff == 0 ) {
            return 0;
        } else if(diff < 0) {
            return -1;

        } else {
            return 1;
        }
    }

    public float getzScore() {
        return zScore;
    }

    public void setzScore(float zScore) {
        this.zScore = zScore;
    }

    public int getMatchedIon() {
        return matchedIon;
    }

    public void setMatchedIon(int matchedIon) {
        this.matchedIon = matchedIon;
    }

    public boolean isIsModified() {
        return isModified;
    }

    public void setIsModified(boolean isModified) {
        this.isModified = isModified;
    }

    public void addMatchedIon(String ionType, float mz) {
        this.matchedIonList.add(new MatchedIon(ionType, mz));
    }

    public List<MatchedIon> getMatchedIonList() {
        return matchedIonList;
    }


    
}
