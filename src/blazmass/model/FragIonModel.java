/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package blazmass.model;

import blazmass.AssignMass;

/**
 *
 * @author robin
 */
public class FragIonModel {
        private int ion;
        private String ionStr;
        private int weight;
        private float[] fragArr;
        private String modSequence;
        private boolean forwardIon;
        
        public FragIonModel(int ion, String ionStr, int weight, float[] fragArr, boolean forwardIon) {
            this.ion = ion;
            this.ionStr = ionStr;
            this.weight = weight;
            this.fragArr = fragArr;            
            this.forwardIon = forwardIon;
        }
        
        public FragIonModel(int ion, String ionStr, int weight, float[] fragArr, boolean forwardIon, String modSequence) {
            this(ion, ionStr, weight, fragArr, forwardIon);
            this.modSequence = modSequence;
        }
        
        public int getIon() {
            return ion;
        }

        public void setIon(int ion) {
            this.ion = ion;
        }

        public float[] getFragArr() {
            return fragArr;
        }

        public void setFragArr(float[] fragArr) {
            this.fragArr = fragArr;
        }

        public int getWeight() {
            return weight;
        }

        public void setWeight(int weight) {
            this.weight = weight;
        }

    public String getModSequence() {
        return modSequence;
    }

    public void setModSequence(String modSequence) {
        this.modSequence = modSequence;
    }

    public String getIonStr() {
        return ionStr;
    }

    public void setIonStr(String ionStr) {
        this.ionStr = ionStr;
    }

    public boolean isForwardIon() {
        return forwardIon;
    }

    public void setForwardIon(boolean forwardIon) {
        this.forwardIon = forwardIon;
    }

    
    
        
}
