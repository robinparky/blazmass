/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package blazmass.io;

/**
 *
 * @author rpark2
 */
public class MassIntensityModel implements Comparable {
    private float intensity;
    private int index;

    public MassIntensityModel(float intensity, int index) {
        this.intensity = intensity;
        this.index = index;
    }
    
    public float getIntensity() {
        return intensity;
    }

    public void setIntensity(float intensity) {
        this.intensity = intensity;
    }

    public int getIndex() {
        return index;
    }

    public void setIndex(int index) {
        this.index = index;
    }
    
    public int compareTo(Object o) throws ClassCastException {

        MassIntensityModel m = (MassIntensityModel)o;

        if(m.getIntensity()<this.intensity)
            return 1;
        else if(m.getIntensity()>this.intensity)
            return -1;
        else
            return 0;        
    }
}
