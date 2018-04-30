/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package blazmass.io;

/**
 *
 * @author rpark
 */

public class ModResidue {

    private char residue;
    private float massShift;

    public ModResidue(char residue, float massShift)
    {
        this.residue = residue;
        this.massShift = massShift;
    }

    @Override
    public String toString() {
        return "ModResidue{" + "residue=" + residue + ", massShift=" + massShift + '}';
    }

    
    
    public char getResidue() {
        return residue;
    }

    public void setResidue(char residue) {
        this.residue = residue;
    }

    public void setMassShift(float massShift) {
        this.massShift = massShift;
    }

    public float getMassShift() {
        return massShift;
    }

}
