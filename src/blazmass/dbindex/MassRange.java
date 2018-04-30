package blazmass.dbindex;

/**
 * Represents a mass range to be passed to sequences query
 *
 * @author Adam
 */
public class MassRange {

    private float precMass;
    private float tolerance;

    /**
     * Create new mass range
     * @param precMass precursor mass
     * @param tolerance mass dependant pre-calculated tolerance for that mass to use
     * example: resulting mass range will be <precMass-tolerance, precMass+tolerance>
     */
    public MassRange(float precMass, float tolerance) {
        this.precMass = precMass;
        this.tolerance = tolerance;
    }

    public float getPrecMass() {
        return precMass;
    }

    public float getTolerance() {
        return tolerance;
    }

    @Override
    public String toString() {
        return "MassRange{" + "precMass=" + precMass + ", tolerance=" + tolerance + '}';
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 29 * hash + Float.floatToIntBits(this.precMass);
        hash = 29 * hash + Float.floatToIntBits(this.tolerance);
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final MassRange other = (MassRange) obj;
        if (Float.floatToIntBits(this.precMass) != Float.floatToIntBits(other.precMass)) {
            return false;
        }
        if (Float.floatToIntBits(this.tolerance) != Float.floatToIntBits(other.tolerance)) {
            return false;
        }
        return true;
    }
    
    
    
    
}
