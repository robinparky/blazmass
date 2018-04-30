/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package util;

/**
 *
 * @author rpark2
 */
public class MathUtil {
    
    public static int round(double d) {

        if (d > 0) {
            return (int) (d + 0.5d);
        } else {
            return (int) (d - 0.5d);
        }
    }
    
}
