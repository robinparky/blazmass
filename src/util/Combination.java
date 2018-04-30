/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package util;

/**
 *
 * @author rpark
 */
import java.util.*;
import org.paukov.combinatorics.*;

public class Combination {
    

        public static void main(String[] args) {
                //Integer[] arr = {1, 4, 7, };
            List l = new ArrayList<Integer>();
            l.add(1);
            l.add(2);
            l.add(3);

                Generator<Integer> gen = getCombinations(l);

                
                // Print the subsets
                for (ICombinatoricsVector<Integer> subSet : gen) {
                        System.out.println(subSet);
                }
        }

        /*
         * A set A is a subset of a set B if A is "contained" inside B. A and B may coincide. The relationship of one set being a subset of another is called inclusion or sometimes containment.

        Examples:

        The set (1, 2) is a proper subset of (1, 2, 3).
        Any set is a subset of itself, but not a proper subset.
        The empty set, denoted by ∅, is also a subset of any given set X.
        */
       // public static Generator<Integer> getCombinations(Integer[] arr)
        public static Generator<Integer> getCombinations(List<Integer> l)
        {                
            // Create an initial vector/set
                ICombinatoricsVector<Integer> initialSet = Factory.createVector(l);

                // Create an instance of the subset generator
                Generator<Integer> gen = Factory.createSubSetGenerator(initialSet);

                /*
                // Print the subsets
                for (ICombinatoricsVector<Integer> subSet : gen) {
                        System.out.println(subSet);
                }*/               
                return gen;

        }


    
}
