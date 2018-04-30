package blazmass;

import blazmass.io.SearchParams;
import java.util.Set;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author rpark2
 */
public class HighResMassProcessor {

    private int[] mappingArr=null;
    private int binSize=0;
    
    
    public static void assignTheoMass(int[] theorArr, float mass, SearchParams sparam, int cs, HighResMassProcessor hprocess, int value, Set<Integer> theoNonZeroInd) {
              
        if(mass<=0) return;
	//float nmass = (mass + AssignMass.getH())/cs;
        float newMass = (mass+AssignMass.getH()*(cs-1))/cs;

        int index = hprocess.getBinIndex(newMass, sparam);                        
//	System.out.println("====================== \t" + newMass + "\t" + index);
	if(index<=0) return;

        
      //  System.out.println("===" + index + " " + value);
        
        
        theorArr[index] = value;
        theoNonZeroInd.add(index);
        
        if(theorArr[index-1] < value) {
            theorArr[index-1] = value;
            theoNonZeroInd.add(index-1);
        }        
        if(theorArr[index+1] <  value) {
            theorArr[index+1] = value;
            theoNonZeroInd.add(index+1);
        }    
        
/*
        
        
        int index = (int)(mass);
        //if(index==243)
        //System.out.println("==>>" + mass + " " + index + " " + theorMass[index] + " " + minValue);
        if(theorMass[index] < minValue) {
            theorMass[index] = minValue;
        }
        
        */
        
    }
    public static void assignTheoMassTest(int[] theorArr, float[] fragArr, SearchParams sparam, int cs, HighResMassProcessor hprocess, String ion) {
        
        for(float mass:fragArr)
            System.out.println("ion \t" + ion + "\t"+ mass + "\t" + (mass+AssignMass.getH()*(cs-1))/cs);

    }
    public static void assignTheoMassTest(int[] theorArr, float mass, SearchParams sparam, int cs, HighResMassProcessor hprocess, String ion) {
        
            System.out.println("ion \t" + ion + "\t"+ mass + "\t" + (mass+AssignMass.getH()*(cs-1))/cs);

    }    
    

        
        
    public static void assignTheoMass(int[] theorArr, float mass, SearchParams sparam, int cs, HighResMassProcessor hprocess, float weight, Set<Integer> theoNonZeroInd) {
        
        float newMass = (mass+AssignMass.getH()*(cs-1))/cs;
        
        int averagineIndex = (int)(newMass/500);

       //	if(averagineIndex==1)
       // System.out.println("===mass\t" + newMass + "\t" + mass + " " +  cs + " " + (mass + Constants.MADD_DIFF_C12C13)/cs);
        

//This charge state calc is in coreect!!
        switch(averagineIndex) {
            case 0 : 
                     assignValue(theorArr, newMass, 50, sparam, hprocess, weight,theoNonZeroInd);
                     assignValue(theorArr, newMass + Constants.MADD_DIFF_C12C13/cs, 14, sparam, hprocess, weight,theoNonZeroInd);
                     break;
            case 1 : 
                     assignValue(theorArr, newMass, 50, sparam, hprocess, weight,theoNonZeroInd);
                     assignValue(theorArr, newMass + Constants.MADD_DIFF_C12C13/cs, 28, sparam, hprocess, weight, theoNonZeroInd);
                     break;
            case 2 : 
                     assignValue(theorArr, newMass, 50, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + Constants.MADD_DIFF_C12C13/cs, 43, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 2*Constants.MADD_DIFF_C12C13/cs, 22, sparam, hprocess, weight, theoNonZeroInd);
                     break;
            case 3 : 
                     assignValue(theorArr, newMass, 45, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + Constants.MADD_DIFF_C12C13/cs, 50, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 2*Constants.MADD_DIFF_C12C13/cs, 30, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 3*Constants.MADD_DIFF_C12C13/cs, 15, sparam, hprocess, weight, theoNonZeroInd);
                     break;
            case 4 : 
                     assignValue(theorArr, newMass, 36, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + Constants.MADD_DIFF_C12C13/cs, 50, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 2*Constants.MADD_DIFF_C12C13/cs, 39, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 3*Constants.MADD_DIFF_C12C13/cs, 22, sparam, hprocess, weight, theoNonZeroInd);
                     break;
            case 5 : 
                     assignValue(theorArr, newMass, 30, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + Constants.MADD_DIFF_C12C13/cs, 50, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 2*Constants.MADD_DIFF_C12C13/cs, 46, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 3*Constants.MADD_DIFF_C12C13/cs, 30, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 4*Constants.MADD_DIFF_C12C13/cs, 15, sparam, hprocess, weight, theoNonZeroInd);
                     break;
            case 6 : 
                     assignValue(theorArr, newMass, 24, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + Constants.MADD_DIFF_C12C13/cs, 47, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 2*Constants.MADD_DIFF_C12C13/cs, 50, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 3*Constants.MADD_DIFF_C12C13/cs, 37, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 4*Constants.MADD_DIFF_C12C13/cs, 22, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 5*Constants.MADD_DIFF_C12C13/cs, 11, sparam, hprocess, weight, theoNonZeroInd);
                     break;
            case 7 : 
                     assignValue(theorArr, newMass, 19, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + Constants.MADD_DIFF_C12C13/cs, 42, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 2*Constants.MADD_DIFF_C12C13/cs, 50, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 3*Constants.MADD_DIFF_C12C13/cs, 43, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 4*Constants.MADD_DIFF_C12C13/cs, 29, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 5*Constants.MADD_DIFF_C12C13/cs, 15, sparam, hprocess, weight, theoNonZeroInd);
                     break;
            case 8 : 
                     assignValue(theorArr, newMass, 15, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + Constants.MADD_DIFF_C12C13/cs, 38, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 2*Constants.MADD_DIFF_C12C13/cs, 50, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 3*Constants.MADD_DIFF_C12C13/cs, 47, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 4*Constants.MADD_DIFF_C12C13/cs, 35, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 5*Constants.MADD_DIFF_C12C13/cs, 21, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 6*Constants.MADD_DIFF_C12C13/cs, 11, sparam, hprocess, weight, theoNonZeroInd);
                     break;
            case 9 : 
                     assignValue(theorArr, newMass, 12, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + Constants.MADD_DIFF_C12C13/cs, 33, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 2*Constants.MADD_DIFF_C12C13/cs, 49, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 3*Constants.MADD_DIFF_C12C13/cs, 50, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 4*Constants.MADD_DIFF_C12C13/cs, 40, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 5*Constants.MADD_DIFF_C12C13/cs, 27, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 6*Constants.MADD_DIFF_C12C13/cs, 15, sparam, hprocess, weight, theoNonZeroInd);
                     break;
            case 10 : 
                     assignValue(theorArr, newMass, 9, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + Constants.MADD_DIFF_C12C13/cs, 28, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 2*Constants.MADD_DIFF_C12C13/cs, 45, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 3*Constants.MADD_DIFF_C12C13/cs, 50, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 4*Constants.MADD_DIFF_C12C13/cs, 44, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 5*Constants.MADD_DIFF_C12C13/cs, 31, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 6*Constants.MADD_DIFF_C12C13/cs, 19, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 7*Constants.MADD_DIFF_C12C13/cs, 11, sparam, hprocess, weight, theoNonZeroInd);
                     break;
            case 11 : 
                     assignValue(theorArr, newMass + Constants.MADD_DIFF_C12C13/cs, 24, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 2*Constants.MADD_DIFF_C12C13/cs, 42, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 3*Constants.MADD_DIFF_C12C13/cs, 50, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 4*Constants.MADD_DIFF_C12C13/cs, 47, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 5*Constants.MADD_DIFF_C12C13/cs, 36, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 6*Constants.MADD_DIFF_C12C13/cs, 24, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 7*Constants.MADD_DIFF_C12C13/cs, 14, sparam, hprocess, weight, theoNonZeroInd);
                     break;
            case 12 : 
                     assignValue(theorArr, newMass + Constants.MADD_DIFF_C12C13/cs, 20, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 2*Constants.MADD_DIFF_C12C13/cs, 38, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 3*Constants.MADD_DIFF_C12C13/cs, 50, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 4*Constants.MADD_DIFF_C12C13/cs, 50, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 5*Constants.MADD_DIFF_C12C13/cs, 42, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 6*Constants.MADD_DIFF_C12C13/cs, 29, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 7*Constants.MADD_DIFF_C12C13/cs, 18, sparam, hprocess, weight, theoNonZeroInd);
                     break;
            default : 
                     assignValue(theorArr, newMass + 3*Constants.MADD_DIFF_C12C13/cs, 50, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 4*Constants.MADD_DIFF_C12C13/cs, 50, sparam, hprocess, weight, theoNonZeroInd);
                     assignValue(theorArr, newMass + 5*Constants.MADD_DIFF_C12C13/cs, 42, sparam, hprocess, weight, theoNonZeroInd);
                     break;
        }


    }

/*        
    public static void assignTheoMass(int[] theorArr, float mass, SearchParams sparam, int cs, HighResMassProcessor hprocess) {
        
        int averagineIndex = (int)mass/500;
        
//	if(averagineIndex==1)
        System.out.println("===mass\t" + mass + " " +  cs + " " + (mass + Constants.MADD_DIFF_C12C13)/cs);
        

//This charge state calc is in coreect!!
        switch(averagineIndex) {
            case 0 : 
                     assignValue(theorArr, mass/cs, 50, sparam, hprocess);
                     assignValue(theorArr, (mass + Constants.MADD_DIFF_C12C13)/cs, 14, sparam, hprocess);
                     break;
            case 1 : 
                     assignValue(theorArr, mass/cs, 50, sparam, hprocess);
                     assignValue(theorArr, (mass + Constants.MADD_DIFF_C12C13)/cs, 28, sparam, hprocess);
                     break;
            case 2 : 
                     assignValue(theorArr, mass/cs, 50, sparam, hprocess);
                     assignValue(theorArr, (mass + Constants.MADD_DIFF_C12C13)/cs, 43, sparam, hprocess);
                     assignValue(theorArr, (mass + 2*Constants.MADD_DIFF_C12C13)/cs, 22, sparam, hprocess);
                     break;
            case 3 : 
                     assignValue(theorArr, mass/cs, 45, sparam, hprocess);
                     assignValue(theorArr, (mass + Constants.MADD_DIFF_C12C13)/cs, 50, sparam, hprocess);
                     assignValue(theorArr, (mass + 2*Constants.MADD_DIFF_C12C13)/cs, 30, sparam, hprocess);
                     assignValue(theorArr, (mass + 3*Constants.MADD_DIFF_C12C13)/cs, 15, sparam, hprocess);
                     break;
            case 4 : 
                     assignValue(theorArr, mass/cs, 36, sparam, hprocess);
                     assignValue(theorArr, (mass + Constants.MADD_DIFF_C12C13)/cs, 50, sparam, hprocess);
                     assignValue(theorArr, (mass + 2*Constants.MADD_DIFF_C12C13)/cs, 39, sparam, hprocess);
                     assignValue(theorArr, (mass + 3*Constants.MADD_DIFF_C12C13)/cs, 22, sparam, hprocess);
                     break;
            case 5 : 
                     assignValue(theorArr, mass/cs, 30, sparam, hprocess);
                     assignValue(theorArr, (mass + Constants.MADD_DIFF_C12C13)/cs, 50, sparam, hprocess);
                     assignValue(theorArr, (mass + 2*Constants.MADD_DIFF_C12C13)/cs, 46, sparam, hprocess);
                     assignValue(theorArr, (mass + 3*Constants.MADD_DIFF_C12C13)/cs, 30, sparam, hprocess);
                     assignValue(theorArr, (mass + 4*Constants.MADD_DIFF_C12C13)/cs, 15, sparam, hprocess);
                     break;
            case 6 : 
                     assignValue(theorArr, mass/cs, 24, sparam, hprocess);
                     assignValue(theorArr, (mass + Constants.MADD_DIFF_C12C13)/cs, 47, sparam, hprocess);
                     assignValue(theorArr, (mass + 2*Constants.MADD_DIFF_C12C13)/cs, 50, sparam, hprocess);
                     assignValue(theorArr, (mass + 3*Constants.MADD_DIFF_C12C13)/cs, 37, sparam, hprocess);
                     assignValue(theorArr, (mass + 4*Constants.MADD_DIFF_C12C13)/cs, 22, sparam, hprocess);
                     assignValue(theorArr, (mass + 5*Constants.MADD_DIFF_C12C13)/cs, 11, sparam, hprocess);
                     break;
            case 7 : 
                     assignValue(theorArr, mass/cs, 19, sparam, hprocess);
                     assignValue(theorArr, (mass + Constants.MADD_DIFF_C12C13)/cs, 42, sparam, hprocess);
                     assignValue(theorArr, (mass + 2*Constants.MADD_DIFF_C12C13)/cs, 50, sparam, hprocess);
                     assignValue(theorArr, (mass + 3*Constants.MADD_DIFF_C12C13)/cs, 43, sparam, hprocess);
                     assignValue(theorArr, (mass + 4*Constants.MADD_DIFF_C12C13)/cs, 29, sparam, hprocess);
                     assignValue(theorArr, (mass + 5*Constants.MADD_DIFF_C12C13)/cs, 15, sparam, hprocess);
                     break;
            case 8 : 
                     assignValue(theorArr, mass/cs, 15, sparam, hprocess);
                     assignValue(theorArr, (mass + Constants.MADD_DIFF_C12C13)/cs, 38, sparam, hprocess);
                     assignValue(theorArr, (mass + 2*Constants.MADD_DIFF_C12C13)/cs, 50, sparam, hprocess);
                     assignValue(theorArr, (mass + 3*Constants.MADD_DIFF_C12C13)/cs, 47, sparam, hprocess);
                     assignValue(theorArr, (mass + 4*Constants.MADD_DIFF_C12C13)/cs, 35, sparam, hprocess);
                     assignValue(theorArr, (mass + 5*Constants.MADD_DIFF_C12C13)/cs, 21, sparam, hprocess);
                     assignValue(theorArr, (mass + 6*Constants.MADD_DIFF_C12C13)/cs, 11, sparam, hprocess);
                     break;
            case 9 : 
                     assignValue(theorArr, mass/cs, 12, sparam, hprocess);
                     assignValue(theorArr, (mass + Constants.MADD_DIFF_C12C13)/cs, 33, sparam, hprocess);
                     assignValue(theorArr, (mass + 2*Constants.MADD_DIFF_C12C13)/cs, 49, sparam, hprocess);
                     assignValue(theorArr, (mass + 3*Constants.MADD_DIFF_C12C13)/cs, 50, sparam, hprocess);
                     assignValue(theorArr, (mass + 4*Constants.MADD_DIFF_C12C13)/cs, 40, sparam, hprocess);
                     assignValue(theorArr, (mass + 5*Constants.MADD_DIFF_C12C13)/cs, 27, sparam, hprocess);
                     assignValue(theorArr, (mass + 6*Constants.MADD_DIFF_C12C13)/cs, 15, sparam, hprocess);
                     break;
            case 10 : 
                     assignValue(theorArr, mass/cs, 9, sparam, hprocess);
                     assignValue(theorArr, (mass + Constants.MADD_DIFF_C12C13)/cs, 28, sparam, hprocess);
                     assignValue(theorArr, (mass + 2*Constants.MADD_DIFF_C12C13)/cs, 45, sparam, hprocess);
                     assignValue(theorArr, (mass + 3*Constants.MADD_DIFF_C12C13)/cs, 50, sparam, hprocess);
                     assignValue(theorArr, (mass + 4*Constants.MADD_DIFF_C12C13)/cs, 44, sparam, hprocess);
                     assignValue(theorArr, (mass + 5*Constants.MADD_DIFF_C12C13)/cs, 31, sparam, hprocess);
                     assignValue(theorArr, (mass + 6*Constants.MADD_DIFF_C12C13)/cs, 19, sparam, hprocess);
                     assignValue(theorArr, (mass + 7*Constants.MADD_DIFF_C12C13)/cs, 11, sparam, hprocess);
                     break;
            case 11 : 
                     assignValue(theorArr, (mass + Constants.MADD_DIFF_C12C13)/cs, 24, sparam, hprocess);
                     assignValue(theorArr, (mass + 2*Constants.MADD_DIFF_C12C13)/cs, 42, sparam, hprocess);
                     assignValue(theorArr, (mass + 3*Constants.MADD_DIFF_C12C13)/cs, 50, sparam, hprocess);
                     assignValue(theorArr, (mass + 4*Constants.MADD_DIFF_C12C13)/cs, 47, sparam, hprocess);
                     assignValue(theorArr, (mass + 5*Constants.MADD_DIFF_C12C13)/cs, 36, sparam, hprocess);
                     assignValue(theorArr, (mass + 6*Constants.MADD_DIFF_C12C13)/cs, 24, sparam, hprocess);
                     assignValue(theorArr, (mass + 7*Constants.MADD_DIFF_C12C13)/cs, 14, sparam, hprocess);
                     break;
            case 12 : 
                     assignValue(theorArr, (mass + Constants.MADD_DIFF_C12C13)/cs, 20, sparam, hprocess);
                     assignValue(theorArr, (mass + 2*Constants.MADD_DIFF_C12C13)/cs, 38, sparam, hprocess);
                     assignValue(theorArr, (mass + 3*Constants.MADD_DIFF_C12C13)/cs, 50, sparam, hprocess);
                     assignValue(theorArr, (mass + 4*Constants.MADD_DIFF_C12C13)/cs, 50, sparam, hprocess);
                     assignValue(theorArr, (mass + 5*Constants.MADD_DIFF_C12C13)/cs, 42, sparam, hprocess);
                     assignValue(theorArr, (mass + 6*Constants.MADD_DIFF_C12C13)/cs, 29, sparam, hprocess);
                     assignValue(theorArr, (mass + 7*Constants.MADD_DIFF_C12C13)/cs, 18, sparam, hprocess);
                     break;
            default : 
                     assignValue(theorArr, (mass + 3*Constants.MADD_DIFF_C12C13)/cs, 50, sparam, hprocess);
                     assignValue(theorArr, (mass + 4*Constants.MADD_DIFF_C12C13)/cs, 50, sparam, hprocess);
                     assignValue(theorArr, (mass + 5*Constants.MADD_DIFF_C12C13)/cs, 42, sparam, hprocess);
                     break;
        }


    }
*/

    public static void assignTheoMass(int[] theorArr, float[] fragArr, SearchParams sparam, int cs, HighResMassProcessor hprocess, float weight, Set<Integer> theoNonZeroInd) {
        for(float mass:fragArr)
            assignTheoMass(theorArr, mass, sparam, cs, hprocess, weight, theoNonZeroInd);             
    }

    public static void assignTheoMass(int[] theorArr, float[] fragArr, SearchParams sparam, int cs, HighResMassProcessor hprocess, int value, Set<Integer> theoNonZeroInd) {
        for(float mass:fragArr)
            assignTheoMass(theorArr, mass, sparam, cs, hprocess, value, theoNonZeroInd);             
    }
    
    public static void assignValue(int theorArr[], float mass, int minValue, SearchParams sparam, HighResMassProcessor hprocess, float weight, Set<Integer> theoNonZeroInd) {
        //int index = (int)(mass*sparam.getFragmentIonToleranceBinScale());         
       // int index = getBinIndex(mass)       
        
       // System.out.println("mass\t" + mass);
        
        minValue = (int)(minValue*weight);
        
	if(mass<=0) return;
     
        int index = hprocess.getBinIndex(mass, sparam);              
         
	if(index<=0) return;
        //System.out.println("mass change " + mass + " " + index + " " + minValue + " " + theorArr[index]);
        if(theorArr[index] < minValue) {
            theorArr[index] = minValue;
            theoNonZeroInd.add(index);
        }        

        if(theorArr[index-1] < minValue) {
            theorArr[index-1] = minValue;
            theoNonZeroInd.add(index-1);
        }        
        if(theorArr[index+1] < minValue) {
            theorArr[index+1] = minValue;
            theoNonZeroInd.add(index+1);
        }        
    }    
    
    public int getBinSize(SearchParams sParam) {
        if(mappingArr==null) getBinMappingArr(sParam);
        
        return binSize;        
    }
    
    public int getBinIndex(float mass, SearchParams sParam) {
        if(mappingArr==null) getBinMappingArr(sParam);
        
        mass += 0.0005;
        mass *= 1000;
        int massIndex = (int)mass;
        
        int binIndex = mappingArr[massIndex];
        
        return binIndex;
    }
    
    public int[] getBinMappingArr(SearchParams sParam) {
        
        if(mappingArr!=null) return mappingArr;
        
        
        double startMass=50;
        double endMass = 8000;
        double convertFactor = sParam.getFragmentIonTolerance()/1000000;

        mappingArr = new int[(int)endMass*1000];
       
        double startRange = startMass-startMass*convertFactor;
        double endRange = startMass+startMass*convertFactor;
        
        int startInt = (int)((startRange+0.0005)*1000);
        int endInt = (int)((endRange+0.0005)*1000);
            
        //System.out.println(convertFactor);
        while(true) {           
            for(int i=startInt;i<=endInt;i++)   
                mappingArr[i] = binSize;
            
            startRange = endRange + 0.001;
            endRange = endRange + endRange*convertFactor;            
            startInt = endInt+1;
            endInt = (int)((endRange+0.0005)*1000);
            
            binSize++;
            
            if(endRange>endMass) break;
            //if(count>5) break;
        }
        
        
     //   for(int i=mappingArr.length-1;i<0;i--)
            //if(mappingArr[i]>0) {
     //           System.out.println(i  + "\t" + mappingArr.length + " " + mappingArr[i]);
                
    //        }
        
        
        /*
        
        //System.out.println(mappingArr[301980]);
        System.out.println(mappingArr[6000000-7000]);
        System.out.println(mappingArr.length);
        System.out.println(mappingArr[310009]);
        for(int i=310009;i<(6000000-7000);i++)
            if(mappingArr[i]<=0) System.out.print("=");
        
        */

/*
for(int i=0;i<mappingArr.length;i++)
System.out.println(i+ "\t" + printMappingRange(i));        
        System.exit(0);
*/
        return mappingArr;
    }
    

    public String printMappingRange(int index) {
        String indexValue="";
        for(int i=0;i<this.mappingArr.length;i++) {
            if( this.mappingArr[i] == index ) {
                indexValue += i;
             
                for(int j=i+1;j<this.mappingArr.length;j++) {
                    if(mappingArr[j] != index) {
                        indexValue += "\t" + (j-1);
                        
                        return indexValue;
                    }
                    
                }
                    
                   
            }
                
        }
        return null;
    }
}
