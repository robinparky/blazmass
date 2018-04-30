package blazmass.io;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Collections;


/**
 * Represents MS2 scan
 *
 * @author Adam
 */
public class MS2Scan {

    private static final int DEFAULT_SIZE = 500;
    private List<Float> masses;
    private List<Float> intensities;
    private int size; //num masses/intensities
    private int isScan1;
    private int isScan2;
    private List<Integer> chargeStates;
    private List<Float> precursorMasses;
    private String scanType;

    public MS2Scan() {
        size = 0;
        masses = new ArrayList<Float>(DEFAULT_SIZE);
        intensities = new ArrayList<Float>(DEFAULT_SIZE);
        chargeStates = new ArrayList<Integer>();
        precursorMasses = new ArrayList<Float>();
    }

    void addMassIntensity(float mass, float intensity) {
        masses.add(mass);
        intensities.add(intensity);
        ++size;
    }

    void setIsScan1(int isScan1) {
        this.isScan1 = isScan1;
    }

    void setIsScan2(int isScan2) {
        this.isScan2 = isScan2;
    }

    void addChargeState(int chargeState) {
        chargeStates.add(chargeState);
    }

    void addPrecursorMass(float precursorMass) {
        precursorMasses.add(precursorMass);
    }

    public List<Float> getMasses() {
        return masses;
    }

    public List<Float> getIntensities() {
        return intensities;
    }

    public int getSize() {
        return size;
    }

    public int getIsScan1() {
        return isScan1;
    }

    public int getIsScan2() {
        return isScan2;
    }

    public List<Integer> getChargeStates() {
        return chargeStates;
    }

    public List<Float> getPrecursorMasses() {
        return precursorMasses;
    }

    @Override
    public String toString() {
        StringBuffer sb = new StringBuffer();

        sb.append("S\t").append(isScan1).append("\t").append(isScan2).append("\n");
        sb.append("Z\t").append(chargeStates).append("\t").append(precursorMasses).append("\n");

        if (false) {
            //do not print intensities by default
            for (int i = 0; i < size; ++i) {
                String mass = String.format("%.4f", masses.get(i));
                String intens = String.format("%.1f", intensities.get(i));
                sb.append(mass).append(" ").append(intens).append("\n");
            }
        }
        return sb.toString();

    }

    public String getScanType() {
        return scanType;
    }

    public void setScanType(String scanType) {
        this.scanType = scanType;
    }
    
    public void processHCDScan() {

        /*

        ArrayList<Peak> newpeaks = new ArrayList<Peak>(DEFAULTNUMPEAKS);
        ArrayList<Peak> toberemoved = new ArrayList();


        double maxmz = getMaxM2z();
        
        int peakindex = 0;    
        */
        
        float max = this.masses.get(this.masses.size()-1);
        int peakSize = masses.size(); 
        int peakindex = 0;     
     //   System.out.println("===" + max);
        List<Float> tempMasses = new ArrayList<Float>();
        List<Float> tempIntensities = new ArrayList<Float>();
        
        for(int i = 0; i*100 < max; i++) {
            ArrayList<MassIntensityModel> temp = new ArrayList(100);
            double uppermzlimit = (i+1)*100;
            
            while(peakindex < peakSize) {
                float mass = masses.get(peakindex);

              //  System.out.println("aa" + mass + " " + uppermzlimit);
                
                if(mass <  uppermzlimit) {
                    
                  //  System.out.println(mass);
                    temp.add(new MassIntensityModel(this.intensities.get(peakindex), peakindex));
                } else {
                    break;
                }
                
                peakindex++;
            }
          
//System.out.println("New Peak List:");
            //sort temp by intensity
            
          //  for(int ii=0;ii<temp.size();ii++) {
          //      System.out.println("------------\t" + ii + " " + temp.get(ii).getIntensity() + " " + temp.get(ii).getIndex());
          //  }
            
            Collections.sort(temp);
          //  for(int ii=0;ii<temp.size();ii++) {
          //      System.out.println("------------>>\t" + ii + " " + temp.get(ii).getIntensity() + " " + temp.get(ii).getIndex());
          //  }
          //  System.out.println("aaaaaa");
             
            int[] indexArr = new int[18];
            int tempCount=0;

            for(int j=0;j<=18;j++) {
            //for(int j = temp.size() - 18; j >=0 && j < temp.size(); j++) {
                if(j>=temp.size()) break;
                //System.out.println("==\t" + temp.get(j).getM2z());
                int tmpIndex = temp.get(j).getIndex();
                indexArr[tempCount++] = tmpIndex;                
//System.out.println(temp.get(j).getM2z() + " " + temp.get(j).getIntensity());
            }
            
            Arrays.sort(indexArr);
            
            for(int j=0;j<indexArr.length;j++) {
                
                if(indexArr[j]<=0) return;
                
            //    System.out.println(indexArr[j]);
                tempMasses.add(masses.get(indexArr[j]));
                tempIntensities.add(intensities.get(indexArr[j]));                 
            }
           
        } 
        
        for(int i=0;i<tempMasses.size();i++) {
            System.out.println(tempMasses.get(i));
        }
//System.out.println("Number of peaks before: " + peaks.size());
        this.masses = tempMasses;
        this.intensities = intensities;
        
    }
    
}
