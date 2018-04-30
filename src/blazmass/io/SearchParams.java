/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package blazmass.io;

/**
 *
 * @author rpark
 */

import blazmass.AssignMass;
import java.util.Date;
import java.util.*;
import blazmass.Enzyme;
import blazmass.dbindex.DBIndexer;
import blazmass.dbindex.DBIndexer.IndexType;
import blazmass.dbindex.Util;
        
public class SearchParams {
    
    private String program;
    private String parametersFile;
    private String parameters;
    private float peptideMassTolerance;
    
    private float fragmentIonTolerance;
    private int fragmentIonToleranceInt;
    private float fragmentIonToleranceBinScale;
    
    private int maxMissedCleavages;
    private int massTypeParent;
    private int massTypeFragment;   
    private int numPeptideOutputLnes;
    private int removePrecursorPeak;

    private boolean useMonoParent;
    private boolean useMonoFragment;
    private boolean diffSearch;
    private int enzymeOffset;
    private String databaseName;
    private String indexDatabaseName;
    
    private int maxNumDiffMod=3;
    private boolean variableTolerance=false;
    private float variablePeptideMassTolerance;
    private boolean usePPM =false;
    private float relativePeptideMassTolerance;
    private int isotopes=0;
    private float n15Enrichment=0.0f;
    private float matchPeakTolerance;
    private int maxInternalCleavageSites=3;
    private float hplusparent;
    private float hparent;
    private float oHparent;
    private float[] pdAAMassParent;
    private float[] pdAAMassFragment;    
    private int numIonSeriesUsed;
    private String enzymeBreakAA;
    private String enzymeNoBreakAA;    
    private float binWidth;
    private int[] ionSeries = new int[9];
    private int[] ionToUse;
    
    float diffMass1, diffMass2, diffMass3;
    String diffChar1, diffChar2, diffChar3;
    private boolean useEnzyme;
    private char[] enzymeArr;
    private static SearchParams sparams = null;
    
    private int neutralLossAions;
    private int neutralLossBions;
    private int neutralLossYions;
    private String enzymeName;
    private String enzymeResidues;
    private String enzymeCut;
    private String enzymeNocutResidues;
    private static float minPrecursorMass=500.0f;
    private static float maxPrecursorMass=6000.0f;
    private int minFragPeakNum=8;

    private boolean useIndex = true;
    private DBIndexer.IndexType indexType = DBIndexer.IndexType.INDEX_NORMAL; //default
    private boolean inMemoryIndex = false;
    private int indexFactor = 6;
    
    private Enzyme enzyme = new Enzyme();
    private static StringBuffer staticParams = new StringBuffer();
    private boolean highResolution=false;
    private int maxChargeState=6;
    float[] weightArr = new float[12];

    public SearchParams() {
        staticParams = new StringBuffer();
        sparams = null;

    }
    
    private List<ModResidue> modList = new ArrayList<ModResidue>();
    //private Set<Float> modShiftSet = new HashSet<Float>();
    private List<List<Double>> modGroupList = new ArrayList<List<Double>>();
    private boolean precursorHighResolution=true;
    
    private int scoreWin=10;
    private boolean neturalLossIsotope=false;
    
    private String mongodbServer;
    
    
    public static SearchParams getInstance() {
        if (sparams == null)
            sparams = new SearchParams();

        return sparams;
    }

    public IndexType getIndexType() {
        return indexType;
    }

    public void setIndexType(IndexType indexType) {
        this.indexType = indexType;
    }
    
    

    public boolean isUseIndex() {
        return useIndex;
    }

    public void setUseIndex(boolean useIndex) {
        this.useIndex = useIndex;
    }

    public boolean isInMemoryIndex() {
        return inMemoryIndex;
    }

    public void setInMemoryIndex(boolean inMemoryIndex) {
        this.inMemoryIndex = inMemoryIndex;
    }

    public int getIndexFactor() {
        return indexFactor;
    }

    public void setIndexFactor(int indexFactor) {
        this.indexFactor = indexFactor;
    }
    
    
    
    
    public void setMassTypeFragment(int massTypeFragment) {
        this.massTypeFragment = massTypeFragment;
    }

    public void setMassTypeParent(int massTypeParent) {
        this.massTypeParent = massTypeParent;
    }

    public void setMaxMissedCleavages(int maxMissedCleavages) {
        this.maxMissedCleavages = maxMissedCleavages;
    }

    public void setFragmentIonTolerance(float fragmentIonTolerance) {
            
        this.fragmentIonTolerance = fragmentIonTolerance;
        
        
        //System.out.println("============================" + fragmentIonTolerance);
        
        this.fragmentIonToleranceBinScale = 1000.0f/fragmentIonTolerance;
        //System.out.println("============================" + fragmentIonTolerance +  " " + this.fragmentIonToleranceBinScale);
       
        /*
        if(fragmentIonTolerance<100f)
            this.neturalLossIsotope=true;
        else
            this.neturalLossIsotope=false;
        */
        
     //   if(fragmentIonTolerance>300.0f || fragmentIonTolerance<=0.0f) 
     //       this.highResolution=false;
        
    }

    public void setPeptideMassTolerance(float peptideMassTolerance) {

        if(peptideMassTolerance>900)
                this.precursorHighResolution=false;

        this.peptideMassTolerance = peptideMassTolerance;
    }

    public void setParameters(String parameters) {
        this.parameters = parameters;
    }

    public void setParametersFile(String parametersFile) {
        this.parametersFile = parametersFile;
    }

    public void setProgram(String program) {
        this.program = program;
    }

    public String getProgram() {
        return program;
    }


    public String getParametersFile() {
        return parametersFile;
    }

    public String getParameters() {
        return parameters;
    }

    public float getPeptideMassTolerance() {
        return peptideMassTolerance;
    }

    public float getFragmentIonTolerance() {
        return fragmentIonTolerance;
    }

    public int getMaxMissedCleavages() {
        return maxMissedCleavages;
    }

    public int getMassTypeParent() {
        return massTypeParent;
    }

    public int getMassTypeFragment() {
        return massTypeFragment;
    }

    public int getNumPeptideOutputLnes() {
        return numPeptideOutputLnes;
    }

    public void setNumPeptideOutputLnes(int numPeptideOutputLnes) {
        this.numPeptideOutputLnes = numPeptideOutputLnes;
    }

    public int getRemovePrecursorPeak() {
        return removePrecursorPeak;
    }

    public void setRemovePrecursorPeak(int removePrecursorPeak) {
        this.removePrecursorPeak = removePrecursorPeak;
    }

    public float getBinWidth() {
        return binWidth;
    }

    public void setBinWidth(float binWidth) {
        this.binWidth = binWidth;
    }

    public String getDatabaseName() {
        return databaseName;
    }

    public void setDatabaseName(String databaseName) {
        this.databaseName = databaseName;
    }

    public boolean isDiffSearch() {
        return diffSearch;
    }

    public void setDiffSearch(boolean diffSearch) {
        this.diffSearch = diffSearch;
    }

    public String getEnzymeBreakAA() {
        return enzymeBreakAA;
    }

    public void setEnzymeBreakAA(String enzymeBreakAA) {
        this.enzymeArr = enzymeBreakAA.toCharArray();
        for(char ch:this.enzymeArr)        
            enzyme.addEnzyme(ch);        
                
        this.enzymeBreakAA = enzymeBreakAA;
    }

    public String getEnzymeNoBreakAA() {
        return enzymeNoBreakAA;
    }

    public void setEnzymeNoBreakAA(String enzymeNoBreakAA) {
        this.enzymeNoBreakAA = enzymeNoBreakAA;
    }

    public int getEnzymeOffset() {
        return enzymeOffset;
    }

    public void setEnzymeOffset(int enzymeOffset) {
        this.enzymeOffset = enzymeOffset;
    }

    public float getHparent() {
        return hparent;
    }

    public void setHparent(float hparent) {
        this.hparent = hparent;
    }

    public float getHplusparent() {
        return hplusparent;
    }

    public void setHplusparent(float hplusparent) {
        this.hplusparent = hplusparent;
    }

    public int[] getIonSeries() {
        return ionSeries;
    }

    public void setIonSeries(int[] ionSeries) {
        
        int count=0;
        List<Integer> l = new ArrayList<Integer>();
        
        for(int i=0;i<ionSeries.length;i++) {            
            
            if(ionSeries[i]>0) {
                //System.out.println("++" + i);
                l.add(i);
                
                /*
                switch (i) {                    
                    /*a* case 0: l.add(0); break;
                    /*b* case 1: l.add(1); break;
                    /*c* case 2: l.add(2); break;
                    /*x* case 6: l.add(6); break;
                    /*y* case 7: l.add(7); break;
                    /*z* case 8: l.add(8); break;
                }*/
            }
            
            count++;
        }

        int[] arr = new int[l.size()];
        count =0;
        for(Iterator<Integer> itr=l.iterator();itr.hasNext(); ) {
            arr[count++] = itr.next();            
        }
        
        this.ionToUse = arr;        
        this.ionSeries = ionSeries;
    }

    public int getIsotopes() {
        return isotopes;
    }

    public void setIsotopes(int isotopes) {
        this.isotopes = isotopes;
    }

    public float getMatchPeakTolerance() {
        return matchPeakTolerance;
    }

    public void setMatchPeakTolerance(float matchPeakTolerance) {
        this.matchPeakTolerance = matchPeakTolerance;
    }

    public int getMaxInternalCleavageSites() {
        return maxInternalCleavageSites;
    }

    public void setMaxInternalCleavageSites(int maxInternalCleavageSites) {
        this.maxInternalCleavageSites = maxInternalCleavageSites;
    }

    public int getMaxNumDiffMod() {
        return maxNumDiffMod;
    }

    public void setMaxNumDiffMod(int maxNumDiffMod) {
        this.maxNumDiffMod = maxNumDiffMod;
    }

    public float getN15Enrichment() {
        return n15Enrichment;
    }

    public void setN15Enrichment(float n15Enrichment) {
        this.n15Enrichment = n15Enrichment;
    }

    public int getNumIonSeriesUsed() {
        return numIonSeriesUsed;
    }

    public void setNumIonSeriesUsed(int numIonSeriesUsed) {
        this.numIonSeriesUsed = numIonSeriesUsed;
    }

    public float getoHparent() {
        return oHparent;
    }

    public void setoHparent(float oHparent) {
        this.oHparent = oHparent;
    }

    public float[] getPdAAMassFragment() {
        return pdAAMassFragment;
    }

    public void setPdAAMassFragment(float[] pdAAMassFragment) {
        this.pdAAMassFragment = pdAAMassFragment;
    }

    public float[] getPdAAMassParent() {
        return pdAAMassParent;
    }

    public void setPdAAMassParent(float[] pdAAMassParent) {
        this.pdAAMassParent = pdAAMassParent;
    }

    public float getRelativePeptideMassTolerance() {
        return relativePeptideMassTolerance;
    }

    public void setRelativePeptideMassTolerance(float relativePeptideMassTolerance) {
        this.relativePeptideMassTolerance = relativePeptideMassTolerance;
    }

    public boolean isUseMonoFragment() {
        return useMonoFragment;
    }

    public void setUseMonoFragment(boolean useMonoFragment) {
        this.useMonoFragment = useMonoFragment;
    }

    public boolean isUseMonoParent() {
        return useMonoParent;
    }

    public void setUseMonoParent(boolean useMonoParent) {
        this.useMonoParent = useMonoParent;
    }

    public float getVariablePeptideMassTolerance() {
        return variablePeptideMassTolerance;
    }

    public void setVariablePeptideMassTolerance(float variablePeptideMassTolerance) {
        this.variablePeptideMassTolerance = variablePeptideMassTolerance;
    }

    public boolean isVariableTolerance() {
        return variableTolerance;
    }

    public void setVariableTolerance(boolean variableTolerance) {
        this.variableTolerance = variableTolerance;
    }

    public boolean isUsePPM() {
        return usePPM;
    }

    public void setUsePPM(boolean usePPM) {
        this.usePPM = usePPM;
    }

    public String getDiffChar1() {
        return diffChar1;
    }

    public void setDiffChar1(String diffChar1) {
        this.diffChar1 = diffChar1;
    }

    public String getDiffChar2() {
        return diffChar2;
    }

    public void setDiffChar2(String diffChar2) {
        this.diffChar2 = diffChar2;
    }

    public String getDiffChar3() {
        return diffChar3;
    }

    public void setDiffChar3(String diffChar3) {
        this.diffChar3 = diffChar3;
    }

    public float getDiffMass1() {
        return diffMass1;
    }

    public void setDiffMass1(float diffMass1) {
        this.diffMass1 = diffMass1;
    }

    public float getDiffMass2() {
        return diffMass2;
    }

    public void setDiffMass2(float diffMass2) {
        this.diffMass2 = diffMass2;
    }

    public float getDiffMass3() {
        return diffMass3;
    }

    public void setDiffMass3(float diffMass3) {
        this.diffMass3 = diffMass3;
    }

    public boolean isUseEnzyme() {
        return useEnzyme;
    }

    public void setUseEnzyme(boolean useEnzyme) {
        this.useEnzyme = useEnzyme;
    }

    public char[] getEnzymeArr() {
        return enzymeArr;
    }

    public void setEnzymeArr(char[] enzymeArr) {
        this.enzymeArr = enzymeArr;
    }


    public int getNeutralLossAions() {
        return neutralLossAions;
    }

    public void setNeutralLossAions(int neutralLossAions) {
        this.neutralLossAions = neutralLossAions;
    }

    public int getNeutralLossBions() {
        return neutralLossBions;
    }

    public void setNeutralLossBions(int neutralLossBions) {
        this.neutralLossBions = neutralLossBions;
    }

    public int getNeutralLossYions() {
        return neutralLossYions;
    }

    public void setNeutralLossYions(int neutralLossYions) {
        this.neutralLossYions = neutralLossYions;
    }

    public int[] getIonToUse() {
        return ionToUse;
    }

    public void setIonToUse(int[] ionToUse) {
        this.ionToUse = ionToUse;
    }

    public static SearchParams getSparams() {
        return sparams;
    }

    public static void setSparams(SearchParams sparams) {
        SearchParams.sparams = sparams;
    }

    public void addModResidue(ModResidue r) {
        modList.add(r);
    }
    
    public List<ModResidue> getModList() {
        return modList;
    }

    public void setModList(List<ModResidue> modList) {
        this.modList = modList;
    }

    public List<List<Double>> getModGroupList() {
        return modGroupList;
    }

    public void setModGroupList(List<List<Double>> modGroupList) {
        this.modGroupList = modGroupList;
    }

    public void addModGroupList(List<Double> modGroup) {
        this.modGroupList.add(modGroup);
    }

    public String getIndexDatabaseName() {
        return indexDatabaseName;
    }

    public void setIndexDatabaseName(String indexDatabaseName) {
        this.indexDatabaseName = indexDatabaseName;
    }

    public String getEnzymeName() {
        return enzymeName;
    }

    public void setEnzymeName(String enzymeName) {
        this.enzymeName = enzymeName;
    }
    
    public String getEnzymeCut() {
        return enzymeCut;
    }

    public void setEnzymeCut(String enzymeCut) {
        this.enzymeCut = enzymeCut;
    }

    public String getEnzymeResidues() {
        return enzymeResidues;
    }

    public void setEnzymeResidues(String enzymeResidues) {
        this.enzymeResidues = enzymeResidues;
    }

    public String getEnzymeNocutResidues() {
        return enzymeNocutResidues;
    }

    public void setEnzymeNocutResidues(String enzymeNocutResidues) {
        this.enzymeNocutResidues = enzymeNocutResidues;
    }

    public String getFullFileNameWithNoIndex() {

        String uniqueIndexName = databaseName + "_";

        //generate a unique string based on current params that affect the index
        final StringBuilder uniqueParams = new StringBuilder();
        //uniqueParams.append(sparam.getEnzyme().toString());
        //uniqueParams.append(sparam.getEnzymeNumber());
        uniqueParams.append(getEnzymeOffset());
        uniqueParams.append("_").append(getEnzymeResidues());
        if(null == getEnzymeNocutResidues())
            uniqueParams.append("_NA");
        else
            uniqueParams.append("_").append(getEnzymeNocutResidues());

        //uniqueParams.append(" ").append(getIndexFactor() ) ;

        uniqueParams.append("_cleav_");
        uniqueParams.append(getMaxInternalCleavageSites());
        uniqueParams.append("_");
        uniqueParams.append(getMaxMissedCleavages());
        uniqueParams.append("_");
        String staticParamStr = SearchParams.getStaticParams().toString().trim().replaceAll("\\.", "-");
        uniqueParams.append("static").append("_").append(staticParamStr.trim());

//        System.out.println("--->>" + uniqueParams.toString());

        return uniqueIndexName + uniqueParams.toString();
    }

    public String getFullIndexFileName() {    
        
        String uniqueIndexName = databaseName + "_";

        //generate a unique string based on current params that affect the index
        final StringBuilder uniqueParams = new StringBuilder();
        //uniqueParams.append(sparam.getEnzyme().toString());
        //uniqueParams.append(sparam.getEnzymeNumber());                
        uniqueParams.append(getEnzymeOffset());
        uniqueParams.append(" ").append(getEnzymeResidues());
        uniqueParams.append(" ").append(getEnzymeNocutResidues());
        
        //uniqueParams.append(" ").append(getIndexFactor() ) ;

        uniqueParams.append(", Cleav: ");
        uniqueParams.append(getMaxInternalCleavageSites());
        uniqueParams.append(getMaxMissedCleavages());

        uniqueParams.append(", Static: ").append(SearchParams.getStaticParams());
                
        
        /*uniqueParams.append(getMaxNumDiffMod());
        uniqueParams.append("\nMods:");
        for (final ModResidue mod : getModList()) {
            uniqueParams.append(mod.toString()).append(" ");
        }
        uniqueParams.append("\nMods groups:");
        for (final List<Float> modGroupList : getModGroupList()) {
            for (final Float f : modGroupList) {
                uniqueParams.append(f).append(" ");
            }
        }*/

//System.out.println("===" + uniqueParams.toString());

        final String uniqueParamsStr = uniqueParams.toString();

        //logger.log(Level.INFO, "Unique params: " + uniqueParamsStr);

        String uniqueParamsStrHash = Util.getMd5(uniqueParamsStr);
        System.out.println("param===========" + uniqueParamsStr + "\t" + uniqueParamsStrHash);

        return uniqueIndexName + uniqueParamsStrHash;
    }


    public static void addStaticParam(char ch, float f) {
        //    public static void addMass(int i, float mass) {
        SearchParams.staticParams.append(ch).append(f);
    }

    public static void addStaticParam(String str, float f) {
        //    public static void addMass(int i, float mass) {

        SearchParams.staticParams.append(str.trim()).append(f);
    }

    public static StringBuffer getStaticParams() {        

        return staticParams;
    }

    public static void setStaticParams(StringBuffer staticParams) {
        SearchParams.staticParams = staticParams;
    }    

    public float getMinPrecursorMass() {
        return minPrecursorMass;
    }

    public void setMinPrecursorMass(float minPrecursorMass) {
        this.minPrecursorMass = minPrecursorMass;
    }

    public float getMaxPrecursorMass() {
        return maxPrecursorMass;
    }

    public void setMaxPrecursorMass(float maxPrecursorMass) {
        this.maxPrecursorMass = maxPrecursorMass;
    }

    public int getMinFragPeakNum() {
        return minFragPeakNum;
    }

    public void setMinFragPeakNum(int minFragPeakNum) {
        this.minFragPeakNum = minFragPeakNum;
    }

    public boolean isHighResolution() {
        return highResolution;
    }

    public void setHighResolution(boolean highResolution) {
        this.highResolution = highResolution;
    }

    public int getFragmentIonToleranceInt() {
        return fragmentIonToleranceInt;
    }

    public void setFragmentIonToleranceInt(int fragmentIonToleranceInt) {
                
        this.fragmentIonToleranceInt = fragmentIonToleranceInt;
    }    

    public int getMaxChargeState() {
        return maxChargeState;
    }

    public void setMaxChargeState(int maxChargeState) {
        this.maxChargeState = maxChargeState;
    }

    public float getFragmentIonToleranceBinScale() {
        return fragmentIonToleranceBinScale;
    }

    public void setFragmentIonToleranceBinScale(float fragmentIonToleranceBinScale) {
        this.fragmentIonToleranceBinScale = fragmentIonToleranceBinScale;
    }

    public boolean isPrecursorHighResolution() {
        return precursorHighResolution;
    }

    public void setPrecursorHighResolution(boolean precursorHighResolution) {
        this.precursorHighResolution = precursorHighResolution;
    }

    public int getScoreWin() {
        return scoreWin;
    }

    public void setScoreWin(int scoreWin) {
        this.scoreWin = scoreWin;
    }

    public boolean isNeturalLossIsotope() {
        return neturalLossIsotope;
    }

    public void setNeturalLossIsotope(boolean neturalLossIsotope) {
        this.neturalLossIsotope = neturalLossIsotope;
    }

    public float[] getWeightArr() {
        return weightArr;
    }

    public void setWeightArr(float[] weightArr) {
        this.weightArr = weightArr;
    }

    public Enzyme getEnzyme() {
        return enzyme;
    }

    public void setEnzyme(Enzyme enzyme) {
        this.enzyme = enzyme;
    }

    public String getMongodbServer() {
        return mongodbServer;
    }

    public void setMongodbServer(String mongodbServer) {
        this.mongodbServer = mongodbServer;
    }
    

    public float getWeight(int i) {
        return this.weightArr[i];
    }
    
    
    
}

