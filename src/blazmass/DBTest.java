/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package blazmass;

//import blazmass.dbindex.DBIndexStore;
//import blazmass.dbindex.DBIndexStoreSQLiteMult;
import blazmass.dbindex.DBIndexer;
import blazmass.dbindex.DBIndexerException;
import blazmass.io.FastaReader;
import blazmass.io.SearchParamReader;
import blazmass.io.SearchParams;
//import gnu.trove.list.array.TIntArrayList;
//import gnu.trove.map.hash.TIntObjectHashMap;
//import gnu.trove.map.hash.TObjectIntHashMap;
//
//
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.logging.Level;
//import org.apache.commons.collections4.Trie;
//import org.mapdb.LongHashMap;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.TIntObjectMap;
/**
 *
 * @author rpark
 */
public class DBTest {
    
    
  //  private DBIndexStore indexStore;
    private static SearchParams sparam;
    //private HashMap<> dbMap = new LongHashMap(6000000);
    //private List<StringBuffer> totalList = new ArrayList<>() ;
    private TIntObjectMap dbMap = new TIntObjectHashMap();
    //private Object[] pepArr = new Object[6000001];
  //private StringBuffer[] pepArr = new StringBuffer[6000001];
    public static void main(String[] args) throws Exception {
        DBTest test = new DBTest();
        
        test.run();
        
    }
    
    public void run() throws Exception {
        //indexStore = new DBIndexStoreSQLiteMult();
        //String path = "/home/rampuria/data/prolucid";
        String path = "/home/rpark/test_data/prolucid";
        String paramFile = "search.xml";
                  SearchParamReader preader;

            preader = new SearchParamReader(path, paramFile);

        sparam = preader.getParam();


        //indexer in indexing mode
        final DBIndexer indexer = new DBIndexer(sparam, DBIndexer.IndexerMode.INDEX);
       
        double indexedProteins=0;
        int count =0;
        double percent=0.0;
        FileInputStream fis = new FileInputStream(sparam.getDatabaseName());
                    for (Iterator<blazmass.io.Fasta> itr = FastaReader.getFastas(fis); itr.hasNext();) {
                //int size = Iterators.size(itr);
                blazmass.io.Fasta fasta = itr.next();
               // protCache.addProtein(fasta.getSequestLikeAccession(), fasta.getSequence());
                cutSeq(fasta);
                      //  System.out.println("Percentage completed : "in);
//                System.out.print("Printing the last buffer....");
//                indexStore.lastBuffertoDatabase();
                      

                ++indexedProteins;
                percent = (indexedProteins/47022)*100;
                
                System.out.println("Percentage completed: "+percent);
                /*if (statusWriter != null) {
                    statusWriter.append(totalProteins).append("\t").append(Integer.toString(indexedProteins)).append("\n");
                    if (indexedProteins % 100 == 0) {
                        statusWriter.flush();
                    }
                }
                */
            }  
            //indexer.init();
          //  indexer.run();
          System.out.print("");
       
    }
    
        private void cutSeq(final blazmass.io.Fasta fasta) throws IOException {
        final String protAccession = fasta.getSequestLikeAccession();
        final String protSeq = fasta.getSequence();
       
        cutSeq(protAccession, protSeq);
        

    }

    /**
     * Cut a Fasta protein sequence according to params spec and index the
     * protein and generated sequences
     *
     * Should be called once per unique Fasta sequence
     *
     * @param fasta fasta to index
     * @throws IOException
     */
    private void cutSeq(final String protAccession, final String protSeq) throws IOException {

        //Enzyme enz = sparam.getEnzyme();
        //System.out.print(".");
        final int length = protSeq.length();

        final int MAX_SEQ_LENGTH = 1000;
        //AssignMass aMass = AssignMass.getInstance(true);

        final char[] pepSeq = new char[MAX_SEQ_LENGTH]; //max seq length
        int curSeqI = 0;

        final int maxMissedCleavages = sparam.getMaxMissedCleavages();
        int maxIntCleavage = sparam.getMaxInternalCleavageSites();

          //  long proteinId = indexStore.addProteinDef(++protNum, protAccession, protSeq);
//	    System.out.println(fasta.getSequestLikeAccession());
            //System.out.println(fasta.getDefline());
         //   DBPeptide pep = null;
            for (int start = 0; start < length; ++start) {
                int end = start;

                //clear the preallocated seq byte array
                //Arrays.fill(seq, 0, curSeqI > 0?curSeqI-1:0, (byte) 0); //no need, we copy up to curSeqI nowu
                curSeqI = 0;

                //float precMass = Constants.H2O_PROTON_SCALED_DOWN;
                float precMass = Constants.H2O_PROTON;
                precMass += AssignMass.getcTerm();
                precMass += AssignMass.getnTerm();
                

                // System.out.println("===>>" + precMass + "\t" + Constants.MAX_PRECURSOR);
                //System.out.println("==" + j + " " + length + " " + (j < length));

                //int testC=0;
                int pepSize = 0;

                int intMisCleavageCount = -1;

               // TIntObjectHashMap<DBPeptide> dbMap = new TIntObjectHashMap<>();
                  
                //while (precMass <= Constants.MAX_PRECURSOR_MASS && end < length) {
                if(precMass<3000||precMass>3600) continue;

                while (precMass <= sparam.getMaxPrecursorMass() && end < length) 
                {
                    pepSize++;

                    final char curIon = protSeq.charAt(end);
                    pepSeq[curSeqI++] = curIon;
                    precMass += AssignMass.getMass(curIon);

                    if (Enzyme.isEnzyme(protSeq.charAt(end))) {
                        intMisCleavageCount++;
                    }

                    final int cleavageStatus = Enzyme.checkCleavage(protSeq, start, end, sparam.getEnzymeNocutResidues());

                    //System.out.println("---" + String.valueOf(Arrays.copyOf(pepSeq, curSeqI)) + "\t" + intMisCleavageCount + "\t" + maxIntCleavage);

                    if (intMisCleavageCount > maxIntCleavage) {
                        break;
                    }

                    
                    //if (precMass > Constants.MAX_PRECURSOR_MASS) {
                    if (precMass > sparam.getMaxPrecursorMass()) {
                        break;
                    }
                    
                    if (pepSize >= Constants.MIN_PEP_LENGTH && precMass >= sparam.getMinPrecursorMass() ) { //Constants.MIN_PRECURSOR ) {
                        if (cleavageStatus >= maxMissedCleavages) {
                        //if (cleavageStatus == 2) {
                            //qualifies based on params

                            final String peptideSeqString = String.valueOf(Arrays.copyOf(pepSeq, curSeqI));

                            //check if index will accept it
                          //  final DBIndexStore.FilterResult filterResult = indexStore.filterSequence(precMass, peptideSeqString);

                            //if (filterResult.equals(DBIndexStore.FilterResult.SKIP_PROTEIN_START) ) {
                                //bail out earlier as we are no longer interested in this protein starting at start
                             //   break; //move to new start position
                           // }
                            //else if (filterResult.equals(DBIndexStore.FilterResult.INCLUDE) ) {
                                final int resLeftI = start >= Constants.MAX_INDEX_RESIDUE_LEN ? start - Constants.MAX_INDEX_RESIDUE_LEN : 0;
                                final int resLeftLen = Math.min(Constants.MAX_INDEX_RESIDUE_LEN, start);
                                StringBuilder sbLeft = new StringBuilder(Constants.MAX_INDEX_RESIDUE_LEN);
                                for (int ii = 0; ii < resLeftLen; ++ii) {
                                    sbLeft.append(protSeq.charAt(ii + resLeftI));
                                }
                                final int resRightI = end + 1;
                                final int resRightLen = Math.min(Constants.MAX_INDEX_RESIDUE_LEN, length - end - 1);
                                StringBuilder sbRight = new StringBuilder(Constants.MAX_INDEX_RESIDUE_LEN);
                                if (resRightI < length) {
                                    for (int jj = 0; jj < resRightLen; ++jj) {
                                        sbRight.append(protSeq.charAt(jj + resRightI));
                                    }
                                }


                                //add -- markers to fill Constants.MAX_INDEX_RESIDUE_LEN length
                                final int lLen = sbLeft.length();
                                for (int c = 0; c < Constants.MAX_INDEX_RESIDUE_LEN - lLen; ++c) {
                                    sbLeft.insert(0, '-');
                                }
                                final int rLen = sbRight.length();
                                for (int c = 0; c < Constants.MAX_INDEX_RESIDUE_LEN - rLen; ++c) {
                                    sbRight.append('-');
                                }

                                final String resLeft = sbLeft.toString();
                                final String resRight = sbRight.toString();
                                
                                
                               int key = (int)(precMass*1000);
                               
                              // DBPeptide tmp = totalList[key];
                              // DBPeptide tmp = dbMap.get(key);
                            //System.out.println(key);
                               Object obj = this.dbMap.get(key);



                               if(obj != null){
                                  /*  tmp.addStarttoList(start);
                                    tmp.addcurseqItoList(curSeqI);
                                    tmp.addresLefttoList(resLeft);
                                    tmp.addresRighttoList(resRight);
                                    tmp.addSeqtoList(peptideSeqString);*/
                                   TByteArrayList byteArrayList = (TByteArrayList)obj;

                                    //StringBuffer sb = new StringBuffer(1000);
                                   byteArrayList.add(";".getBytes());
                                   byteArrayList.add(String.valueOf(start).getBytes());
                                   byteArrayList.add(",".getBytes());
                                   byteArrayList.add(String.valueOf(curSeqI).getBytes());
                                   byteArrayList.add(resLeft.getBytes());
                                   byteArrayList.add(resRight.getBytes());
                                   byteArrayList.add(peptideSeqString.getBytes());

                                   /*
                                   sb.append(";").append(start).append(",")
                                            .append(curSeqI).append(",")
                                            .append(resLeft)
                                            .append(resRight)
                                           .append(peptideSeqString).append(",");
*/
                                   //byteArrayList.add(sb.toString().getBytes());
                                   // totalList[key]=sb;
                                   //dbMap.put(key,sb);
                                   
                                }
                                else {
                                  /* pep = new DBPeptide();
                                    pep.addStarttoList(start);
                                    pep.addcurseqItoList(curSeqI);
                                    pep.addresLefttoList(resLeft);
                                    pep.addresRighttoList(resRight);
                                    pep.addSeqtoList(peptideSeqString);
                                    totalList[key]=pep;*/
                                    
                                    //String str =start+","+curSeqI+","+resLeft+","+resRight+","+peptideSeqString;
                                   TByteArrayList byteArrayList = new TByteArrayList();


                                   //StringBuffer sb = new StringBuffer(1000);

                                   byteArrayList.add(String.valueOf(start).getBytes());
                                   byteArrayList.add(",".getBytes());
                                   byteArrayList.add(String.valueOf(curSeqI).getBytes());
                                   byteArrayList.add(resLeft.getBytes());
                                   byteArrayList.add(resRight.getBytes());
                                   byteArrayList.add(peptideSeqString.getBytes());

                                   dbMap.put(key, byteArrayList);
                                   /*
                                   sb.append(start).append(",")
                                           .append(curSeqI).append(",")
                                           .append(resLeft)
                                           .append(resRight)
                                           .append(peptideSeqString).append(",");



                                   this.pepArr[key] = sb;
                                   */

                                   //dbMap.put(key, pep);
                                   //dbMap.put(key, sb);
                            //   }
                               // System.out.println(precMass + peptideSeqString);

                               // if(peptideSeqString.contains("WHLKTEIES"))

                              //  System.out.println(precMass + peptideSeqString);
                                //indexStore.addSequence(precMass, start, curSeqI, peptideSeqString, resLeft, resRight, proteinId);
                               //System. out.println((int)(precMass*1000) + "\t" + String.valueOf(Arrays.copyOf(pepSeq, curSeqI)));
                            } //end if add sequence
                        }
                    }


                    if (intMisCleavageCount == maxIntCleavage) {
                        break;
                    }

                    ++end;

                }
                
//
            }

/*
        Runtime runtime = Runtime.getRuntime();

        //NumberFormat format = NumberFormat.getInstance();


        long maxMemory = runtime.maxMemory();
        long allocatedMemory = runtime.totalMemory();
        long freeMemory = runtime.freeMemory();

        System.out.println("free memory: " + freeMemory / 1024);
        System.out.println("allocated memory: " + allocatedMemory / 1024);
        System.out.println("max memory: " + maxMemory / 1024);
        System.out.println("total free memory: " + (freeMemory + (maxMemory - allocatedMemory)) / 1024);
*/

    }
}
