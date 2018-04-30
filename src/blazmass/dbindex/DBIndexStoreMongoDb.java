package blazmass.dbindex;

import blazmass.io.FastaDefReader;
import blazmass.io.FastaReader;
import blazmass.io.SearchParams;
import blazmass.util.ProcessUtil;
import com.mongodb.BasicDBObject;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;


import com.mongodb.MongoClient;
import com.mongodb.DB;
import com.mongodb.DBCollection;
import com.mongodb.DBCursor;
import com.mongodb.DBObject;
import java.io.FileWriter;
import java.net.UnknownHostException;
import java.util.Collection;
import java.util.Iterator;

/**
 *
 * Index store implementation based on mongo db engine. To support large indexes
 * and search of large databases.
 *
 * @author Adam
 */

/*
 * @author HARSHIL SHAH
 * 
 * It might cause problem when having multiple connection in mongoDB with the DBNAME 
 * as in the constructor i have get the database .. I am not sure with this,.
 * 
 * 
 */
public class DBIndexStoreMongoDb implements DBIndexStore {

    //Mongo
    private DB db;
    private String MongoServerName;
    private MongoClient mongoClient;
    private DBCollection sequences;
    private ProteinCache protCache;

    private enum FIELDS {

        MASS, PROT_ID, SEQUENCE, LR, RR, LEN, OFFSET, BUCKET
    };
    private   String COLLECTION = "blazmass_sequences";
//     private static final String COLLECTION = "test31";
   
    private String dbName;
    private final String INDEX_DIR_SUFFIX = ".idx2";
    private String indexDir;
    //indexing
    //searching
    private static final Logger logger = Logger.getLogger(DBIndexStoreMongoDb.class.getName());
    private long indexed = 0;
    private long totalDatabaseEntry=0;
    private long counter=1,maxBufferLimit=100000;
    private boolean inited = false;
    private static int DB_PORT = 27017;
    private final ProcessUtil mongoDbProcess = new ProcessUtil();
    private FileWriter mongoLogWriter;
    private Map<Integer , ArrayList<BufferSequences> > hashBuffer = new HashMap<>() ;
    private BufferSequences buffersequence =null;
    private ArrayList<BufferSequences> listofBufferSequencese = null;
    private List massList = new ArrayList();
    private FastaDefReader fastadefreader= new FastaDefReader();
    
    DBIndexStoreMongoDb()
    {}
    DBIndexStoreMongoDb (SearchParams sparam)
    {
        
        File f= new File(sparam.getDatabaseName());
        String filename= f.getName();
          this.dbName=filename.replaceAll(".fasta", "");
          f=new File(sparam.getFullIndexFileName());     
          filename= f.getName();
         this.COLLECTION =filename.substring(filename.indexOf("fasta")+6, filename.length());
      //  this.dbName= sparam.getFullIndexFileName();
        this.MongoServerName = sparam.getMongodbServer();
             //     this.MongoServerName = "localhost";
//        System.out.println("The file name is :" + sparam.getDatabaseName());
        fastadefreader.setDefs(sparam.getDatabaseName());
        setProteinCache(fastadefreader.getproteinCache());
    }

    //creates mongdo db dir for the index 
    //relative to current dir
    //and starts mongod
    //mongod must be either in ./bin relative to current dir
    //or in PATH
    private void startDB(String databaseDir) throws DBIndexStoreException {
/*
        //create db folder if does not exist, relative to current dir

        File dbDir = new File(databaseDir);

        if (!dbDir.exists()) {
            logger.log(Level.INFO, "Creating db folder: " + dbDir.getAbsolutePath());
            boolean created = dbDir.mkdirs();
            if (!created) {
                throw new DBIndexStoreException("Error creating db folder: " + dbDir.getAbsolutePath());
            }
        }

        String mongoProcName = "mongod";
        if (System.getProperty("os.name").toLowerCase().contains("win")) {
            mongoProcName += ".exe";
        }
        
        else {
            try {
                //ensure previous instances are out (for now linux only)
                //run pkill -9 mongod 
                logger.log(Level.INFO, "Try to stop previous mongod instance");
                String [] cmd  = { "pkill", "-9", mongoProcName};
                Runtime.getRuntime().exec(cmd);
            } catch (IOException ex) {
                logger.log(Level.WARNING, "Could stop previous mongod instance", ex);
            }
        }

        //locate mongod, it could be in 3 diferent locations, finish with trying PATH
        File mongodFile = new File("../bin/" + mongoProcName);
        boolean usePath = false;
        if (mongodFile.exists() && mongodFile.canExecute()) {
            logger.log(Level.INFO, "Found mongod at: " + mongodFile.getAbsolutePath());
        } else {
            mongodFile = new File("./bin/" + mongoProcName);
            if (mongodFile.exists() && mongodFile.canExecute()) {
                logger.log(Level.INFO, "Found mongod at: " + mongodFile.getAbsolutePath());
            } else {
                usePath = true;
            }
            logger.log(Level.INFO, "Will try mongod from PATH" + mongodFile.getAbsolutePath());
        }

        final String mongoPath = usePath ? mongoProcName : mongodFile.getAbsolutePath();


        try {
            mongoLogWriter = new FileWriter("mongod.log");
            mongoDbProcess.execute(mongoLogWriter, false, mongoPath,
                    "--port", Integer.toString(DB_PORT),
                    "--bind_ip", "localhost",
                    "--dbpath", dbDir.getAbsolutePath(),
                    //  "--directoryperdb",
                    "--cpu",
                    "--maxConns", "12",
                    "--noauth" //"--nojournal" //will cause corrutpion if not shutdown cleanly
                    //"--sysinfo"
                    );
            Thread.sleep(2000);
        } catch (IOException ex) {
            logger.log(Level.SEVERE, "Cannot start mongodb instance, is mongod in your PATH?", ex);
        } catch (InterruptedException ex) {
            logger.log(Level.SEVERE, "Cannot start mongodb instance, is mongod in your PATH?", ex);
        }
*/
    }

    private void stopDB() {
        //db.shutdownServer({timeoutSecs : 5})



        //   if (db != null) {
        //  db.command("{shutdown : 1, timeoutSecs : 5}");
        //     db = null;
        //}

        if (mongoClient != null) {

            mongoClient.close();
            mongoClient = null;
        }

        mongoDbProcess.stop();
    }

    @Override
    protected void finalize() throws Throwable {

        try {
            //in case not yet stopped
            stopDB();
        } finally {
            super.finalize();
        }
    }

    @Override
    public void init(String databaseID) throws DBIndexStoreException {
        this.indexDir = databaseID + INDEX_DIR_SUFFIX;

      //  startDB(indexDir);

        try {
           // mongoClient = new MongoClient("localhost", DB_PORT);
//             mongoClient = new MongoClient("thrall.scripps.edu", DB_PORT);
              mongoClient = new MongoClient(MongoServerName, DB_PORT);
//             System.out.println("Connected to the mongoDB database...");
              System.out.println("Connected to MongoD_client:"+MongoServerName+"\t"+ "databse\t" + dbName + "\tcollection\t" + COLLECTION);
        //    fastadefreader.setDefs(this.);
             
             
        } catch (UnknownHostException ex) {
            logger.log(Level.SEVERE, "Cannot connect to MongoDb instance at port: " + DB_PORT, ex);
            throw new DBIndexStoreException("Cannot connect to MongoDb instance at port: " + DB_PORT, ex);
        }
        db = mongoClient.getDB(dbName);

        sequences = db.getCollection(COLLECTION);

        if (sequences != null) {
            inited = true;
        }
                  sequences.createIndex(new BasicDBObject("MASS", 1));

    }

    @Override
    public void startAddSeq() throws DBIndexStoreException {
        logger.log(Level.INFO, "Starting indexing");


        //db.command("{\"beginTransaction\", \"isolation\": \"mvcc\"}");
        // sequences.beginTransaction();


    }

    @Override
    public void stopAddSeq() throws DBIndexStoreException {

        //db.command("commitTransaction");
        sequences.ensureIndex(FIELDS.MASS.toString());
//        
//        if (sequences != null) {
//            logger.log(Level.INFO, "Starting index optimization");
//            //sequences.ensureNumderIndex(FIELDS.MASS.toString());
//
//            sequences.ensureIndex(FIELDS.MASS.toString());
//
//            logger.log(Level.INFO, "Starting index commit");
//            //sequences.commitTransaction();
//
//            logger.log(Level.INFO, "Starting index sync to disk");
//
//           // ejdb.sync();
//            //ejdb.close();
//
//            
//            
//            logger.log(Level.INFO, "Done indexing");
//        }

        this.stopDB();
    }

    @Override
    public boolean indexExists() throws DBIndexStoreException {
        if (!inited) {
            return false;
        }

        File dir = new File(indexDir);
        if (!dir.exists() && !dir.isDirectory()) {
            return false;
        }

        return this.getNumberSequences() > 0;
    }

    @Override
    public FilterResult filterSequence(float precMass, String sequence) {
        //no filtering, we add every sequence
        return FilterResult.INCLUDE;
    }

    public void lastBuffertoDatabase()
    {
        addData adddata= new addData();
        adddata.BuffertoDatabase();
            massList.clear();
            hashBuffer.clear();
            sequences.createIndex(new BasicDBObject("MASS", 1));
    }
    
    
    @Override
    public void addSequence(float precMass, int sequenceOffset, int sequenceLen,
            String sequence, String resLeft, String resRight, long proteinId)
            throws DBIndexStoreException {
        
        if (indexed % 1000 == 0) {
//            System.out.println(indexed);
        }
        addData adddata = new addData();
        String seq;
//        if(sequence.equals("NPRTSPEKKATQTI"))
//            System.out.println(precMass + "NPRTSPEKKATQTI");
        int mass = (int)(precMass*1000 + 0.5);
        indexed++;
        if(!((counter % maxBufferLimit) == 0))
//        if((counter <maxBufferLimit))
        {
//            listofBufferSequencese.get
//            buffersequence = new BufferSequences(sequence, sequenceLen,resLeft,resRight,proteinId);
            
            adddata.datatoBuffer(mass,sequence,resLeft,resRight,sequenceLen,proteinId);
            counter++;
           // hashBuffer.put(mass, null);
        }
        else
        {
            counter=1;
            
            adddata.BuffertoDatabase();
            massList.clear();
            hashBuffer.clear();
            sequences.createIndex(new BasicDBObject("MASS", 1));

//            int bucket = (int) precMass / 1000;
//   //MASS, PROT_ID, SEQUENCE, LR, RR, LEN, OFFSET, BUCKET
//           BasicDBObject newSequence = new BasicDBObject()
//                   .append(FIELDS.MASS.toString(), mass)
//                   .append(FIELDS.PROT_ID.toString(), Long.toString(proteinId))
//                   .append(FIELDS.OFFSET.toString(), sequenceOffset)
//                   .append(FIELDS.SEQUENCE.toString(),sequence)
//                   .append(FIELDS.LR.toString(),resLeft)
//                   .append(FIELDS.RR.toString(),resRight)
//                   .append(FIELDS.LEN.toString(), sequenceLen);
//           ////.append(FIELDS.SEQUENCE.toString(), sequence) do not store 
//           //.append(FIELDS.LR.toString(), resLeft)
//           //.append(FIELDS.RR.toString(), resRight)
//           //.append(FIELDS.BUCKET.toString(), bucket);
//
//              sequences.insert(newSequence);
              
        }

    }
    
    public void test(float precMass,float tolerance, DBCursor cursor) throws DBIndexStoreException
        {
 /*       ResultsCollector collector = new ResultsCollector();
        List<IndexedSequence> allSeqs = null;
        try {
            collector.collect(cursor);
            allSeqs = collector.getSequences();
            //note: massRangeResults might contain duplicate sequences
            //we should either filter out, or faster: add step during indexing to remove duplicates
            //and merge sequences, and preserve protein id information.  Check sqlite index implementation for example.
           System.out.println("Printing after removing all the repeating occurances");
            for(IndexedSequence s : allSeqs)
                System.out.println(s.toString());
            
            
            System.out.println(allSeqs.size());
            //System.out.println(allSeqs);
        } catch (IOException ex) {
            logger.log(Level.SEVERE, "Cannot collect search results ", ex);

        } finally {
            cursor.close();
        }


*/
    }
    @Override
    public List<IndexedSequence> getSequences(float precMass, float tolerance) throws DBIndexStoreException {

        List<IndexedSequence> allSeqs = null;
/*

        float massLow = precMass - tolerance;
        if (massLow < 0f) {
            massLow = 0f;
        }
        float massHigh = precMass + tolerance;


        BasicDBObject query = new BasicDBObject(FIELDS.MASS.toString(), new BasicDBObject("$gte", massLow).
                append("$lte", massHigh));
        DBCursor cursor = sequences.find(query);

        ResultsCollector collector = new ResultsCollector();

        try {
            collector.collect(cursor);
            allSeqs = collector.getSequences();
            //note: massRangeResults might contain duplicate sequences
            //we should either filter out, or faster: add step during indexing to remove duplicates
            //and merge sequences, and preserve protein id information.  Check sqlite index implementation for example.

        } catch (IOException ex) {
            logger.log(Level.SEVERE, "Cannot collect search results ", ex);

        } finally {
            cursor.close();
        }

 */
        return allSeqs;
       
    }

    @Override
    public List<IndexedSequence> getSequences(List<MassRange> ranges) throws DBIndexStoreException {

        List<IndexedSequence> allSeqs = new LinkedList<IndexedSequence>();
/*        MongoClient mongoClient1 = null;
            try {
                mongoClient1 = new MongoClient( MongoServerName, DB_PORT);
            } catch (UnknownHostException ex) {
                Logger.getLogger(DBIndexStoreMongoDb.class.getName()).log(Level.SEVERE, null, ex);
            }
                  System.out.println("Connected to MongoD_client1");
          //  DB db = mongoClient1.getDB(dbName);
          //  DBCollection dbcollection = db.getCollection(COLLECTION);
 */         
            DBCollection dbcollection = sequences;
            System.out.println("Connected to MongoD_client1:"+MongoServerName+"\t"+ "databse\t" + dbName + "collection\t" + COLLECTION);

        //convert mass ranges to intverals and merge overlapping mass ranges
        ArrayList<Interval> intervals = new ArrayList<Interval>();
        for (MassRange range : ranges) {
            Interval ith = Interval.massRangeToInterval(range);
            intervals.add(ith);
        }
        List<Interval> mergedIntervals = MergeIntervals.mergeIntervals(intervals);

        
        for (Interval massInterval : mergedIntervals)
        {

            ResultsCollector collector = new ResultsCollector();

            float minMass = massInterval.getStart()*1000;
            float maxMass = massInterval.getEnd()*1000;
            
             
      
            
            BasicDBObject query = new BasicDBObject("MASS", new BasicDBObject("$gte", (int)(minMass)).
                    append("$lte", (int)(maxMass)));
            DBCursor cursor = dbcollection.find(query);
          //  System.out.println("The Start is " + minMass + " end is " + maxMass);

            try {
               // collector.collect(cursor);
                
              //  System.out.println("\nThe size of cursor "+cursor.size());
                
                
                
                List<IndexedSequence> massRangeResults = collector.getSequences(cursor);
                //note: massRangeResults might contain duplicate sequences
                //we should either filter out, or faster: add step during indexing to remove duplicates
                //and merge sequences, and preserve protein id information.  Check sqlite index implementation for example.
          //       allSeqs = collector.getSequences();


//                System.out.println("Printing after removing all the repeating occurances");
//                 for (IndexedSequence s : massRangeResults) {
//                    System.out.println(s.toString());
//                }
//                 System.out.println("The reduced size of cursor is "+ massRangeResults.size());
               allSeqs.addAll(massRangeResults);

          //      System.out.println("The size of the allSEQS is : "+ allSeqs.size());
                
                
                
                
            } catch (Exception ex) {
                logger.log(Level.SEVERE, "Cannot collect search results ", ex);

            } finally {
                cursor.close();
            }

        } //end for each range

        return allSeqs;

    }

    @Override
    public long addProteinDef(long num, String accession, String protSequence) throws DBIndexStoreException {
        //using protein cache entirely
        //nothing here
        return num;
    }

    @Override
    public void setProteinCache(ProteinCache proteinCache) {
        this.protCache = proteinCache;
    }

    @Override
    public boolean supportsProteinCache() {
        return true;
    }

    @Override
    public List<IndexedProtein> getProteins(IndexedSequence sequence) throws DBIndexStoreException {
        //get prot id from index and prot ids from cache
        List<IndexedProtein> ret = new ArrayList<IndexedProtein>();

        for (Integer protId : sequence.getProteinIds()) 
        {
           IndexedProtein ip = new IndexedProtein(protCache.getProteinDef(protId),(protId));
            ret.add(ip);
        }
//        for(int i= 0;i<sequence.getProteinIds().size();i++)
//        {
//            IndexedProtein ip = new IndexedProtein(protCache.getProteinDef(i), Long.valueOf(i));
//            ret.add(ip);
//        }
        

        return ret;


    }

    @Override
    public long getNumberSequences() throws DBIndexStoreException {
        if (!inited) {
            return 0;
        }

        long len = sequences.getCount();

        return len;
    }

    @Override
    public ResidueInfo getResidues(IndexedSequence peptideSequence, IndexedProtein protein) throws DBIndexStoreException {
        //in this implementation residues are already stored in the index
        return new ResidueInfo(peptideSequence.getResLeft(), peptideSequence.getResLeft());
    }

    /**
     * test driver
     *
     * @param args
     */
    public static void main(String[] args) {
        DBIndexStore store = new DBIndexStoreMongoDb();

        String protDef1 = "4R79.2 CE19650 WBGene00007067 Ras family status:Partially_confirmed TR:Q9XXA4 protein_id:CAA20282.1";
        String protDef2 = "Reverse_4R79.2  CE19650 WBGene00007067 Ras family status:Partially_confirmed TR:Q9XXA4 protein_id:CAA20282.1";

        try {
            logger.log(Level.INFO, "Testing init");
            File idx = new File("./test.fasta.idx2");
            if (idx.exists() && idx.isDirectory()) {
                for (File f : idx.listFiles()) {
                    f.delete();
                }
                idx.delete();
            }

            ProteinCache protCache = ProteinCache.getInstance();
            store.init("./test.fasta");
            store.setProteinCache(protCache);

            //test a transaction
            logger.log(Level.INFO, "Testing adds 1");
            store.startAddSeq();

            long protId = store.addProteinDef(0, protDef1, "ABCDEFGHIJKL");
            protCache.addProtein("4R79.2", "ABCDEFGHIJKL");

            store.addSequence(1f, 0, 1, "A", null, null, protId);
            store.addSequence(2f, 0, 2, "AB", null, null, protId);
            store.addSequence(3f, 0, 3, "ABC", null, null, protId);
            store.addSequence(4f, 0, 4, "ABCD", null, null, protId);
            store.addSequence(6000.42323f, 0, 5, "ABCDE", null, null, protId);
            store.addSequence(6999.42323f, 0, 6, "ABCDEZ", null, null, protId);
            store.addSequence(3, 6, 3, "GHI", null, null, protId);
            //store.addSequence(3, 6, 3, "GHI", null, null, protId);
            //store.stopAddSeq();
        } catch (DBIndexStoreException ex) {
            logger.log(Level.SEVERE, null, ex);
        }

        try {
            //test more adds with the same statement in a new transaction
            logger.log(Level.INFO, "Testing adds 2");
            //store.startAddSeq();

            long protId = store.addProteinDef(1, protDef2, "GHIJKLMNOPR");
            ProteinCache protCache = ProteinCache.getInstance();
            protCache.addProtein("Reverse_4R79.2", "GHIJKLMN");

            store.addSequence(3, 1, 3, "HIJ", null, null, protId);
            store.addSequence(5, 2, 5, "IJKLM", null, null, protId);
            //store.addSequence(7, 2, 7, "IJKLMNO", null, null, protId);
            //store.addSequence(7.1f, 2, 8, "IJKLMNOO", null, null, protId);
            store.addSequence(3, 0, 3, "GHI", null, null, protId); //test dup
            store.addSequence(3, 0, 3, "GHI", null, null, protId); //test dup
            store.stopAddSeq();
        } catch (DBIndexStoreException ex) {
            logger.log(Level.SEVERE, null, ex);
        }

        try {
            //test gets
            logger.log(Level.INFO, "Testing gets 1");
            List<IndexedSequence> res1 = store.getSequences(10, 8.9f);
            for (IndexedSequence seq : res1) {
                System.out.println(seq);
                List<IndexedProtein> proteins = store.getProteins(seq);
                for (IndexedProtein protein : proteins) {
                    ResidueInfo res = store.getResidues(seq, protein);
                    System.out.println(protein);
                    System.out.println(res);
                }
            }
        } catch (DBIndexStoreException ex) {
            logger.log(Level.SEVERE, null, ex);
        }


        try {
            //test gets
            System.out.println("Testing range gets 1");
            MassRange r1 = new MassRange(2, 1f);
            MassRange r2 = new MassRange(6, 1f);
            MassRange r3 = new MassRange(6, 1f); //test redundant
            MassRange r4 = new MassRange(6, 1.2f); //test overlapping
            List<MassRange> ranges = new ArrayList<MassRange>();
            ranges.add(r2);
            ranges.add(r1);
            ranges.add(r3);
            ranges.add(r4);
            List<IndexedSequence> res1 = store.getSequences(ranges);
            for (IndexedSequence seq : res1) {
                System.out.println(seq);
                List<IndexedProtein> proteins = store.getProteins(seq);
                for (IndexedProtein protein : proteins) {
                    ResidueInfo res = store.getResidues(seq, protein);
                    System.out.println(protein);
                    System.out.println(res);
                }
            }
        } catch (DBIndexStoreException ex) {
            logger.log(Level.SEVERE, null, ex);
        }


        if (true) {
            return;
        }


        try {
            store = new DBIndexStoreMongoDb();
            store.init("EBI-IPI_Human_IPI_Human_3_85_06-29-2011_reversed.fasta");
            for (int i = 0; i < 1000; ++i) {
                logger.log(Level.INFO, "Testing getting sequence " + i);
                List<IndexedSequence> res = store.getSequences(i, 3f);
                if (res.size() > 0) {
                    IndexedSequence seq = res.get(0);
                    logger.log(Level.INFO, "Got " + res.size() + " sequences, first sequence: " + seq);
                    List<IndexedProtein> proteins = store.getProteins(seq);
                    logger.log(Level.INFO, "Got " + res.size() + " sequences, first sequence's proteins: ");
                    for (IndexedProtein protein : proteins) {
                        System.out.println(protein);
                    }
                }

            }

        } catch (DBIndexStoreException ex) {
            logger.log(Level.SEVERE, null, ex);
        }

        try {
            logger.log(Level.INFO, "Testing getting multiple sequence proteins ");
            store = new DBIndexStoreMongoDb();
            store.init("test_02.fasta");

            List<IndexedSequence> res = store.getSequences(1000, 3000);
            for (IndexedSequence seq : res) {
                System.out.println(seq);
                List<IndexedProtein> prots = store.getProteins(seq);
                for (IndexedProtein prot : prots) {
                    System.out.println("\t" + prot);
                }
            }

        } catch (DBIndexStoreException ex) {
            logger.log(Level.SEVERE, null, ex);
        }

    }

    public class ResultsCollector {
        //to collapse multiple sequences into single one, with multproteins

        private final Map<String, List<IndexedSeqInternal>> temp = new HashMap<String, List<IndexedSeqInternal>>();

        public void collect(DBCursor cursor) throws IOException {

            DBObject currentCursor;
            List<IndexedSeqInternal> sequences = new ArrayList();
            List<BasicDBObject> sequenceList = new ArrayList<>();
            while (cursor.hasNext()) 
            {
                //System.out.println("\t" + currentCursor);
                currentCursor = cursor.next();
//                System.out.println("Printing total occurences of the specified mass");
//                System.out.println(currentCursor.get("MASS"));
                float massF = Float.parseFloat(currentCursor.get("MASS").toString());
                // sequences = (List<IndexedSeqInternal>) currentCursor.get("SEQUENCES");

                sequenceList = (List<BasicDBObject>) currentCursor.get("SEQUENCES");
                for (BasicDBObject currentsequence : sequenceList) {
                    String peptideSequence = (String) currentsequence.get("SEQ");
                    int offsetI = 0;
                    long protId = 1;
                    int lenI = (int) (long)(currentsequence.get("LENGTH"));
                    String protDesc = (String) currentsequence.get("RESLEFT");
                    //System.out.print("prot id: " + protIdI + " offset: " + offsetI + " len: " + lenI);
                    final IndexedSeqInternal tempSequence = new IndexedSeqInternal(massF, (int) offsetI, (int) lenI, (int) protId, peptideSequence, protDesc);
                    sequences = temp.get(peptideSequence);
                    if (sequences == null) 
                    {
                        sequences = new LinkedList<IndexedSeqInternal>();
                        temp.put(peptideSequence, sequences);
                    }
                    sequences.add(tempSequence);
                }
            }
        }

        List<IndexedSequence> getSequences(DBCursor cursor) 
        {
          
            DBObject currentCursor;
            List<IndexedSequence> sequences = new ArrayList();
            IndexedSequence tempSequence ;
            List<BasicDBObject> sequenceList = new ArrayList<>();
            while (cursor.hasNext()) 
            {
                currentCursor = cursor.next();
                float massF = Float.parseFloat(currentCursor.get("MASS").toString())/1000;
                sequenceList = (List<BasicDBObject>) currentCursor.get("SEQUENCES");
                Iterator itr = sequenceList.iterator();
                while (itr.hasNext())
                {
                    BasicDBObject currentsequence = (BasicDBObject) itr.next();
                    String peptideSequence = (String) currentsequence.get("SEQ");
                    List protID = new ArrayList();
                    List temp= new ArrayList();
                    temp=(List) currentsequence.get("PROT_ID");
                    Iterator itr1=temp.iterator();
                    while (itr1.hasNext())
                    {
                        protID.add(itr1.next());
                    }
                    //protID.addAll((Collection<Long>) currentsequence.get("PROT_ID"));
                    int length = (int) (currentsequence.get("LENGTH"));
                    String resLeft = (String) currentsequence.get("RESLEFT");
                    String resRight = (String) currentsequence.get("RESRIGHT");
                    //System.out.print("prot id: " + protIdI + " offset: " + offsetI + " len: " + lenI);
                    tempSequence = new IndexedSequence(massF, peptideSequence, resLeft, resRight, length, protID);
                    sequences.add(tempSequence);
                }

              
            }
              return sequences;
    }

    /*
    private class ResultsCollector {
        //to collapse multiple sequences into single one, with multproteins

        private final Map<String, List<IndexedSeqInternal>> temp = new HashMap<String, List<IndexedSeqInternal>>();

        public void collect(DBCursor cursor) throws IOException {


            for (DBObject currentCursor : cursor) {
                //System.out.println("\t" + currentCursor);

                double massF = (Double) currentCursor.get(FIELDS.MASS.toString());
                int offsetI = (Integer) currentCursor.get(FIELDS.OFFSET.toString());
                long protIdI = (Long) currentCursor.get(FIELDS.PROT_ID.toString());
                int lenI = (Integer) currentCursor.get(FIELDS.LEN.toString());

                String peptideSequence = null;
                //System.out.print("prot id: " + protIdI + " offset: " + offsetI + " len: " + lenI);

                peptideSequence = protCache.getPeptideSequence((int)protIdI, (int)offsetI, (int)lenI);
                //System.out.println("pep: " + peptideSequence);

                //we are cheating and supplying protein id instead of peptide id
                //to set it temporarily, before we merge them into a list
                final IndexedSeqInternal tempSequence = new IndexedSeqInternal((float)massF, (int) offsetI, (int)lenI, (int) protIdI, peptideSequence);
                List<IndexedSeqInternal> sequences = temp.get(peptideSequence);
                if (sequences == null) {
                    sequences = new LinkedList<IndexedSeqInternal>();
                    temp.put(peptideSequence, sequences);
                }
                sequences.add(tempSequence);

            }

        }

        List<IndexedSequence> getSequences() {
            List<IndexedSequence> ret = new LinkedList<IndexedSequence>();
            //group the same peptides from many proteins into single peptide with protein id list
            for (String pepSeqKey : temp.keySet()) {

                List<IndexedSeqInternal> sequences = temp.get(pepSeqKey);

                Set<Integer> proteinIds = new HashSet<Integer>();
                for (IndexedSeqInternal tempSeq : sequences) {
                    proteinIds.add((int) tempSeq.proteinId);
                }
                IndexedSeqInternal firstSeq = sequences.get(0);
                IndexedSequence mergedSequence = new IndexedSequence(0, firstSeq.mass, pepSeqKey, "", "");
                mergedSequence.setProteinIds(new ArrayList<Integer>(proteinIds));
                //set residues
                final String protSequence = protCache.getProteinSequence(firstSeq.proteinId);
                ResidueInfo residues = Util.getResidues(null, firstSeq.offset, firstSeq.length, protSequence);
                mergedSequence.setResidues(residues);
                ret.add(mergedSequence);
            }

            return ret;
        }
    }
    **/
    }
 public class addData
    {
     // Add the data to the buffer until the max buffer limit reached....
        void datatoBuffer(int mass,String sequence,String resLeft,String resRight,int sequenceLen,long proteinId)
        {
            int index = 0;
            boolean foundSequence=false;
            int seqMatchPosition = 0;
            boolean foundProtein = false;
//             Map<Integer ,ArrayList<BufferSequences>> temphash = new HashMap<>();
                BufferSequences currentSeq=new BufferSequences();
                ArrayList<BufferSequences> temp = new ArrayList<>();
//                if(sequence.equals("NPRTSPEKKATQTI"))
//                         System.out.println("MASS"+ mass);
//                if(mass==1494780)
//                   System.out.println("printing the buffer : SEQ:"  );
            // We found current entry of MASS in the BUFFER
            if(hashBuffer.containsKey(mass))
            {
                
                temp= hashBuffer.get(mass);
                 while(index <temp.size())
                 {
                     currentSeq =(BufferSequences)temp.get(index);
                     //We find the MASS and the specific sequence entry in BUFFER
                      if(currentSeq.getSeq().equals(sequence))
                      {
                          foundSequence = true;
                          seqMatchPosition=index;
                          if(currentSeq.getProtId().contains(proteinId))
                              foundProtein=true;
                      }
                      index++;
                 }
                 if(foundSequence == true && foundProtein==false)
                 {
//                     if(mass == 1320716)
//                         System.out.println();
                     currentSeq= temp.get(seqMatchPosition);
                     currentSeq.addprotId((int) proteinId);
                     temp.set(seqMatchPosition, currentSeq);
                     hashBuffer.put(mass, temp);
                 }
                 else if(foundSequence== false)
                 {
                     currentSeq= new BufferSequences();
                     currentSeq.setData(sequence, sequenceLen, resLeft, resRight, (int) proteinId);
                     temp.add(currentSeq);
                     //  massList.add(mass);
                     hashBuffer.put(mass, temp);
                 }
            }    
            //WE could not find currentMASS entry in the BUFFER
            else
            {
                currentSeq.setData(sequence, sequenceLen, resLeft, resRight, (int) proteinId);
                temp.add(currentSeq);
                massList.add(mass);
                hashBuffer.put(mass, temp);
            }
            
        }
// Once the buffer limit reached empty the buffer to the database and clear the buffer
// Remember one thing at the last run case the buffer might be having data lless then max bufffer limit
//       so you will have to call the method lastBuffertoDatabse thats is in the dbindexstoremongodb..        
        public void BuffertoDatabase() 
        {
            Iterator itr = massList.iterator();
            String seq;
            int mass;
            String current_id;
            List<BasicDBObject> seqList = new ArrayList();
            BasicDBObject dbnext = new BasicDBObject();
            List<BufferSequences>  currentBufferSeq = new ArrayList<>();
            boolean foundSequence=false,foundProtein=false;
            //loop for the total list of MASS that are in the buffer
            while(itr.hasNext())
            {
                mass=(int) itr.next();
                seq=null;
                foundSequence=false;foundProtein=false;
                BasicDBObject whereQuery = new BasicDBObject().append("MASS",mass);
//                whereQuery.put("MASS,",itr.next());
                DBCursor cursor = (DBCursor) sequences.find(whereQuery);
                BasicDBObject currentCursor =new BasicDBObject();
                //if the MASS is found in the database it will make below condition ture....
                if(cursor.hasNext()) 
                {
                    currentCursor= (BasicDBObject) cursor.next();
                    seqList=(List) currentCursor.get("SEQUENCES");
                  //   dbnext = (BasicDBObject) currentCursor.get("SEQUENCES");
                    // System.out.println(dbnext);
                    
                    currentBufferSeq= hashBuffer.get(mass);
                    
                    for(int k=0;k<currentBufferSeq.size();k++)
                    {
                        BufferSequences currentSeq =currentBufferSeq.get(k);
//                        if(currentSeq.getSeq().equals("NPRTSPEKKATQTI"))
//                         System.out.println("MASS"+ mass);
                        //List tobeAdded = new ArrayList();
                        Map<String,HashSet> tobeAdded = new HashMap();
                        HashSet<String> tobeAddedSeq = new HashSet<>();
                        foundProtein=false;foundSequence=false;
                        for(int j=0;j<seqList.size();j++)
                        {
                            seq= seqList.get(j).getString("SEQ");
                            //did we find the sequence in the databse????
//                            if(currentSeq.getSeq().equals("GRGLDLSS"))
//                                    System.out.println("MASS"+ mass);
                            if (currentSeq.getSeq().equals(seq))
                            {
                                
                                 foundSequence = true;
                            //    seqMatchPosition=index;
                                 List bufferProt=currentSeq.getProtId();
                                 List databseProList = (List) seqList.get(j).get("PROT_ID");
                                 for(int z=0;z<bufferProt.size();z++)
                                     {
                                         foundProtein=false;
                                         for(int t=0;t<databseProList.size();t++)
                                         {
                                             
                                         if((bufferProt.get(z))==(databseProList.get(t)))
                                         {
//                                             System.out.println("Entry for the specific protine is found....");
                                             foundProtein=true;
                                             break;
                                         }
//                                         if(foundProtein==false)
//                                            tobeAdded.add(bufferProt.get(z));                                  }

                                     }
                                 if(foundProtein==false)
                                 {
                                           // tobeAdded.add(bufferProt.get(z));
                                           HashSet ar= new HashSet();
                                           if(!tobeAdded.containsKey(seq))
                                           {  
                                               ar.add(bufferProt.get(z));
                                               tobeAdded.put(seq, ar);
                                               tobeAddedSeq.add(seq);
                                           }
                                           else
                                           {
                                               ar=  tobeAdded.get(seq);
                                               ar.add(bufferProt.get(z));
                                              tobeAdded.put(seq, ar);
                                              tobeAddedSeq.add(seq);
                                           }
                                 }
                                     }
                                 break;
                            }
                        }

                        //This means we found sequence but not the protin id so add the protein id.
                        if(foundSequence == true && foundProtein==false)
                        {
                           
                                if(tobeAdded.containsKey(currentSeq.getSeq()))
                                {
                                current_id = currentCursor.getString("MASS");
                                //Writing the ANDquery...........
                                List<BasicDBObject> obj = new ArrayList<BasicDBObject>();
                                obj.add(new BasicDBObject().append("MASS", Integer.parseInt(current_id)));
                                obj.add(new BasicDBObject().append("SEQUENCES.SEQ", currentSeq.getSeq()));
                                BasicDBObject andquery = new BasicDBObject();
                                andquery.put("$and", obj);
                                // writing the query to push the PROT_ID document.....
                                BasicDBObject newDocument = new BasicDBObject("SEQUENCES.$.PROT_ID", tobeAdded.get(currentSeq.getSeq()).toArray());
                                BasicDBObject updateDoc = new BasicDBObject();
                                updateDoc.put("$pushAll", newDocument);
                              //  System.out.println("New Prot_id added:"+mass);
                                sequences.update(andquery, updateDoc);
                               // System.out.println("New Protein_ID is added:" + current_id);
                                }
                          // continue;
                        }
                        //Means we find the MASS but dont find the sequence 
                        else if (foundSequence== false)
                        {
                           
                            current_id = currentCursor.getString("_id");
                               BasicDBObject bd1=new BasicDBObject().append("SEQUENCES",new BasicDBObject()
                                                        .append("SEQ", currentSeq.getSeq())
                                                         .append("RESLEFT", currentSeq.getResLeft())
                                                          .append("RESRIGHT",currentSeq.getResRight())
                                                           .append("LENGTH",currentSeq.getLength())
                                                            .append("PROT_ID",currentSeq.getProtId()));
                      //         BasicDBObject addProtid = new BasicDBObject().append("$push", new BasicDBObject("SEQUENCES.PROT_ID", currentseq.get(j).getProtId()));
                                BasicDBObject addProtid = new BasicDBObject("$push", bd1);
                               BasicDBObject WhereQuery = new BasicDBObject().append("MASS", mass);
                               sequences.update(WhereQuery, addProtid);
                             //  System.out.println("New Sequence is added:" + current_id);
                             //  continue;
                        }
                    }
                }
                //We didnt find the MASS entry in the databse
                else
                {
                    ArrayList<BufferSequences> currentseqArrayList = new ArrayList();
                    currentseqArrayList= hashBuffer.get(mass);
                     
                    Iterator it1= currentseqArrayList.iterator();
                    List listseq =new ArrayList<>();
                    // LOop for the all the sequence for a specific MASS
                    while(it1.hasNext())
                    {
                        BufferSequences currentbuffer =(BufferSequences) it1.next();
//                        if(currentbuffer.getSeq().equals("NPRTSPEKKATQTI"))
//                         System.out.println("MASS"+ mass);
                        BasicDBObject currentseq= new BasicDBObject().append("SEQ", currentbuffer.getSeq())
                                                         .append("RESLEFT", currentbuffer.getResLeft())
                                                          .append("RESRIGHT", currentbuffer.getResRight())
                                                           .append("LENGTH", currentbuffer.getLength())
                                                            .append("PROT_ID", currentbuffer.getProtId());
                        listseq.add(currentseq);
//                        System.out.println("SEQ:" + currentbuffer.getSeq() + mass);
                        
                    }
                     BasicDBObject addWholeMass = new BasicDBObject().append("MASS", mass)
                                                                    .append("SEQUENCES", listseq);   
                    sequences.insert(addWholeMass);
                    totalDatabaseEntry++;
//                    System.out.println(mass +"MASS Added to the databse:"+ totalDatabaseEntry );
                }
                                       
            
                cursor.close();
                }
             
            }
        }
    }


