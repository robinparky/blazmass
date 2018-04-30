package blazmass.dbindex;

import blazmass.io.SearchParamReader;
import edu.scripps.pms.mspid.SearchParams;
import org.jdom.JDOMException;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by rpark on 3/18/16.
 */
public class IndexQuery {

    private static final Logger logger = Logger.getLogger(IndexQuery.class.getName());
    public static void main(String[] args) throws IOException, JDOMException, DBIndexerException {
        System.out.println("aaaa");

      //  IndexQuery.queryProlucidParam();
      //  if(true) return;

        //String path = "/home/rpark/test_data/search/small";

        String path = "/data/2/rpark/ip2_data/rpark/Jolene_Hela_loading/20141101_HeLa_1ug_BEH60_140_35_IC_DE5_5e3_1_2014_11_12_16_2015_03_19_17_30916/search/projects2016_04_01_16_95760/test";
        String filename = "search.xml";
        //SearchParams sp = new SearchParams(param);
        SearchParamReader paramReader = new SearchParamReader(path, filename);
        blazmass.io.SearchParams param = paramReader.getParam();
        //param.setIndexType(DBIndexer.IndexerMode.SEARCH_INDEXED);


            //final DBIndexer.IndexerMode indexerMode = sParam.isUseIndex()? DBIndexer.IndexerMode.SEARCH_INDEXED: DBIndexer.IndexerMode.SEARCH_UNINDEXED;
        final DBIndexer indexer = new DBIndexer(param);

        List<IndexedSequence> sequenceList = indexer.getSequences(1000f, 2000f);
for(IndexedSequence seq: sequenceList)
        System.out.println(seq);
        /*
        System.out.println(sequenceList.get(0));
        sequenceList = indexer.getSequences(1234.2f, 0.8f);
        System.out.println(sequenceList.get(0));

        sequenceList = indexer.getSequences(1234.2f, 0.8f);
        System.out.println(sequenceList.get(0));

        sequenceList = indexer.getSequences(1234.2f, 0.8f);
        System.out.println(sequenceList.get(0));

        sequenceList = indexer.getSequences(1234.2f, 0.8f);
        System.out.println(sequenceList.get(0));*/


        /*
        for(Iterator<IndexedSequence> itr=sequenceList.iterator(); itr.hasNext(); ) {
            System.out.println(itr.next());
        }
        */
    }

    public static void queryProlucidParam() throws IOException, JDOMException, DBIndexerException {
        String path = "/home/rpark/test_data/search/small";
        String filename = "search.xml";
        SearchParamReader reader = new SearchParamReader(path, filename);

        System.out.println("test");

        final DBIndexer indexer = new DBIndexer( reader.getSearchParams() );

        List<IndexedSequence> sequenceList = indexer.getSequences(1234.2f, 0.8f);

        System.out.println(sequenceList.get(0));
        sequenceList = indexer.getSequences(1234.2f, 0.8f);
        System.out.println(sequenceList.get(0));

    }
}
