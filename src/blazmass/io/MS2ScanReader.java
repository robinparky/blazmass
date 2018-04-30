package blazmass.io;

import blazmass.dbindex.Util;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * MS2 scan reader supports iterator interface and buffering on next
 * SCAN_CACHE_SIZE scans
 *
 * Thread-safe
 *
 * @author Adam
 */
public class MS2ScanReader implements Iterator<MS2Scan> {

    public static final int SCAN_CACHE_SIZE = 256;
    private int numScans;
    private String filePath;
    private final MS2Scan[] cache = new MS2Scan[SCAN_CACHE_SIZE];
    private volatile int curCacheI;
    private volatile int curCacheSize;
    private BufferedReader reader;
    private boolean EOF = false;
    private volatile boolean checkedNext;
    private static final Logger logger = Logger.getLogger(MS2ScanReader.class.getName());
    private static final Object lock = new Object(); //ensure single reader on a resource at a time
    private String curLine;  //current line for peeking
    private static final String SCAN_INDEX_EXT = ".index";

    public MS2ScanReader(String filePath) throws IOException {
        this.filePath = filePath;
        this.numScans = 0;
        this.checkedNext = false;
        this.curCacheSize = 0;
        this.curCacheI = -1;
        this.EOF = false;
        reader = new BufferedReader(new FileReader(filePath));

        synchronized (lock) {
            updateCache();
        }

    }
    
    /**
     * get total number of scans from index file if exists, or -1 if it does not (unknown)
     * @return total number of scans from the index
     */
    public int getNumScansIdx() {
        String indexFilePath = filePath + SCAN_INDEX_EXT;
        File indexFile = new File(indexFilePath);
        if (! indexFile.exists() || ! indexFile.canRead()) {
//            logger.log(Level.WARNING, "Cannot get total number of scans, index file cannot be read: " + indexFilePath);
            return -1;
        }
        
        
        
        int totalNumScans;
        try {
            totalNumScans = Util.countLines(indexFilePath);
        } catch (IOException ex) {
            logger.log(Level.WARNING, "Cannot get total number of scans, index file lines cannot be counted: " + indexFilePath);
            return -1;
        }
        
        return totalNumScans;
    }

    public String getFilePath() {
        return filePath;
    }

    public int getNumScans() {
        return numScans;
    }

    private void updateCache() {
        if (curCacheSize != 0 || EOF) {
            return;
        }

        //curCacheSize is 0
        curCacheI = 0;

        //read more spectra into cache
        MS2Scan curScan = null;
        try {
            if (curLine == null) {
                //first time
                curLine = reader.readLine();
            }

            while (curLine != null) {
                curLine = curLine.trim();
                if (curLine.length() == 0) {
                    continue;
                }

                final char firstChar = curLine.charAt(0);
                if (firstChar == 'H' || firstChar == 'I' || firstChar == 'D') {                    
                    curLine = reader.readLine();
                    if(curLine.startsWith("I\tActivationType\t")) {
                        String[] arr = curLine.split("\t");
                        if(arr.length>=3)
                            curScan.setScanType(arr[2]);
                    }
                    
                    continue;
                } else if (firstChar == 'S') {
                    //new scan, create and add to cache
                    //will fill in rest of data next
                    curScan = new MS2Scan();
                    cache[curCacheSize] = curScan;
                    ++curCacheSize;  //scans read in this updateCache()
                    ++numScans; //total number scans read from this file


                    String[] tokens = curLine.split("\t");

                    if (tokens.length < 3) {
                        throw new IOException("Unexpected number of tokens on S line (need 3): " + curLine);
                    }

                    try {
                        curScan.setIsScan1(Integer.parseInt(tokens[1]));
                        curScan.setIsScan2(Integer.parseInt(tokens[2]));
                    } catch (NumberFormatException e) {
                        throw new IOException("Can't get integers from S line: " + curLine, e);
                    }

                    curLine = reader.readLine();

                } else if (firstChar == 'Z') {
                    if (curScan == null) {
                        throw new IOException("Unexpected Z line before S line");
                    }

                    String[] tokens = curLine.split("\t");
                    if (tokens.length < 3) {
                        throw new IOException("Unexpected number of tokens on Z line (need at least 3): " + curLine);
                    }

                    try {
                        curScan.addChargeState(Integer.parseInt(tokens[1]));
                        curScan.addPrecursorMass(Float.parseFloat(tokens[2]));
                    } catch (NumberFormatException e) {
                        throw new IOException("Can't get values from Z line: " + curLine, e);
                    }

                    curLine = reader.readLine();
                } else {
                    //read mass / intensities
                    if (curScan == null) {
                        throw new IOException("Unexpected line, expected masses: " + curLine + " in scan number: " + numScans);
                    }

                    //read until next S
                    while (curLine != null) {
                        if (curLine.startsWith("S")) {
                            curScan = null;
                            break;
                        }
                        
                        String[] tokens = curLine.split(" ");
                        //if (tokens.length != 2) {  //it can have three columns
                        if (tokens.length < 2) {
                            throw new IOException("Unexpected number of tokens on mass line (need 2): " + curLine + " in scan number: " + numScans);
                        }

                        try {
                            curScan.addMassIntensity(Float.parseFloat(tokens[0]), Float.parseFloat(tokens[1]));
                        } catch (NumberFormatException e) {
                            throw new IOException("Can't get mass / intensities from line: " + curLine + " in scan number: " + numScans, e);
                        }

                        //either next mass or new S line
                        curLine = reader.readLine();
                    }
                    
                    
                    //check if new spectrum
                    if (curScan == null) {
                        //new S line, check if cache is full
                        if (curCacheSize == SCAN_CACHE_SIZE) {
                            break;
                            //curLine with S scan will be used next time
                        }
                    }

                } //end reading masses
            } //and while for each line

            //logger.log(Level.INFO, "Read into cache scans: " + curCacheSize);
            //check if EOF
            if (curCacheSize < SCAN_CACHE_SIZE) {
                EOF = true;
                close();
            }

            
        } catch (IOException e) {
            logger.log(Level.WARNING, "Error reading MS2 file: " + filePath, e);
            EOF = true;
            close();
        }
    }

    @Override
    public boolean hasNext() {
        checkedNext = true;
        return curCacheSize != 0;

    }

    @Override
    public MS2Scan next() {

        if (checkedNext == false) {
            //throw new IllegalStateException("Need to check hasNext() first!");
        }
        checkedNext = false;

        synchronized (lock) {
            //updateCache();
            
            if (curCacheSize == 0) {
                return null;
            }

            final MS2Scan scan = cache[curCacheI++];
            --curCacheSize;
            
            updateCache();

            return scan;
        }
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    private void close() {
        try {
            if (reader != null) {
                reader.close();
                reader = null;
            }
        } catch (IOException ex) {
            logger.log(Level.SEVERE, "Error closing reader", ex);
        }
        EOF = true;
    }

    //test driver
    public static void main(String[] args) {

        if (args.length == 0) {
            System.out.println("Need ms2 file name");
        }
        try {
            MS2ScanReader reader = new MS2ScanReader(args[0]);
            int count = 0;
            while (reader.hasNext()) {
                System.out.println("Scan#: " + (++count));
                MS2Scan scan = reader.next();
                System.out.println(scan);
            }
            System.out.println("total scans: " + reader.getNumScans());
        } catch (IOException ex) {
            Logger.getLogger(MS2ScanReader.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
