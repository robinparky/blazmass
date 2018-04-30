package blazmass.io;

import java.io.UnsupportedEncodingException;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.DataFormatException;
import java.util.zip.Deflater;
import java.util.zip.Inflater;

/**
 * Utils for compressing and decompressing strings
 */
public class CompressZlib {

    private static final String ENCODING = "UTF-8";
    
    private static final Logger logger = Logger.getLogger(CompressZlib.class.getName());
    
    public byte[] compress(byte[] bytesToCompress) {
        Deflater deflater = new Deflater();
        deflater.setInput(bytesToCompress);
        deflater.finish();

        byte[] bytesCompressed = new byte[Short.MAX_VALUE];

        int numberOfBytesAfterCompression = deflater.deflate(bytesCompressed);

        byte[] returnValues = new byte[numberOfBytesAfterCompression];

        System.arraycopy(
                bytesCompressed,
                0,
                returnValues,
                0,
                numberOfBytesAfterCompression);

        return returnValues;
    }

    public byte[] compress(String stringToCompress) {
        byte[] returnValues = null;

        try {

            returnValues = this.compress(
                    stringToCompress.getBytes(ENCODING));
        } catch (UnsupportedEncodingException ex) {
            logger.log(Level.SEVERE, "Error decompressing string", ex);
        }

        return returnValues;
    }

    public byte[] decompress(byte[] bytesToDecompress) {
        Inflater inflater = new Inflater();

        int numberOfBytesToDecompress = bytesToDecompress.length;

        inflater.setInput(
                bytesToDecompress,
                0,
                numberOfBytesToDecompress);

        int compressionFactorMaxLikely = 3;

        int bufferSizeInBytes =
                numberOfBytesToDecompress
                * compressionFactorMaxLikely;

        byte[] bytesDecompressed = new byte[bufferSizeInBytes];

        byte[] returnValues = null;

        try {
            int numberOfBytesAfterDecompression = inflater.inflate(bytesDecompressed);

            returnValues = new byte[numberOfBytesAfterDecompression];

            System.arraycopy(
                    bytesDecompressed,
                    0,
                    returnValues,
                    0,
                    numberOfBytesAfterDecompression);
        } catch (DataFormatException ex) {
            logger.log(Level.SEVERE, "Error decompressing string", ex);
        }

        inflater.end();

        return returnValues;
    }

    public String decompressToString(byte[] bytesToDecompress) {
        byte[] bytesDecompressed = this.decompress(
                bytesToDecompress);

        String returnValue = null;

        try {
            returnValue = new String(
                    bytesDecompressed,
                    0,
                    bytesDecompressed.length,
                    ENCODING);
        } catch (UnsupportedEncodingException ex) {
            logger.log(Level.SEVERE, "Error decompressing string", ex);
        }

        return returnValue;
    }
}
