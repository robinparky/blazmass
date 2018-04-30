package blazmass.io;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Thread safe file result writer with a buffer
 *
 * @author Adam
 */
public class FileResultWriter implements ResultWriter {

    private volatile BufferedWriter fileWriter;
    private String fileName;
    private static final Logger logger = Logger.getLogger(FileResultWriter.class.getName());
    private static final int BUFF_SIZE = 1024 * 128;
//    private static final int BUFF_SIZE = 128;

    public FileResultWriter(String resultFile) throws IOException {
        fileName = resultFile;
        fileWriter = new BufferedWriter(new FileWriter(resultFile), BUFF_SIZE);
    }

    @Override
    public void write(String toWrite) {
        synchronized (this) {
            try {
                fileWriter.write(toWrite);
            } catch (IOException ex) {
                logger.log(Level.SEVERE, "Cannot write result to file: " + fileName, ex);
            }
        }
    }

    @Override
    public void close() {
        closeInternal();
    }

    @Override
    public void flush() {
        if (fileWriter != null) {
            synchronized (this) {
                try {
                    fileWriter.flush();
                } catch (IOException ex) {
                    logger.log(Level.SEVERE, "Error flushing result writer", ex);
                }
            }
        }
    }

    private void closeInternal() {
        if (fileWriter != null) {
            try {
                synchronized (this) {

                    fileWriter.flush();
                    fileWriter.close();
                    fileWriter = null;
                }

            } catch (IOException ex) {
                logger.log(Level.SEVERE, "Cannot close result file: " + fileName, ex);
            }
        }
    }
}
