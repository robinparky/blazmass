package blazmass.io;


/**
 * Interface for search result writers
 * 
 * @author Adam
 */
public interface ResultWriter {

    public void write(String toWrite);

    public void close();
    
    public void flush();
}
