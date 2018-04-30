package blazmass.io;

import java.io.IOException;

/**
 * Thread safe result writer to stdout
 * 
 * @author Adam
 */
public class StdoutResultWriter implements ResultWriter {

    @Override
    public synchronized void write(String toWrite){
        System.out.print(toWrite);
        System.out.flush();
    }

    @Override
    public void flush() {
        
    }

    
    @Override
    public void close() {
    }
}
