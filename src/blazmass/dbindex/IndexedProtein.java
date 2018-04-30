package blazmass.dbindex;

import blazmass.io.Fasta;

/**
 * A representation of an indexed protein (a fasta header)
 * 
 * @author Adam
 */
public class IndexedProtein {
    private long id;
    private String accession;

    
    IndexedProtein(long id) {
        this.id = id;
    }
    
    IndexedProtein(String fastaDefLine, long id) {
        //this.accession = Fasta.getAccession(fastaDefLine);
        this.accession = Fasta.getSequestLikeAccession(fastaDefLine);
        this.id = id;
    }


    public long getId() {
        return id;
    }

    public void setId(long id) {
        this.id = id;
    }
    
   
    @Override
    public String toString() {
        return "IndexedProtein{" + "accession=" + accession + ", id=" + id + '}';
    }
   

    public String getAccession() {
        return accession;
    }

    public void setAccession(String accession) {
        this.accession = accession;
    }

    
    
}
