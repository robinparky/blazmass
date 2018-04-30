package blazmass;

public class Options {       

	 private  int  peptideLines;
	 private boolean  useEnzyme;
	 private  int  removePrecursor;
	 public Options(int peptideLines,boolean useEnzyme, int removePrecursor ){
		 this.peptideLines=peptideLines;
		 this.useEnzyme=useEnzyme;
		 this.removePrecursor=removePrecursor;
	 }
	 public Options( ){
		 peptideLines=-1;
		 useEnzyme=false;
		 removePrecursor=-1;
	 }

    public int getPeptideLines() {
        return peptideLines;
    }

    public void setPeptideLines(int peptideLines) {
        this.peptideLines = peptideLines;
    }

    public int getRemovePrecursor() {
        return removePrecursor;
    }

    public void setRemovePrecursor(int removePrecursor) {
        this.removePrecursor = removePrecursor;
    }

    public boolean isUseEnzyme() {
        return useEnzyme;
    }

    public void setUseEnzyme(boolean useEnzyme) {
        this.useEnzyme = useEnzyme;
    }

         
}

