import java.awt.Desktop;
import java.awt.Dimension;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.math.BigInteger;
import java.util.Vector;

import javax.media.opengl.GLOffscreenAutoDrawable;
import javax.media.opengl.GLProfile;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JProgressBar;
import javax.swing.SwingWorker;

import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.log4j.BasicConfigurator;

import com.algosome.eutils.blast.BlastParser;
import com.algosome.eutils.blast.GetCommand;
import com.algosome.eutils.blast.PutCommand;
import com.icafe4j.image.ImageIO;
import com.icafe4j.image.ImageParam;
import com.icafe4j.image.ImageType;
import com.icafe4j.image.options.TIFFOptions;
import com.icafe4j.image.tiff.TiffFieldEnum.Compression;
import com.icafe4j.image.tiff.TiffFieldEnum.PhotoMetric;
import com.icafe4j.image.writer.ImageWriter;
import com.metsci.glimpse.canvas.FBOGlimpseCanvas;
import com.metsci.glimpse.examples.Example;
import com.metsci.glimpse.gl.util.GLUtils;
import com.metsci.glimpse.plot.ColorAxisPlot2D;
import com.metsci.glimpse.support.settings.SwingLookAndFeel;

public class main{
    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//Variable initialization
    public static int nmer = 3;
    public static int numShifts = 0;
    public static int numShiftsMinus = 0;
    public static int slidingSize = 5000;
    public static String file1 = "Pa_PA1R.fna";
    public static String file2 = "Pa_PA7.fna";
    public static float[][] regressionValues = null;
    public static Vector<Genes> sequencesX = null;
    public static Vector<Genes> sequencesY = null;
    public static boolean spearman = false;

    public static double numberIterationsX = 0;
    public static double numberIterationsY = 0;
    public static String seqXname = "sequence1";
    public static String seqYname = "sequence2";
    public static String y = "";
    public static String x = "";
    public static double [] maxValues_X;
    public static double [] maxValues_Y;
    public static double [] averageValues_X;
    public static double [] averageValues_Y;
    public static int windowLength = 0;
    public static Vector<float[]> xStorage;
    public static Vector<float[]> yStorage;
    public static boolean labelControl = true;

    public static Vector<Integer> badX;	//these are where we store windows that have more than 20% kmers that are bad
    public static Vector<Integer> badY;
    public static Vector<Integer> interestingX = new Vector<Integer>();
    public static Vector<Integer> interestingY = new Vector<Integer>();
    public static HeatMapExample a = new HeatMapExample();
    
    public static boolean byGene = false;
    public static boolean byExample = false;
    
    //Code starts here
    public static void main(String[] args) throws Exception {        
        //GUI runs, inputting from user. After input, goes to RunOurBaby() 
        FrontGUI window = new FrontGUI();
        window.runGUI();
    }
    
    //Runner method to get kmers from window: ignores kmers with n's
    public static float[] runGetKmers (String sequence) {
    	float[] kmerComp = null;
        try {
        	kmerComp = new float[(int) Math.pow(4, nmer)];
        }
        catch (OutOfMemoryError e) {
            JOptionPane.showMessageDialog(null, "Out of Memory. Please Rerun SPlot with a smaller kmer size.");
            System.exit(0);
        }
        String[] toRun = sequence.split("N");
        for (int i =0; i<toRun.length; i++) {
            if (toRun[i].length()>=nmer) {
                kmerComp = getKmers(toRun[i], kmerComp);
                }
            }
        return kmerComp;
    }

    /*
     * Automated BLAST method
     * Looks through the array of correlation values generated for each window in X and Y
     * If the averageValue of that window and maxValue of that window is below a certain threshold, sends a 
     * BLAST query to NCBI via runBLAST method
     */
    public static void automatedBLAST (double maxCorrelation, double averageCorrelation, boolean xSwitch, boolean ySwitch, 
    		String fileX, String fileY) throws IOException {
        int intIterationsX = (int)numberIterationsX;
        int intIterationsY = (int)numberIterationsY;
        //System.out.println();
        //System.out.println("averageValues_X");
        for(int xv = 0;xv<intIterationsX;xv++){
            averageValues_X[xv] = averageValues_X[xv]/intIterationsX;
            //if thresholds met, add the indice to interesting X for later blasting
            if (averageValues_X[xv]<averageCorrelation && maxValues_X[xv]<maxCorrelation) {
                interestingX.add(xv);
                System.out.print(averageValues_X[xv]);
            }
        }
        
        //System.out.println();
        //System.out.println("averageValues_Y");
        for(int xv = 0;xv<intIterationsY;xv++){
            averageValues_Y[xv] = averageValues_Y[xv]/intIterationsY;
            //if thresholds met, add the indice to interesting Y for later blasting
            if (averageValues_Y[xv]<averageCorrelation && maxValues_Y[xv]<maxCorrelation) {
                interestingY.add(xv);
                //System.out.print(averageValues_Y[xv]);
            }
        }
        //System.out.println();
        fileX = fileX.replace(".txt", "");
        fileY = fileY.replace(".txt", "");
        
        if (ySwitch==true) {
	        for (int i = 0; i<interestingY.size(); i++) {
	            int temp = interestingY.get(i);
	            boolean stillContig = true;
	            while (stillContig){
	                if(interestingY.contains(temp+1)){
	                    interestingY.removeElement(temp+1);
	                    temp++;
	                }
	                else{
	                    stillContig = false;
	                }
	            }
	            String seq = y.substring(interestingY.get(i) * slidingSize, (temp+1) * slidingSize);
	            seq = seq.toUpperCase();
	            String filename = fileY + "_" + interestingY.get(i)+".txt";
	            try {
	                runBlast(seq, filename);
	            } catch (FileNotFoundException e1) {
	                e1.printStackTrace();
	            }
	        }
        }
        
        if (xSwitch==true) {
	        for (int i = 0; i<interestingX.size(); i++) {
	            int temp = interestingX.get(i);
	            boolean stillContig = true;
	            while (stillContig){
	                if(interestingX.contains(temp+1)){
	                    interestingX.removeElement(temp+1);
	                    temp++;
	                }
	                else{
	                    stillContig = false;
	                }
	            }
	            String seq = x.substring(interestingX.get(i) * slidingSize, (temp + 1) * slidingSize);
	            seq = seq.toUpperCase();
	            String filename = fileX + "_" + interestingX.get(i)+".txt";
	            try {
	                runBlast(seq, filename);
	            } catch (FileNotFoundException e1) {
	                e1.printStackTrace();
	            }
	        }
        }
    }
    
    public static void automatedBLASTGenes (double maxCorrelation, double averageCorrelation, boolean xSwitch, boolean ySwitch, 
    		String fileX, String fileY) throws IOException {
        int intIterationsX = (int)numberIterationsX;
        int intIterationsY = (int)numberIterationsY;
        //System.out.println();
        //System.out.println("averageValues_X");
        for(int xv = 0;xv<intIterationsX;xv++){
            averageValues_X[xv] = averageValues_X[xv]/intIterationsX;
            //if thresholds met, add the indice to interesting X for later blasting
            if (averageValues_X[xv]<averageCorrelation && maxValues_X[xv]<maxCorrelation) {
                interestingX.add(xv);
                //System.out.print(averageValues_X[xv]);
            }
        }
        
       // System.out.println();
        //System.out.println("averageValues_Y");
        for(int xv = 0;xv<intIterationsY;xv++){
            averageValues_Y[xv] = averageValues_Y[xv]/intIterationsY;
            //if thresholds met, add the indice to interesting Y for later blasting
            if (averageValues_Y[xv]<averageCorrelation && maxValues_Y[xv]<maxCorrelation) {
                interestingY.add(xv);
               // System.out.print(averageValues_Y[xv]);
            }
        }
        //System.out.println();
        fileX = fileX.replace(".txt", "");
        fileY = fileY.replace(".txt", "");
        
        if (ySwitch==true) {
	        for (int i = 0; i<interestingY.size(); i++) {
	            String seq = sequencesY.get(interestingY.get(i)).sequence;
	            seq = seq.toUpperCase();
	            String filename = fileY + "_" + interestingY.get(i) + ".txt";
	            try {
	                runBlast(seq, filename);
	            } catch (FileNotFoundException e1) {
	                e1.printStackTrace();
	            }
	        }
        }
        
        if (xSwitch==true) {
	        for (int i = 0; i<interestingX.size(); i++) {
	            String seq = sequencesX.get(interestingX.get(i)).sequence;
	            seq = seq.toUpperCase();
	            String filename = fileX + "_" + interestingX.get(i)+".txt";
	            try {
	                runBlast(seq, filename);
	            } catch (FileNotFoundException e1) {
	                e1.printStackTrace();
	            }
	        }
        }
    }
    
    /*
     * Similar to autoBLAST,
     * but instead, this method writes all interesting sequences to file.
     */
    public static void automatedFile (double maxCor, double averageCor, boolean xSwitch, boolean ySwitch, 
    		String xFile, String yFile) throws FileNotFoundException {
        int intIterationsX = (int)numberIterationsX;
        int intIterationsY = (int)numberIterationsY;
        Vector<Double> maxCorStorageX = new Vector<Double>();
        Vector<Double> avgCorStorageX = new Vector<Double>();
        Vector<Double> maxCorStorageY = new Vector<Double>();
        Vector<Double> avgCorStorageY = new Vector<Double>();
        //System.out.println();
        //System.out.println("averageValues_X");
        for(int xv = 0;xv<intIterationsX;xv++){
            averageValues_X[xv] = averageValues_X[xv]/intIterationsX;
            if (averageValues_X[xv]<averageCor && maxValues_X[xv]<maxCor) {
                interestingX.add(xv);
                maxCorStorageX.add(maxValues_X[xv]);
                avgCorStorageX.add(averageValues_X[xv]);
                //System.out.print(averageValues_X[xv]);
            }
        }
        //System.out.println();
       //System.out.println("averageValues_Y");
        for(int xv = 0;xv<intIterationsY;xv++){
            averageValues_Y[xv] = averageValues_Y[xv]/intIterationsY;
            if (averageValues_Y[xv]<averageCor && maxValues_Y[xv]<maxCor) {
                interestingY.add(xv);
                maxCorStorageY.add(maxValues_Y[xv]);
                avgCorStorageY.add(averageValues_Y[xv]);
                //System.out.print(averageValues_Y[xv]);
            }
        }
        //System.out.println();
        String filename = "";
        
        if (ySwitch==true) {
	        filename = yFile;
	        PrintStream outt = new PrintStream(new File(filename));
	        for (int i = 0; i<interestingY.size(); i++) {
	            int temp = interestingY.get(i);
	            boolean stillContig = true;
	            while (stillContig){
	                if(interestingY.contains(temp+1)){
	                    interestingY.removeElement(temp+1);
	                    maxCorStorageY.removeElement(temp+1);
	                    avgCorStorageY.removeElement(temp+1);
	                    temp++;
	                }
	                else{
	                    stillContig = false;
	                }
	            }
	            String seq = y.substring(interestingY.get(i) * slidingSize, (temp+1) * slidingSize);
	            seq = seq.toUpperCase();
	            outt.println(">"+seqYname+":"+interestingY.get(i)*windowLength+"-"+(interestingY.get(i)*windowLength+slidingSize)+"|"+maxCorStorageY.get(i)+"|"+avgCorStorageY.get(i)+"|");
	            int numIter = seq.length()/70;
	            int j = 0;
	            for (j = 0; j<numIter; j++) {
	                outt.println(seq.substring(j*70, (j+1)*70));
	            }
	            if (j*70<seq.length()) {
	                outt.println(seq.substring((j)*70, seq.length()));
	            }
	        }
	        outt.close();
	        //then we want to open the file for the user, because we are nice guys
	        if (Desktop.isDesktopSupported()) {
	        	try {
					Desktop.getDesktop().open(new File(filename));
				} catch (IOException e) {
				
					e.printStackTrace();
				};
	        }
        }
        
        if (xSwitch==true) {
        filename = xFile;
        PrintStream out = new PrintStream(new File(filename));
        for (int i = 0; i<interestingX.size(); i++) {
            int temp = interestingX.get(i);
            boolean stillContig = true;
            while (stillContig && byGene==false){
                if(interestingX.contains(temp+1)){
                    interestingX.removeElement(temp+1);
                    maxCorStorageX.removeElement(temp+1);
                    avgCorStorageX.removeElement(temp+1);
                    temp++;
                }
                else{
                    stillContig = false;
                }
            }
            String seq = x.substring(interestingX.get(i) * slidingSize, (temp + 1) * slidingSize);
            seq = seq.toUpperCase();
            out.println(">"+seqXname+":"+interestingX.get(i)*windowLength+"-"+(interestingX.get(i)*windowLength+slidingSize)+"|"+maxCorStorageX.get(i)+"|"+avgCorStorageX.get(i)+"|");
            int numIter = seq.length()/70;
            int j = 0;
            for (j = 0; j<numIter; j++) {
                out.println(seq.substring(j*70, (j+1)*70));
            }
            if (j*50<seq.length()) {
                out.println(seq.substring((j)*70, seq.length()));
            }
        }
        out.close();
        
      //then we want to open the file for the user, because we are nice guys
        if (Desktop.isDesktopSupported()) {
        	try {
				Desktop.getDesktop().open(new File(filename));
			} catch (IOException e) {
				
				e.printStackTrace();
			};
        }
        }
        }
    
    public static void automatedFileByGene (double maxCor, double averageCor, boolean xSwitch, boolean ySwitch, 
    		String xFile, String yFile) throws FileNotFoundException {
        int intIterationsX = (int)numberIterationsX;
        int intIterationsY = (int)numberIterationsY;
        Vector<Double> maxCorStorageX = new Vector<Double>();
        Vector<Double> avgCorStorageX = new Vector<Double>();
        Vector<Double> maxCorStorageY = new Vector<Double>();
        Vector<Double> avgCorStorageY = new Vector<Double>();
       // System.out.println();
        //System.out.println("averageValues_X");
        for(int xv = 0;xv<intIterationsX;xv++){
            averageValues_X[xv] = averageValues_X[xv]/intIterationsX;
            if (averageValues_X[xv]<averageCor && maxValues_X[xv]<maxCor) {
                interestingX.add(xv);
                maxCorStorageX.add(maxValues_X[xv]);
                avgCorStorageX.add(averageValues_X[xv]);
               // System.out.print(averageValues_X[xv]);
            }
        }
       // System.out.println();
       // System.out.println("averageValues_Y");
        for(int xv = 0;xv<intIterationsY;xv++){
            averageValues_Y[xv] = averageValues_Y[xv]/intIterationsY;
            if (averageValues_Y[xv]<averageCor && maxValues_Y[xv]<maxCor) {
                interestingY.add(xv);
                maxCorStorageY.add(maxValues_Y[xv]);
                avgCorStorageY.add(averageValues_Y[xv]);
                //System.out.print(averageValues_Y[xv]);
            }
        }
        //System.out.println();
        String filename = "";
        
        if (ySwitch==true) {
	        filename = yFile;
	        PrintStream outt = new PrintStream(new File(filename));
	        for (int i = 0; i<interestingY.size(); i++) {
	            int temp = interestingY.get(i);
	            /*boolean stillContig = true;
	            while (stillContig){
	                if(interestingY.contains(temp+1)){
	                    interestingY.removeElement(temp+1);
	                    temp++;
	                }
	                else{
	                    stillContig = false;
	                }
	            }*/
	            String seq = sequencesY.get(interestingY.get(i)).sequence;
	            seq = seq.toUpperCase();
	            outt.println(">"+sequencesY.get(interestingY.get(i)).ID+":"+maxCorStorageY.get(i)+"|"+avgCorStorageY.get(i)+"|");
	            int numIter = seq.length()/70;
	            int j = 0;
	            for (j = 0; j<numIter; j++) {
	                outt.println(seq.substring(j*70, (j+1)*70));
	            }
	            if (j*70<seq.length()) {
	                outt.println(seq.substring((j)*70, seq.length()));
	            }
	        }
	        outt.close();
	        //then we want to open the file for the user, because we are nice guys
	        if (Desktop.isDesktopSupported()) {
	        	try {
					Desktop.getDesktop().open(new File(filename));
				} catch (IOException e) {
				
					e.printStackTrace();
				};
	        }
        }
        
        if (xSwitch==true) {
        filename = xFile;
        PrintStream out = new PrintStream(new File(filename));
        for (int i = 0; i<interestingX.size(); i++) {
            /*int temp = interestingX.get(i);
            boolean stillContig = true;
            while (stillContig){
                if(interestingX.contains(temp+1)){
                    interestingX.removeElement(temp+1);
                    temp++;
                }
                else{
                    stillContig = false;
                }
            }*/
            String seq = sequencesX.get(interestingX.get(i)).sequence;
            seq = seq.toUpperCase();
            out.println(">"+sequencesX.get(interestingX.get(i)).ID+":"+maxCorStorageX.get(i)+"|"+avgCorStorageX.get(i)+"|");
            int numIter = seq.length()/70;
            int j = 0;
            for (j = 0; j<numIter; j++) {
                out.println(seq.substring(j*70, (j+1)*70));
            }
            if (j*50<seq.length()) {
                out.println(seq.substring((j)*70, seq.length()));
            }
        }
        out.close();
        
      //then we want to open the file for the user, because we are nice guys
        if (Desktop.isDesktopSupported()) {
        	try {
				Desktop.getDesktop().open(new File(filename));
			} catch (IOException e) {
				
				e.printStackTrace();
			};
        }
        }
        maxCorStorageX.clear();
        avgCorStorageX.clear();
        maxCorStorageY.clear();
        avgCorStorageY.clear();
        interestingY.clear();
        interestingX.clear();
        
        
        
        }

    
    
    
    
    
    
    
    
    
    /*
     * Uses JeUtil package to send query to NCBI's Servers and gets 
     * files back with BLAST results
     */
    public static void runBlast(String seq, String filename) throws IOException {
        BasicConfigurator.configure();
        //logger.info("Blast utility test");
        PutCommand put = new PutCommand();
        PrintStream out = new PrintStream(new File(filename));
        put.setQuery(seq);
        put.setProgram("blastn");
        put.setDatabase("nr");
        
        GetCommand get = new GetCommand(new BlastParser(){
            @Override
            public void parseBlastOutput(String output) {
                out.println(output);
            }            
        });
        get.setFormatType("Text");
        Blast blast = new Blast(put, get);
        blast.run();
        out.close();
    }
    
    /*
     * Replaces those funky nucleotides with Ns
     */
    public static String replaceNucs(String sequence) {
        sequence = sequence.replaceAll("M", "N");
        sequence = sequence.replaceAll("K", "N");
        sequence = sequence.replaceAll("R", "N");
        sequence = sequence.replaceAll("S", "N");
        sequence = sequence.replaceAll("W", "N");
        sequence = sequence.replaceAll("Y", "N");
        sequence = sequence.replaceAll("B", "N");
        sequence = sequence.replaceAll("D", "N");
        sequence = sequence.replaceAll("H", "N");
        sequence = sequence.replaceAll("V", "N");
        sequence = sequence.toUpperCase();
        return sequence;
    }
    
    public static void RunOurJeans() throws Exception{
        //This chunk imports files and replaces weird nucleotides
        System.out.println("Importing file");
        File file;
        file = new File(file1);
        Vector<String> storage = new Vector<String>();
        sequencesX = new Vector<Genes>();
        sequencesY = new Vector<Genes>();
        if (seqXname.length()<1) {
        	seqXname = file.getName().substring(0, file.getName().lastIndexOf("."));
        }
        InputStream fis = null;
        try {
            fis = new FileInputStream(file);
        } catch (FileNotFoundException e1) {
            JOptionPane.showMessageDialog(null, "Sequence1 not found. Please rerun SPlot");
            System.exit(0);
        }
        BufferedReader reader = new BufferedReader(new InputStreamReader(fis));
        String line;
        while ((line = reader.readLine()) != null) {
            storage.add(line.trim());
        }
        String sequence = "";
        String id = "";
        for (int i = 0; i<storage.size();) {
            if (storage.get(i).contains(">")) {
                id = storage.get(i);
                i++;
                try {
                    while (!storage.get(i).contains(">")){
                        sequence += storage.get(i);
                        if (sequence.length()>100000) {
                            JOptionPane.showMessageDialog(null, "Seems like you're running a genome in Gene-Mode. Please re-run SPlot.");
                            System.exit(0);
                        }
                        i++;
                    }
                } catch(ArrayIndexOutOfBoundsException exception) {
                    sequencesX.add(new Genes(id, sequence));
                    sequence = "";
                    id = "";
                    break;
                }
            	if (sequence.contains("U")) {
            		sequence.replaceAll("U", "T");
            	}
                //System.out.println("1: " + i);
                sequence = replaceNucs(sequence);
                sequencesX.add(new Genes(id, sequence));
                sequence = "";
                id = "";
            }
        }
       // System.out.println("Pie");
        reader.close();
        fis.close();
        
        Vector<String> storage2 = new Vector<String>();
        File fileTwo;
        fileTwo = new File(file2);
        InputStream is = null;
        try {
            is = new FileInputStream(fileTwo);
        } catch (FileNotFoundException e1) {
            JOptionPane.showMessageDialog(null, "Sequence2 not found. Please rerun SPlot");
            System.exit(0);
        }
        BufferedReader reader2 = new BufferedReader(new InputStreamReader(is));
        while ((line = reader2.readLine()) != null) {
            storage2.add(line.trim());
        }
        sequence = "";
        id = "";
        for (int i = 0; i<storage.size();) {
            if (storage2.get(i).contains(">")) {
                id = storage2.get(i);
                i++;
                try {
                    while (!storage2.get(i).contains(">")){
                        sequence += storage2.get(i);
                        if (sequence.length()>100000) {
                            JOptionPane.showMessageDialog(null, "Seems like you're running a genome in Gene-Mode. Please re-run SPlot.");
                            System.exit(0);
                        }
                        i++;
                    }
                } catch(ArrayIndexOutOfBoundsException exception) {

                    sequencesY.add(new Genes(id, sequence));
                    sequence = "";
                    id = "";
                    break;
                }
                //System.out.println("2: " + i);
            	if (sequence.contains("U")) {
            		sequence.replaceAll("U", "T");
            	}
                sequence = replaceNucs(sequence);
                sequencesY.add(new Genes(id, sequence));
                sequence = "";
                id = "";
            }
        }
      if (seqYname.length()<1) {
      	seqYname = fileTwo.getName().substring(0, fileTwo.getName().lastIndexOf("."));
      }
      reader2.close();
      is.close();
      int intNumGenesX = sequencesX.size();
      int intNumGenesY = sequencesY.size();

        //calculate bitshifting backwards/forwards numbers
        numShifts = 64 - (2*(nmer-1));
        numShiftsMinus = numShifts - 2;
        
        
        //Vectors and such to be used for autoID of RUCPS
        maxValues_X = new double[intNumGenesX];
        maxValues_Y = new double[intNumGenesY];
        averageValues_X = new double[intNumGenesX];
        averageValues_Y = new double[intNumGenesY];
        
        // storage for kmerCounts
        xStorage = new Vector<float[]>();
        yStorage = new Vector<float[]>();

        /*
         * We start filtering things here
         */
        int sum = 0;
        badX = new Vector<Integer>();
        badY = new Vector<Integer>();
        
        
        //adds windows that do not have many Rs (less than 20% Rs, to be processed later)
        for (int iii = 0; iii<intNumGenesX; iii++) {
            sum = 0;
            int len = sequencesX.get(iii).sequence.length();
            float[] toAdd = runGetKmers(sequencesX.get(iii).sequence);
            for (int i = 0; i<toAdd.length; i++) {
                sum += (int)toAdd[i];
            }
            if (sum < len * .8) {
            	badX.add(iii);
            	//System.out.println("this is bad in X: " + iii + "," + sum);
            }
            for (int i = 0; i<toAdd.length; i++) {
                toAdd[i] = toAdd[i]/(float)len;								//Add normalized kmer compositions
            }
            xStorage.add(toAdd);
        }
        System.out.println("xStorage Done");
        
        //adds windows that do not have many Rs (less than 20% Rs, to be processed later)
        for (int iii = 0; iii<intNumGenesY; iii++) {
            sum = 0;
            int len = sequencesY.get(iii).sequence.length();
            float[] toAdd = runGetKmers(sequencesY.get(iii).sequence);
            for (int i = 0; i<toAdd.length; i++) {
                sum += (int)toAdd[i];
            }
            if (sum < len * .8) {
            	badX.add(iii);
            	//System.out.println("this is bad in Y: " + iii + "," + sum);
            }
            for (int i = 0; i<toAdd.length; i++) {
                toAdd[i] = toAdd[i]/(float)len;
            }
            yStorage.add(toAdd);
        }
        
        regressionValues = new float[intNumGenesX][intNumGenesY];
        //So Glimpse doesn't hurt itself when it tries to make heatmap-> Glimpse size needs to be defined
        numberIterationsX = intNumGenesX;
        numberIterationsY = intNumGenesY;
        //Begin the bit shifting magic
        final int intIterX = intNumGenesX;
        final int intIterY = intNumGenesY;
        
        SwingWorker aWorker = new SwingWorker() 
        {
        	 private JFrame frame = new JFrame();
        	 private JDialog dialog = new JDialog(frame, "SPlot Loading...", true);
        	 private JProgressBar progressBar = new JProgressBar();
        	  
		@Override
		protected Object doInBackground() throws Exception {
			progressBar.setString("SPlot Loading...");
 		    progressBar.setStringPainted(true);
 		    progressBar.setPreferredSize(new Dimension(600, 100));
 		    dialog.getContentPane().add(progressBar);
 		    dialog.pack();
 		    dialog.setModal( false );
 		    dialog.setVisible(true);
 		    progressBar.setMaximum(intIterX);
 	        System.out.println("Begin Bit Shifting");
 	        SpearmansCorrelation correlationMaker = new SpearmansCorrelation();
 	        
 	        for (int ii = 0; ii<intNumGenesX; ii++) {
 	            for (int jj = 0; jj<intNumGenesY; jj++) {
 	                //for each y and x, get a regression value and store it in our 2D array
 	                float rvalue;
 	            	if (spearman == true) {
 	                	rvalue = (float)correlationMaker.correlation( convertFloatsToDoubles(xStorage.get(ii)) , convertFloatsToDoubles(yStorage.get(jj)));
 	                }		
 	                else {
 	                	rvalue = getR(xStorage.get(ii), yStorage.get(jj));
 	                }
 	                regressionValues[ii][jj] = rvalue;
 	                //tracks lows/highs/averages for automatic analysis
 	                if(rvalue > maxValues_X[ii]){
 	                    maxValues_X[ii] = rvalue;
 	                }
 	                if(rvalue > maxValues_Y[jj]){
 	                    maxValues_Y[jj] = rvalue;
 	                }
 	                averageValues_X[ii] += rvalue;
 	                averageValues_Y[jj] += rvalue;
 	            }
                progressBar.setValue(ii);
            	if (ii%100==0) {
            		System.out.println(ii + "/" + intIterX);
            	}

 	        }       
			frame.dispose();
			return null;
		}
        };
        aWorker.execute();
        aWorker.get();
        System.out.println("Regression Completed");
        
        for (int zz = 0; zz<badX.size(); zz++) {
        	for (int jj = 0; jj<intNumGenesY; jj++) {
        		regressionValues[badX.get(zz)][jj] = 5;
        	}
        }
        
        for (int zz=0; zz < badY.size(); zz++) {
        	for (int jj = 0; jj<intNumGenesX; jj++) {
        		regressionValues[jj][badY.get(zz)] = 5;
        	}
        }

        //initiate graphing
        System.out.println("Graphing");
        Example.showWithSwing(a);
        
        //run GUI in separate thread
        Thread z = new Thread (new Runnable() {
            @Override
            public void run() {
                try {
                    BackGUI processing = new BackGUI();
                    processing.runGUI();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        });
        z.start();        
    }
    /*
     * This is where stuff happens
     */
    public static void RunOurBaby() throws Exception{
        //This chunk imports files and replaces weird nucleotides
        System.out.println("Importing file");
        
        File file;
        file = new File(file1);
        System.out.println("file1: "+ file1);
        System.out.println("file2: "+ file2);
        if (seqXname.length()<1) {
        	seqXname = file.getName().substring(0, file.getName().lastIndexOf("."));
        }
        InputStream fis = null;
        BufferedReader reader = null;
        try {
        	if(byExample == false){ 
        		fis = new FileInputStream(file);
        		reader = new BufferedReader(new InputStreamReader(fis));}
        	else {
        		fis = new FileInputStream("NC_022806.1_P.a_PA1R.fna");
        		reader = new BufferedReader(new InputStreamReader(fis));}
        } catch (FileNotFoundException e1) {
            JOptionPane.showMessageDialog(null, "Sequence1 not found. Please rerun SPlot");
            System.exit(0);
        }
        
        StringBuilder out = new StringBuilder();
        String line;
        line = reader.readLine();		
        if(!line.contains(">")) out.append(line.trim());
        else System.out.println("this is a fasta: " + line);
        		//TODO make a checker for fasta format, if they are idiots we dont want crashes
        while ((line = reader.readLine()) != null) {
            out.append(line.trim());
        }
        x = out.toString();   							//Prints the string content read from input stream
        reader.close();
        fis.close();
        x = replaceNucs(x);
        
        if (seqYname.length()<1) {
          	seqYname = file.getName().substring(0, file.getName().lastIndexOf("."));
          }
        
        
        File yfile;
        yfile = new File(file2);
        if (seqYname.length()<1) {
        	seqYname = yfile.getName().substring(0, yfile.getName().lastIndexOf("."));
        }
        InputStream yfis = null;
        BufferedReader yreader = null;
        try {
        	if(byExample == false){ 
        		yfis = new FileInputStream(yfile);
        		yreader = new BufferedReader(new InputStreamReader(yfis));}
        	else {
        		yfis = new FileInputStream("NC_009656.1_P.a_PA7.fna");
        		yreader = new BufferedReader(new InputStreamReader(yfis));}
        } catch (FileNotFoundException e1) {
            JOptionPane.showMessageDialog(null, "Sequence2 not found. Please rerun SPlot");
            System.exit(0);
        }
        
        StringBuilder yout = new StringBuilder();
        String yline;
        yline = yreader.readLine();		
        if(!yline.contains(">")) yout.append(yline.trim());
        else System.out.println("this is a fasta: " + yline);
        		//TODO make a checker for fasta format, if they are idiots we dont want crashes
        while ((yline = yreader.readLine()) != null) {
            yout.append(yline.trim());
        }
        y = yout.toString();   							//Prints the string content read from input stream
        yreader.close();
        y = replaceNucs(y);
        if (seqYname.length()<1) {
          	seqYname = yfile.getName().substring(0, yfile.getName().lastIndexOf("."));
        }
        /*
        File fileTwo;
        fileTwo = new File(file2);
        InputStream is = null;
        BufferedReader reader2 = null;
        try {
        	
        	if(byExample == false){ 
        		is = new FileInputStream(fileTwo); 
        	}
        	else {
        		is = Main.class.getResourceAsStream("/com/Pa_PA7.txt"); 
        	}
        } catch (FileNotFoundException e1) {
            JOptionPane.showMessageDialog(null, "Sequence2 not found. Please rerun SPlot");
            System.exit(0);
        }
        reader2 = new BufferedReader(new InputStreamReader(is));
        StringBuilder out2 = new StringBuilder();
        String line2;
        
        line2 = reader2.readLine();		
        if(!line2.contains(">")) out2.append(line2.trim());
        else System.out.println("this is a fasta: " + line2);
        
        while ((line2 = reader2.readLine()) != null) {
            out2.append(line2);
        }
        y = out2.toString();
        reader2.close();
        y = replaceNucs(y);
        */
        
        //Some prints to know how long are the things we are dealing with
        System.out.println("x.length: "+x.length());
        System.out.println("y.length: "+y.length());
        
        //calculate number of iterations for loops
        numberIterationsX = (x.length()/slidingSize) - (int)Math.ceil((windowLength/slidingSize)-1);
        int intIterationsX = (int)numberIterationsX;
        numberIterationsY = (y.length()/slidingSize) - (int)Math.ceil((windowLength/slidingSize)-1);
        System.out.println(numberIterationsY);
        int intIterationsY = (int)numberIterationsY;
        //calculate bitshifting backwards/forwards numbers
        numShifts = 64 - (2*(nmer-1));
        numShiftsMinus = numShifts - 2;
        
        //Vectors and such to be used for autoID of RUCPS
        maxValues_X = new double[intIterationsX];
        maxValues_Y = new double[intIterationsY];
        averageValues_X = new double[intIterationsX];
        averageValues_Y = new double[intIterationsY];
        Vector<float[]> xStorage = new Vector<float[]>();
        Vector<float[]> yStorage = new Vector<float[]>();

        /*
         * We start filtering things here
         */
        
        int sum = 0;
        
        intIterationsX--;
        intIterationsY--;//naive but i don't care
        float[] toAddX;
        float[] toAddY;
        badX = new Vector<Integer>();
        
        //adds windows that do not have many Rs (less than 20% Rs, to be processed later)
        for (int iii = 0; iii<intIterationsX; iii++) {
            sum = 0;
            toAddX = runGetKmers(x.substring(slidingSize*iii, (slidingSize*(iii))+windowLength));
            
            /*try {
                toAdd = runGetKmers(x.substring(slidingSize*iii, (slidingSize*(iii))+windowLength)); //Switched slidingSize with windowLength, for namesake
                } catch (IndexOutOfBoundsException e) {
                    break;
                }
            */
            for (int i = 0; i<toAddX.length; i++) {
                sum += (int)toAddX[i];
            }
            if (sum < windowLength * .8) {
                badX.add(iii);
                //System.out.println("this is bad in X: " + iii);
                
            }
            xStorage.add(toAddX);//delete this to take out bad windows... added this recently

        }
        System.out.println("xStorage complete");
        //adds windows that do not have many Rs (less than 20% Rs, to be processed later)
        
        badY = new Vector<Integer>();
        for (int iiiy = 0; iiiy<intIterationsY; iiiy++) {
            sum = 0;
            toAddY = runGetKmers(y.substring(slidingSize*iiiy, slidingSize*(iiiy) + windowLength));
            
            /*try {
              toAdd = runGetKmers(y.substring(slidingSize*iii, slidingSize*(iii) + windowLength));
            } catch (IndexOutOfBoundsException e) {
                break;
            }*/
            
            for (int i = 0; i<toAddY.length; i++) {
                sum += (int)toAddY[i];
            }
            if (sum < windowLength * .8) {
                badY.add(iiiy);
               // System.out.println("this is bad in Y: " + iiiy);
                
            }
            yStorage.add(toAddY);
        }
        System.out.println("yStorage complete");
        
        //System.out.println(badX);
        //System.out.println(badY);
        //adjust number of iterations after throwing out bad windows
        //TODO dont throw out windows, but rather paint them white
        
        
        //badX and badY became holders for which windows are bad, but the windows still exist
        //intIterationsX = intIterationsX - badX;
        //intIterationsY = intIterationsY - badY;
        numberIterationsX = (double)intIterationsX;
        numberIterationsY = (double)intIterationsY;
        try {
        	regressionValues = new float[intIterationsX][intIterationsY];
        }
        catch (OutOfMemoryError e) {
            JOptionPane.showMessageDialog(null, "Out of Memory. Please Rerun SPlot with a larger window size.");
            System.exit(0);
        }
        final int intIterY = intIterationsY;
        final int intIterX = intIterationsX;
        //Begin the bit shifting magic

        SwingWorker aWorker = new SwingWorker() 
        {
        	 private JFrame frame = new JFrame();
        	 private JDialog dialog = new JDialog(frame, "SPlot Loading...", true);
        	 private JProgressBar progressBar = new JProgressBar();
        	  
		@Override
		protected Object doInBackground() throws Exception {
			progressBar.setString("SPlot Loading...");
 		    progressBar.setStringPainted(true);
 		    progressBar.setPreferredSize(new Dimension(600, 100));
 		    dialog.getContentPane().add(progressBar);
 		    dialog.pack();
 		    dialog.setModal( false );
 		    dialog.setVisible(true);
 		    progressBar.setMaximum(intIterX);
 	        SpearmansCorrelation correlationMaker = new SpearmansCorrelation();

			for (int ii = 0; ii<intIterX; ii++) {
                for (int jj = 0; jj<intIterY; jj++) {
                   	float rvalue;
                	if (spearman == true) {
                        rvalue =  (float)correlationMaker.correlation( convertFloatsToDoubles(xStorage.get(ii)) , convertFloatsToDoubles(yStorage.get(jj)));
                   	}
                    //for each y and x, get a regression value and store it in our 2D array
                	else {
                		rvalue= getR(xStorage.get(ii), yStorage.get(jj));
                	}
                    
                    regressionValues[ii][jj] = rvalue;
                    //tracks lows/highs/averages for automatic analysis
                    if(rvalue > maxValues_X[ii]){
                        maxValues_X[ii] = rvalue;
                    }
                    if(rvalue > maxValues_Y[jj]){
                        maxValues_Y[jj] = rvalue;
                    }
                    averageValues_X[ii] += rvalue;
                    averageValues_Y[jj] += rvalue;
                
                	
                }
                progressBar.setValue(ii);
            	if (ii%100==0) {
            		System.out.println(ii + "/" + intIterX);
            	}
            	} 
			frame.dispose();
			return null;
		}
        };
        aWorker.execute();
        aWorker.get();

        
 
        //this makes the correlational values of the rows and columns of badX and badY all equal to 5, signifying they are bad rows
        for(int zz = 0; zz < badX.size(); zz++){
        	for (int jj = 0; jj<intIterationsY; jj++) {
        		regressionValues[badX.get(zz)][jj] = 5;
        	}
        }
        for(int zz = 0; zz < badY.size(); zz++){
        	for (int jj = 0; jj<intIterationsX; jj++) {
        		regressionValues[jj][badY.get(zz)] = 5;
        	}
        }
        //
        
        
        System.out.println("Regression Completed");
        

        //initiate graphing
        System.out.println("Graphing");
        
        Example.showWithSwing(a);
        
        //run back end GUI in separate thread
        Thread z = new Thread (new Runnable() {
            @Override
            public void run() {
                try {
                    //System.out.println("0");
                    BackGUI processing = new BackGUI();
                    processing.runGUI();
                    //System.out.println("1");
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        });
        z.start();        
        }
        
 
    
    
    
    //Glimpse code to take a picture of our pretty SPlot
    public static void takePicture (HeatMapExample a, String name) throws Exception {
       // a.cursorPainter.setVisible(false);
        GLProfile glProfile = GLUtils.getDefaultGLProfile( );
        
        // generate a GLContext by constructing a small offscreen framebuffer
        final GLOffscreenAutoDrawable glDrawable = GLUtils.newOffscreenDrawable( glProfile );

        // create an offscreen GlimpseCanvas which shares an OpenGL context with the above drawable
        // (its size is 1000 by 1000 pixels)
        final FBOGlimpseCanvas canvas = new FBOGlimpseCanvas( glDrawable.getContext( ), 1000, 1000);

        // set the Glimpse look and feed of the canvas just like we would for an onscreen canvas
        canvas.setLookAndFeel( new SwingLookAndFeel( ) );
        
        // use one of the previous examples to build a simple plot to draw
        ColorAxisPlot2D layout = a.getLayout( );
        // add the layout to the offscreen canvas
        layout.getCrosshairPainter().setVisible(false);
        
        canvas.addLayout( layout );
        // draw the canvas to a BufferedImage and write the image to a file
        //        BufferedImage image = canvas.toBufferedImage( );


        
        //testing
        ImageType imageType = ImageType.TIFF;
        FileOutputStream fo = new FileOutputStream(name + "." + imageType.getExtension());
        ImageWriter writer = ImageIO.getWriter(imageType);
        ImageParam.ImageParamBuilder builder = ImageParam.getBuilder();
        TIFFOptions tiffOptions = new TIFFOptions();
        tiffOptions.setApplyPredictor(true);
        tiffOptions.setTiffCompression(Compression.JPG);
        tiffOptions.setJPEGQuality(60);
        tiffOptions.setPhotoMetric(PhotoMetric.SEPARATED);
        tiffOptions.setWriteICCProfile(true);
        builder.imageOptions(tiffOptions);
        writer.setImageParam(builder.build());
        ImageIO.write(canvas.toBufferedImage(), fo, imageType);
        fo.close();

        // File myNewTIFF_File =  new File(name + ".tiff");
 /*
        /ImageIO.write(image, "TIF", myNewTIFF_File);
        //ImageIO.write( image, "TIFF", new File( name + ".tif" ) );
        ImageIO.write( image, "PNG", new File( name + ".png" ) );
*/
        //a.cursorPainter.setVisible(true);

    }
    
    //Returns the Pearson Correlation coefficient
    public static float getR (float[] x, float[] y) {
        float xAverage = 0;
        float yAverage = 0;
        for (int i =0; i<x.length; i++) {
            xAverage += x[i];
            yAverage += y[i];
        }
        xAverage = xAverage/x.length;
        yAverage = yAverage/y.length;
        float xy = 0;
        float xSq = 0;
        float ySq = 0;
        for (int i = 0; i<x.length; i++) {
            x[i] = (x[i] - xAverage);
            y[i] = (y[i] - yAverage);
            xy += (x[i]*y[i]);
            xSq += x[i]*x[i];
            ySq += y[i]*y[i];
        }
        return (float) (xy/Math.sqrt(xSq*ySq));
    }

//This method is the runner for each windows' bitshifts    
    public static float[] getKmers(String sequence, float[] kmerComp) {
        //initialize first set of kmers
        Long temp = null;
        Long full = Long.parseUnsignedLong("0");
        int i = 0;
        for (i=0; i<nmer; i++) {
            temp  = nucToNum(sequence.charAt(i));
            try {
            	full = full + temp;
            } catch (NullPointerException e) {
                if (byGene == true) {
                	JOptionPane.showMessageDialog(null, "It seems like there are some invalid characters in your input for by-Gene mode... Please check your input files");
                 }
                 if (byGene == false) {
                 	JOptionPane.showMessageDialog(null, "It seems like there are some invalid characters in your input for Genome mode... Please check your input files");
                 }
               System.exit(0);
            }
            if (i<nmer-1) {
                full = full<<2;
            }
        }
        //add it and its reverse kmer to count array
        kmerComp[full.intValue()] += 1;
        kmerComp[reverser(full)] += 1;

        //delete first nucleotide and add to the end of it
        //add it and its reverse complement to count array
        while (i<sequence.length()) {
            temp = nucToNum(sequence.charAt(i));
            full = fancyShift(full);
            try {
	            full = full + temp;
	            kmerComp[full.intValue()] += 1;
	            kmerComp[reverser(full)] += 1;
	            i++;
            }
            catch(NullPointerException e) {
            	break;
            }
        }
        //return kmerComposition array
        return kmerComp;
    }
    
    //shifts kmer off the cliff and back
    public static Long fancyShift (Long a) {
        a = a << numShifts;
        a = a >>> numShiftsMinus;
        return a;
    }
    
    //reverses kmer and makes the reverse complement
    public static int reverser(Long xx) {
        xx = ~xx;
        xx = Long.reverse(xx);
        xx = xx >>> numShiftsMinus;
        String rc = Long.toBinaryString(xx);
        int length = rc.length();
        if (length%2==1) {
            rc = '0' + rc;
        }
        char[] twosies = rc.toCharArray();
        String newString = "";
        for (int i = 0; i<twosies.length; i=i+2) {
            if (twosies[i] == '0' && twosies[i+1] == '1') {
                twosies[i] = '1';
                twosies[i+1] = '0';
            }
            else if (twosies[i] == '1' && twosies[i+1] == '0') {
                twosies[i] = '0';
                twosies[i+1] = '1';
            }
            newString += twosies[i];
            newString += twosies[i+1];
        }
        Long l = parseLong(newString, 2);
        return l.intValue();
    }
    
    private static long parseLong(String s, int base) {
        return new BigInteger(s, base).longValue();
    }
    //Silly initializations for optimization
    public static Long aLong = Long.parseUnsignedLong("0");
    public static Long cLong = Long.parseUnsignedLong("1");
    public static Long gLong = Long.parseUnsignedLong("2");
    public static Long tLong = Long.parseUnsignedLong("3");

    //turns a nucleotide character to a binary representation
    public static Long nucToNum (char a) {
        switch (a) {
            case 'A':
                return aLong;
            case 'C':
                return cLong;
            case 'G':
                return gLong;
            case 'T':
                return tLong;
            default:
                return null;
        }
    }
    
    public static double[] convertFloatsToDoubles(float[] input)
    {
        if (input == null)
        {
            return null; // Or throw an exception - your choice
        }
        double[] output = new double[input.length];
        for (int i = 0; i < input.length; i++)
        {
            output[i] = input[i];
        }
        return output;
    }


}
