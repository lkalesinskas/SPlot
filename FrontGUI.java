import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;

import javax.swing.AbstractAction;
import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTextField; 
 
public class FrontGUI extends Frame implements ActionListener{ 
    public JFrame frame; 
    public JTextField slidingText; 
    public JTextField windowText; 
    public JTextField kmerText; 
    public JTextField firstSeqBox; 
    public JTextField secondSeqBox; 
    public JTextField seqPath2; 
    public JTextField seqPath1; 
    public ButtonGroup fastaGroup = new ButtonGroup(); 
    public JRadioButton byGenomeButton; 
    public JRadioButton byGeneButton;
    public JRadioButton exampleButton;
    public ButtonGroup labelOnOffGroup = new ButtonGroup();
    public JRadioButton labelOn;
    public JRadioButton labelOff;
    public JComboBox<String> cb;
     
      public void placeComponents(JPanel panel) throws Exception { 
            panel.setLayout(null); 
             
            //image stuff 
            /*
            BufferedImage logo = ImageIO.read("SPlotLogo_3.png");
            Image logo2 = logo.getScaledInstance(450,150,Image.SCALE_DEFAULT);//way too big at first 
            JLabel logoPic = new JLabel(new ImageIcon(logo2)); 
            logoPic.setBounds(10, 10, 450, 150); 
            panel.add(logoPic); 
			*/
            /*
            BufferedImage logo = null; 
            logo = ImageIO.read(getClass().getResource("/SPlotLogo_3.png"));   
            
            Image logo2 = logo.getScaledInstance(450,150,Image.SCALE_DEFAULT);//way too big at first 
            JLabel logoPic = new JLabel(new ImageIcon(logo2)); 
            logoPic.setBounds(10, 10, 450, 150); 
            panel.add(logoPic); 
            */
            
            JLabel title = new JLabel("S-plot2"); 
            title.setBounds(30, 10, 450, 150); 
            title.setFont(title.getFont().deriveFont(45.0f)); 
            panel.add(title);
                        
            JLabel briefEx = new JLabel("<html>Welcome to SPlot!<br>SPlot is a program that visualizes similarities and differences between sequences based upon k-mer compositions. A sliding window approach is used. </html>"); 
            briefEx.setBounds(20,150,460,80); 
            briefEx.setFont(briefEx.getFont().deriveFont(12.0f)); 
            panel.add(briefEx); 
             
            JLabel description1 = new JLabel("<html>SPlot has two modes for comparison. SPlot by Genome compares single FASTA genome/chromosome files using a uniform window size. SPlot by Gene requires a multi-FASTA file in which each individual sequence is considered an independent window.</html>"); 
            description1.setBounds(20,235,460,80); 
            description1.setFont(description1.getFont().deriveFont(12.0f)); 
            panel.add(description1); 
             
            byGenomeButton = new JRadioButton(); 
            byGenomeButton.setActionCommand("byGene"); 
            byGenomeButton.setText("SPlot by Genome"); 
            byGenomeButton.setBounds(20,320,150,25); 
            byGenomeButton.setVisible(true); 
            //byGenomeButton.setSelected(true); 
            panel.add(byGenomeButton); 
            fastaGroup.add(byGenomeButton); 
             
            byGeneButton = new JRadioButton(); 
            byGeneButton.setActionCommand("byGenome"); 
            byGeneButton.setText("SPlot by Gene"); 
            byGeneButton.setBounds(200,320,150,25); 
            byGeneButton.setVisible(true); 
            byGeneButton.setVisible(true); 
            panel.add(byGeneButton); 
            fastaGroup.add(byGeneButton); 
             
            JLabel kmerLabel = new JLabel("k-mer Size:"); //here we initiate the different pieces of the front end gui 
            kmerLabel.setBounds(342,550,150,25);            //these are just the text fields 
            panel.add(kmerLabel); 
             
            kmerText = new JTextField(20); 
            kmerText.setBounds(430,550,50,25); 
            panel.add(kmerText); 
             
            JLabel description3 = new JLabel("<html>Parameters:</html>"); 
            description3.setBounds(20,490,460,80); 
            description3.setFont(description3.getFont().deriveFont(16.0f)); 
            panel.add(description3); 
 
            JLabel slidingLabel = new JLabel("Offset:"); 
            slidingLabel.setBounds(368,640,150,25); 
            panel.add(slidingLabel); 
 
            slidingText = new JTextField(20); 
            slidingText.setBounds(430,640,50,25); 
            panel.add(slidingText); 
             
            JLabel windowLabel = new JLabel("Window Size:"); 
            windowLabel.setBounds(330,610,150,25); 
            panel.add(windowLabel); 
 
            windowText = new JTextField(20); 
            windowText.setBounds(430,610,50,25); 
            panel.add(windowText); 
                 
            JLabel firstSeq = new JLabel("X-Axis Sequence Label:"); 
            firstSeq.setBounds(20,350,150,25); 
            panel.add(firstSeq); 
 
            firstSeqBox = new JTextField(20); 
            firstSeqBox.setBounds(180,350,300,25); 
            panel.add(firstSeqBox); 
             
            JLabel secondSeq = new JLabel("Y-Axis Sequence Label:"); 
            secondSeq.setBounds(20,380,150,25); 
            panel.add(secondSeq); 
 
            secondSeqBox = new JTextField(20); 
            secondSeqBox.setBounds(180,380,300,25); 
            panel.add(secondSeqBox); 
 
            JButton button = new JButton("SUBMIT");    //the submit button takes in all of the information inputed into the gui and starts the program 
            button.setBounds(175, 750, 150, 40); 
            panel.add(button); 
            button.addActionListener(this); 
             
            JLabel xdescription = new JLabel("<html>X-Axis Sequence:</html>"); 
            xdescription.setBounds(20,415,280,25); 
            xdescription.setFont(xdescription.getFont().deriveFont(12.0f)); 
            panel.add(xdescription);
            
            
            
            seqPath1 = new JTextField(20); 
            seqPath1.setBounds(200,440,280,25); 
            panel.add(seqPath1); 
             
            JButton firstSequenceButton = new JButton(new AbstractAction("Browse...") {    //this directs you to a window to choose the file location 
                public void actionPerformed(ActionEvent e) { 
                    JFileChooser choosed = new JFileChooser(); 
                    File workingDirectory = new File(System.getProperty("user.dir"));
                    choosed.setCurrentDirectory(workingDirectory);
                    int returnValue = choosed.showOpenDialog(null); 
                    if (returnValue == JFileChooser.APPROVE_OPTION) { 
                         seqPath1.setText(choosed.getSelectedFile().getPath()); 
                    } 
                } 
            }); 
            firstSequenceButton.setBounds(10, 440, 180, 25); 
            panel.add(firstSequenceButton); 
            
            JLabel ydescription = new JLabel("<html>Y-Axis Sequence:</html>"); 
            ydescription.setBounds(20,465,280,25); 
            ydescription.setFont(xdescription.getFont().deriveFont(12.0f)); 
            panel.add(ydescription);
            
            
            
            
            seqPath2 = new JTextField(20); 
            seqPath2.setBounds(200,490,280,25); 
            panel.add(seqPath2); 
             
            JButton SecondSequenceButton = new JButton(new AbstractAction("Browse...") { 
                public void actionPerformed(ActionEvent e) { 
                    JFileChooser chooser = new JFileChooser(); 
                    File workingDirectory = new File(System.getProperty("user.dir"));
                    chooser.setCurrentDirectory(workingDirectory);
                    int returnVal = chooser.showOpenDialog(null); 
                    if (returnVal == JFileChooser.APPROVE_OPTION) { 
                         seqPath2.setText(chooser.getSelectedFile().getPath()); 
                    } 
                } 
            }); 
            SecondSequenceButton.setBounds(10, 490, 180, 25); 
            panel.add(SecondSequenceButton); 
             
            JLabel description2 = new JLabel("<html> SPlot by Genome Parameters: </html>"); 
            description2.setBounds(20,560,460,80); 
            description2.setFont(description2.getFont().deriveFont(12.0f)); 
            panel.add(description2); 
            
            //here is for the example run
            JLabel exampleLabel = new JLabel("<html>Run Sample Data</html>");
            exampleLabel.setFont(xdescription.getFont().deriveFont(16.0f)); 
            exampleLabel.setBounds(20,780,460,40); 
            panel.add(exampleLabel); 
            
            //here is for the example run
            JLabel exampleRunInfo = new JLabel("<html>Pseudomonas aeruginosa PA7 (NC_009656) vs. Pseudomonas aeruginosa PA1R (NC_022806)</html>");
            exampleRunInfo.setFont(xdescription.getFont().deriveFont(12.0f)); 
            exampleRunInfo.setBounds(20,810,460,40); 
            panel.add(exampleRunInfo);
            
            
            exampleButton = new JRadioButton(); 
            exampleButton.setActionCommand("Example"); 
            exampleButton.setText("Example"); 
            exampleButton.setBounds(200,840,150,25); 
            exampleButton.setVisible(true); 
            exampleButton.setVisible(true); 
            panel.add(exampleButton); 
            fastaGroup.add(exampleButton);
            
            labelOff = new JRadioButton();
            labelOff.setActionCommand("Off");
            labelOff.setText("Off");
            labelOff.setBounds(185,670,80,25); 
            labelOff.setVisible(true);
            panel.add(labelOff);
            labelOnOffGroup.add(labelOff);
            
            labelOn = new JRadioButton();
            labelOn.setActionCommand("On");
            labelOn.setText("On");
            labelOn.setBounds(275,670,80,25);
            labelOn.setSelected(true);
            labelOn.setVisible(true);
            panel.add(labelOn);
            labelOnOffGroup.add(labelOn);
            
            JLabel labelInfo = new JLabel("<html>Coordinate Label:</html>");
            labelInfo.setFont(xdescription.getFont().deriveFont(12.0f)); 
            labelInfo.setBounds(20,660,250,40); 
            panel.add(labelInfo); 
           
            JLabel lbl = new JLabel("Select Correlation Method:");
            lbl.setFont(xdescription.getFont().deriveFont(12.0f));
            lbl.setBounds(20, 700, 150, 40);
            lbl.setVisible(true);
            panel.add(lbl);
            
            String[] choices = { "Pearson","Spearman"};
            cb = new JComboBox<String>(choices);
            cb.setBounds(200, 710, 100, 20);
            cb.setVisible(true);
            
            panel.add(cb);
        } 
       
      public void geneVsGenomeAction(ActionEvent ae){ 
            String action = ae.getActionCommand(); 
            if (action.equals("byGene")){ 
                System.out.println("byGene"); 
            }else if (action.equals("byGenome")){ 
                System.out.println("byGenome"); 
            }else {
            	System.out.println("example");
            }
        } 
       
        public void runGUI() throws Exception { 
          frame = new JFrame("S-plot2"); 
          JPanel panel = new JPanel(new BorderLayout());
          panel.setBorder(BorderFactory.createLineBorder(Color.red));
          panel.setPreferredSize(new Dimension(600, 1000));
          JScrollPane pane = new JScrollPane(panel, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);   
          pane.getVerticalScrollBar().setUnitIncrement(10);
          frame.add(pane, BorderLayout.CENTER); 
          frame.setSize(600, 1000); 
          frame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE); 
          frame.addWindowListener(new WindowAdapter() {
        	    @Override
        	    public void windowClosing(WindowEvent we)
        	    { 
        	        String ObjButtons[] = {"Yes","No"};
        	        int PromptResult = JOptionPane.showOptionDialog(null,"Are you sure you want to exit?","S-plot2",JOptionPane.DEFAULT_OPTION,JOptionPane.WARNING_MESSAGE,null,ObjButtons,ObjButtons[1]);
        	        if(PromptResult==JOptionPane.YES_OPTION)
        	        {
        	            System.exit(0);
        	        }
        	    }
        	});
          placeComponents(panel); 
          frame.setVisible(true); 
        } 
        public void runMessage (String message) { 
            JOptionPane.showMessageDialog(null, message);        //some error handling to deal with issues before the program is run and it crashes 
        } 
        @Override 
        public void actionPerformed(ActionEvent e) { 
        	
    		if (labelOff.isSelected()) {
        		main.labelControl = false;
        	}
    		if (cb.getSelectedItem().equals("Spearman")) {
    			main.spearman=true;
    		}
    		
        	//Error handling for starts here
        	boolean runnable = true;				//this variable states when the user parameters are runnable, it turns false if an error is found
        	if(exampleButton.isSelected()){			//all of the parameters for the example

        		main.byExample = true;
	        	main.byGene = false;
	            main.slidingSize = 5000;  
	            main.windowLength = 5000; 
	            main.nmer = 6; 
	            main.seqXname = "P. aeruginosa PA1R";  
	            main.seqYname = "P. aeruginosa PA7"; 
	          
	        //error handling for byGene: files and kmer size needed only: names optional
        	}else if (byGeneButton.isSelected()) {
        			
                	main.byGene = true; 
                	if (kmerText.getText().isEmpty() || Integer.parseInt(kmerText.getText()) <=0) { 
    	                runMessage("Please enter kmer Size."); 
    	                runnable = false;
    	            }else if (seqPath1.getText().isEmpty()) { 
    	                runMessage("Please enter First Sequence."); 
    	                runnable = false;
    	            }else if (seqPath2.getText().isEmpty()) { 
    	                runMessage("Please enter Second Sequence."); 
    	                runnable = false;
    	            }
                	if(runnable == true){
    		            String x = seqPath1.getText(); 
    		            main.file1 = x; 
    		            
    		            x = seqPath2.getText(); 
    		            main.file2 = x; 
    		            
    		            x = firstSeqBox.getText(); 
    		            main.seqXname = x; 
    		            
    		            x = kmerText.getText(); 
    		            main.nmer = Integer.parseInt(x); 
    		            
    		            x = secondSeqBox.getText(); 
    		            main.seqYname = x; 
                	}
            //byGenome error handling: sequence files, kmer size, window size, and sliding amount needed: names optional
            }else if (byGenomeButton.isSelected()){ 
	                main.byGene = false; 
	                if (slidingText.getText().isEmpty() || Integer.parseInt(slidingText.getText()) <=0) { 
		                runMessage("Please enter Sliding amount."); 
    	                runnable = false;
	                }else if (windowText.getText().isEmpty() || Integer.parseInt(windowText.getText()) <=0){ 
			            runMessage("please enter Window Size.");
			            runnable = false;
		            }else if (Integer.parseInt(slidingText.getText())>Integer.parseInt(windowText.getText())){
		            	runMessage("Offset cannot be greater than Window Size.");
		            	runnable = false;
		            }else if (kmerText.getText().isEmpty() || Integer.parseInt(kmerText.getText()) <=0) { 
		                runMessage("Please enter k-mer Size."); 
		                runnable = false;
		            }else if (seqPath1.getText().isEmpty()) { 
		                runMessage("Please enter First Sequence."); 
		                runnable = false;
		            }else if (seqPath2.getText().isEmpty()) { 
		                runMessage("Please enter Second Sequence."); 
		                runnable = false;
		            }
		            else if(Math.pow(4, Integer.parseInt(kmerText.getText()))>Integer.parseInt(windowText.getText())) {
		            	runMessage("4^kmer size must be less than the window length");
		            	runnable = false;
		            }
	                if(runnable == true){ 
		            String x = slidingText.getText(); 
		            main.slidingSize = Integer.parseInt(x); 
		            x = windowText.getText(); 
		            main.windowLength = Integer.parseInt(x); 
		            x = kmerText.getText(); 
		            main.nmer = Integer.parseInt(x); 
		            x = seqPath1.getText(); 
		            main.file1 = x; 
		            x = seqPath2.getText(); 
		            main.file2 = x; 
		            x = firstSeqBox.getText(); 
		            main.seqXname = x; 
		            x = secondSeqBox.getText(); 
		            main.seqYname = x; 
		            }
	        }
        	if(runnable == true){
	            frame.dispose(); 
	                Thread t = new Thread (new Runnable() {        //a new thread is opened to start the actual splot code,  
	                    @Override                                //opening the graphics gui within the front end gui isnt liked, hence the new thread 
	                    public void run() { 
	                        try { 
	                            if (main.byGene == false){System.out.println("by genome");main.RunOurBaby();} 
	                            else {System.out.println("by gene");main.RunOurJeans();} 
	                        } catch (Exception e) { 
	                            e.printStackTrace(); 
	                        } 
	                    } 
	                }); 
	                t.start(); 
        	}
        } 
 
}
