import java.awt.Desktop;
import java.awt.Frame;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.imageio.ImageIO;
import javax.swing.AbstractAction;
import javax.swing.ButtonGroup;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPasswordField;
import javax.swing.JProgressBar;
import javax.swing.JRadioButton;
import javax.swing.JTextField;
import javax.swing.filechooser.FileNameExtensionFilter;

import org.apache.log4j.BasicConfigurator;

import com.algosome.eutils.blast.BlastParser;
import com.algosome.eutils.blast.GetCommand;
import com.algosome.eutils.blast.PutCommand;

public class BackGUI extends Frame implements ActionListener{
	public JFrame frame;
	public JTextField windowText;
	public JTextField kmerText;
	public JTextField startboxinput;
	public JTextField endboxinput;
	public JTextField maxboxinput;
	public JTextField pictureboxinput;
	public JTextField averageboxinput;
	public JTextField seqPath2;
	public JTextField seqPath1;
	public ButtonGroup group = new ButtonGroup();
	//public ButtonGroup bigGroup = new ButtonGroup();
	public JRadioButton seq1Button;
	public JRadioButton seq2Button;
	public JRadioButton seq1er;
	public JRadioButton seq2er; 
	

	//TODO make lots of error handling BLEHHHHHHHHH
	//automated cutoffs				DONE
	//start and stop window			DONE
	//picture name					DONE
	//Internet connectivity check
	
	public void runMessage (String message) { 
        JOptionPane.showMessageDialog(null, message);        //some error handling to deal with issues before the program is run and it crashes 
    }
	public boolean containsIllegals(String toExamine) {		//this creates a check for some nasty characters so that they don't show up in a file path
	    Pattern pattern = Pattern.compile("[{}\\[\\]|\"\\^]");
	    Matcher matcher = pattern.matcher(toExamine);
	    return matcher.find();
	}
	  public void placeComponents(JPanel panel) {
	        panel.setLayout(null);
	        
	        JLabel seqSelect = new JLabel("<html>Select Sequence: </html>");
	        seqSelect.setBounds(20,125,460,25);
	        seqSelect.setFont(seqSelect.getFont().deriveFont(12.0f)); 
	        panel.add(seqSelect);

	        seq1Button = new JRadioButton(main.seqXname);
	        seq1Button.setActionCommand(main.seqXname);
	        seq1Button.setText(main.seqXname);
	        seq1Button.setBounds(150,140,250,25);
	        seq1Button.setVisible(true);
	        //seq1Button.setSelected(true);
	        panel.add(seq1Button);
	        group.add(seq1Button);
	        
	        seq2Button = new JRadioButton(main.seqYname);
	        seq2Button.setActionCommand(main.seqYname);
	        seq2Button.setText(main.seqYname);
	        seq2Button.setBounds(150,165,250,25);
	        seq2Button.setVisible(true);
	        seq2Button.setVisible(true);
	        panel.add(seq2Button);
	        group.add(seq2Button);
	        
	        JButton button = new JButton("Write to File");
	        button.setBounds(50, 190, 150, 25);
	        button.setToolTipText("Writes specified windows to file");
	        panel.add(button);
	        button.addActionListener(this);
	        

	        
	        JButton autoBLAST = new JButton(new AbstractAction("BLAST Windows"){

				@Override
				public void actionPerformed(ActionEvent arg0) {
					double maxCor = Double.parseDouble(maxboxinput.getText());
					double aveCor = Double.parseDouble(averageboxinput.getText());
					if(maxCor < 1 && maxCor > 0 && maxCor > aveCor && aveCor > 0 && aveCor < 1){
						//TODO check if they have internet connection
						runMessage("This may take a while depending on internet connection and NCBI server speed");
						try {
							if (seq1er.isSelected()==false && seq2er.isSelected()==false) {
								runMessage("You have not selected a sequence to BLAST.");
							}
							else if (seq1er.isSelected() && seq2er.isSelected() && main.byGene==false) {
								File xFile = getFileTxt("Specify Save files for X-Axis Sequence");
								String xFileName = xFile.getAbsolutePath();
								File yFile = getFileTxt("Specify Save files for Y-Axis Sequence");
								String yFileName = yFile.getAbsolutePath();
								main.automatedBLAST (maxCor, aveCor, true, true, xFileName, yFileName);
							}
							else if (seq1er.isSelected() && seq2er.isSelected()==false && main.byGene==false) {
								File xFile = getFileTxt("Specify Save files for X-Axis Sequence");
								String xFileName = xFile.getAbsolutePath();
								String yFileName = null;
								main.automatedBLAST (maxCor, aveCor, true, false, xFileName, yFileName);
							}
							else if (seq1er.isSelected()==false && seq2er.isSelected() && main.byGene==false) {
								File yFile = getFileTxt("Specify Save files for Y-Axis Sequence");
								String yFileName = yFile.getAbsolutePath();
								String xFileName = null;
								main.automatedBLAST (maxCor, aveCor, false, true, xFileName, yFileName);
							}
							else if (seq1er.isSelected() && seq2er.isSelected() && main.byGene==true) {
								File xFile = getFileTxt("Specify Save files for X-Axis Sequence");
								String xFileName = xFile.getAbsolutePath();
								File yFile = getFileTxt("Specify Save files for Y-Axis Sequence");
								String yFileName = yFile.getAbsolutePath();
								main.automatedBLASTGenes (maxCor, aveCor, true, true, xFileName, yFileName);
							}
							else if (seq1er.isSelected() && seq2er.isSelected()==false && main.byGene==true) {
								File xFile = getFileTxt("Specify Save files for X-Axis Sequence");
								String xFileName = xFile.getAbsolutePath();
								String yFileName = null;
								main.automatedBLASTGenes (maxCor, aveCor, true, false, xFileName, yFileName);
							}
							else if (seq1er.isSelected()==false && seq2er.isSelected() && main.byGene==true) {
								File yFile = getFileTxt("Specify Save files for Y-Axis Sequence");
								String yFileName = yFile.getAbsolutePath();
								String xFileName = null;
								main.automatedBLASTGenes (maxCor, aveCor, false, true, xFileName, yFileName);
							}
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
					else{
						runMessage("Invalid Values: 0 < Average Correlation < Max Correlation < 1");
					}
				}
	        	
	        });
	        
	        autoBLAST.setBounds(250, 350, 150, 25);
	        autoBLAST.setToolTipText("BLASTs all windows that have an average correlational value and max correlational value that are less than specified values");
	        panel.add(autoBLAST);
	        
	        
	        
	        JButton autoWrite2File = new JButton(new AbstractAction("Output Windows"){
				@Override
				public void actionPerformed(ActionEvent arg0) {
					try {
						double maxCor = Double.parseDouble(maxboxinput.getText());
						double aveCor = Double.parseDouble(averageboxinput.getText());
						
						if(maxCor < 1 && maxCor > 0 && maxCor > aveCor && aveCor > 0 && aveCor < 1){
							if (seq1er.isSelected()==false && seq2er.isSelected()==false) {
								runMessage("You have not selected a sequence to BLAST.");
							}
							else if (seq1er.isSelected() && seq2er.isSelected() && main.byGene==false) {
								File xFile = getFileFasta("Specify Save file for X-Axis Sequence");
								String xFileName = xFile.getAbsolutePath();
								File yFile = getFileFasta("Specify Save file for Y-Axis Sequence");
								String yFileName = yFile.getAbsolutePath();
								main.automatedFile(maxCor, aveCor, true, true, xFileName, yFileName);
							}
							else if (seq1er.isSelected() && seq2er.isSelected()==false && main.byGene==false) {
								File xFile = getFileFasta("Specify Save file for X-Axis Sequence");
								String xFileName = xFile.getAbsolutePath();
								String yFileName = null;
								main.automatedFile (maxCor, aveCor, true, false, xFileName, yFileName);
							}
							else if (seq1er.isSelected()==false && seq2er.isSelected() && main.byGene==false) {
								File yFile = getFileFasta("Specify Save file for Y-Axis Sequence");
								String yFileName = yFile.getAbsolutePath();
								String xFileName = null;
								main.automatedFile (maxCor, aveCor, false, true, xFileName, yFileName);
							}
							else if (seq1er.isSelected() && seq2er.isSelected() && main.byGene==true) {
								File xFile = getFileFasta("Specify Save file for X-Axis Sequence");
								String xFileName = xFile.getAbsolutePath();
								File yFile = getFileFasta("Specify Save file for Y-Axis Sequence");
								String yFileName = yFile.getAbsolutePath();
								main.automatedFileByGene(maxCor, aveCor, true, true, xFileName, yFileName);
							}
							else if (seq1er.isSelected() && seq2er.isSelected()==false && main.byGene==true) {
								File xFile = getFileFasta("Specify Save file for X-Axis Sequence");
								String xFileName = xFile.getAbsolutePath();
								String yFileName = null;
								main.automatedFileByGene (maxCor, aveCor, true, false, xFileName, yFileName);
							}
							else if (seq1er.isSelected()==false && seq2er.isSelected() && main.byGene==true) {
								File yFile = getFileFasta("Specify Save file for Y-Axis Sequence");
								String yFileName = yFile.getAbsolutePath();
								String xFileName = null;
								main.automatedFileByGene(maxCor, aveCor, false, true, xFileName, yFileName);
							}
						}else{
							runMessage("These values do not work: 0 < average correlation < max correlation < 1");
						}
					} catch (FileNotFoundException e) {
						e.printStackTrace();
					}
				}
	        });
	        autoWrite2File.setBounds(50, 350, 150, 25);
	        autoWrite2File.setToolTipText("Writes all windows that have an average correlation and max correlation lower than specified values");
	        panel.add(autoWrite2File);
	        
	      //image stuff 
	        /*
            BufferedImage logo = null; 
            try { 
                logo = ImageIO.read(getClass().getResource("/SPlotLogo_3.png"));
            } catch (IOException e1) { 
                e1.printStackTrace(); 
            }      
            Image logo2 = logo.getScaledInstance(150,50,Image.SCALE_DEFAULT);//way too big at first 
            JLabel logoPic = new JLabel(new ImageIcon(logo2)); 
            logoPic.setBounds(170, 10, 150, 50); 
            panel.add(logoPic); 
	        */
	        
	        
	        //Basic Functionality descriptions
	        JLabel descriptiveBoxT3 = new JLabel("<html>Save Plot: </html>");
	        descriptiveBoxT3.setBounds(190,10,460,25);
	        descriptiveBoxT3.setFont(descriptiveBoxT3.getFont().deriveFont(16.0f)); 
	        panel.add(descriptiveBoxT3);
	        
	        JLabel descriptiveBoxT4 = new JLabel("<html><hr>Examine Specific Windows:</html>");
	        descriptiveBoxT4.setBounds(125,60,460,45);
	        descriptiveBoxT4.setFont(descriptiveBoxT3.getFont().deriveFont(16.0f)); 
	        panel.add(descriptiveBoxT4);
	        
	        JLabel startbox = new JLabel("Start window:");
	        startbox.setBounds(20,100,150,25);
	        panel.add(startbox);

	        startboxinput = new JTextField(20);
	        startboxinput.setBounds(110,100,60,25);
	        panel.add(startboxinput);
	        
	        JLabel endbox = new JLabel("Stop window:");
	        endbox.setBounds(200,100,150,25);
	        panel.add(endbox);

	        endboxinput = new JTextField(20);
	        endboxinput.setBounds(300,100,60,25);
	        panel.add(endboxinput);
	        
	        JLabel descriptiveBoxT5 = new JLabel("<html><hr>Windows of Unusual Composition: </html>");
	        descriptiveBoxT5.setBounds(105,205,460,55);
	        descriptiveBoxT5.setFont(descriptiveBoxT3.getFont().deriveFont(16.0f)); 
	        panel.add(descriptiveBoxT5);
	        
	        JLabel badWin = new JLabel("<html><hr>Invalid Windows:</html>");
	        badWin.setBounds(180,365,460,55);
	        badWin.setFont(badWin.getFont().deriveFont(16.0f)); 
	        panel.add(badWin);
	        

	        String badXlist = "<html>Bad windows in "+ main.seqXname+": ";
	        String badYlist = "<html>Bad windows in "+ main.seqYname+": ";
		        for(int x : main.badX){
		        	badXlist += x;
		        	badXlist += ", ";
		        }
		        
		        badXlist = badXlist.substring(0,badXlist.length()-2);
		        badXlist += ".</html>";
		        for(int y : main.badY){
		        	badYlist += y;
		        	badYlist += ", ";
		        }
		        badYlist = badYlist.substring(0,badYlist.length()-2);
		        badYlist += ".</html>";
		        
		        if(!main.badX.isEmpty()){
			        JLabel descriptiveBoxT11 = new JLabel(badXlist);
			        descriptiveBoxT11.setBounds(20,375,220,120);
			        panel.add(descriptiveBoxT11);
		        }
		        if(!main.badY.isEmpty()){
			        JLabel descriptiveBoxT12 = new JLabel(badYlist);
			        descriptiveBoxT12.setBounds(260,375,220,120);
			        panel.add(descriptiveBoxT12);
		        }
	        
	        
	        JButton picturebutton = new JButton(new AbstractAction("Save Image..."){
				@Override
				public void actionPerformed(ActionEvent e) {
					JFileChooser saveFile =new JFileChooser();
					saveFile.setDialogTitle("Specify a file to save");
                    File workingDirectory = new File(System.getProperty("user.dir"));    
                    saveFile.setCurrentDirectory(workingDirectory);
                    saveFile.addChoosableFileFilter(new FileNameExtensionFilter(".tiff", "TIFF Image"));
                    int retrival = saveFile.showSaveDialog(null);
                    if (retrival == JFileChooser.APPROVE_OPTION) {
                        try {
                            String pictureName = saveFile.getSelectedFile().getAbsolutePath();
        					if( containsIllegals(pictureName)){
        						runMessage("illegal character(s) in file name.");		//error handling to not make picture to an illegal file location
        					}
        					else{
        						try {
        							main.takePicture(main.a, pictureName);
        							//main.a.cursorPainter.setVisible(true);
        						} catch (Exception e1) {
        							e1.printStackTrace();
        						}
        					}
                        } catch (Exception ex) {
                            ex.printStackTrace();
                        }
                    }

				}
	        	
	        });
	        picturebutton.setBounds(50, 40, 150, 25);
	        picturebutton.setToolTipText("Saves a .TIFF Image of the SPlot");
	        panel.add(picturebutton);

	        JButton matrixbutton = new JButton(new AbstractAction("Save Matrix..."){
				@Override
				public void actionPerformed(ActionEvent e) {
					JFileChooser saveFile =new JFileChooser();
					saveFile.setDialogTitle("Specify a file to save");
                    File workingDirectory = new File(System.getProperty("user.dir"));    
                    saveFile.setCurrentDirectory(workingDirectory);
                    saveFile.addChoosableFileFilter(new FileNameExtensionFilter(".csv", "CSV File"));
                    int retrival = saveFile.showSaveDialog(null);
                    if (retrival == JFileChooser.APPROVE_OPTION) {
                        try {
                            String matrixFileName = saveFile.getSelectedFile().getAbsolutePath();
        					if( containsIllegals(matrixFileName)){
        						runMessage("illegal character(s) in file name.");		//error handling to not make picture to an illegal file location
        					}
        					else{
        						try {
        							String file = saveFile.getSelectedFile().getAbsolutePath();
        							if (file.contains(".csv")==false) {
        								file+=".csv";
        							}
        							PrintStream out = new PrintStream (new File(file));
        							for (int i = 0; i<main.regressionValues.length; i++) {
        								for (int j = 0; j<main.regressionValues[i].length; j++) {
        									out.print(main.regressionValues[i][j]+",");
        								}
        								out.println();
        							}
        							out.close();
        						} catch (Exception e1) {
        							e1.printStackTrace();
        						}
        					}
                        } catch (Exception ex) {
                            ex.printStackTrace();
                        }
                    }

				}
	        	
	        });
	        matrixbutton.setBounds(250, 40, 150, 25);
	        matrixbutton.setToolTipText("Saves a .CSV of all correlational values between windows");
	        panel.add(matrixbutton);
	        
	        //this is for auto analysis restrictions
	        JLabel maxbox = new JLabel("Max Correlation:");
	        maxbox.setBounds(20,250,150,25);
	        panel.add(maxbox);

	        maxboxinput = new JTextField(20);
	        maxboxinput.setBounds(130,250,60,25);
	        maxboxinput.setText(".9");
	        panel.add(maxboxinput);
	        
	        JLabel averagebox = new JLabel("Average Correlation:");
	        averagebox.setBounds(250,250,150,25);
	        panel.add(averagebox);

	        averageboxinput = new JTextField(20);
	        averageboxinput.setBounds(380,250,60,25);
	        averageboxinput.setText(".4");
	        panel.add(averageboxinput);
	        
	        JLabel seqSelect2 = new JLabel("<html>Select Sequence: </html>");
	        seqSelect2.setBounds(20,275,460,25);
	        seqSelect2.setFont(seqSelect2.getFont().deriveFont(12.0f)); 
	        panel.add(seqSelect2);

	        seq1er = new JRadioButton(main.seqXname);
	        seq1er.setActionCommand(main.seqXname);
	        seq1er.setText(main.seqXname);
	        seq1er.setBounds(150,290,250,25);
	        seq1er.setVisible(true);
	        //seq1Button.setSelected(true);
	        panel.add(seq1er);
	        //bigGroup.add(seq1er);
	        
	        seq2er = new JRadioButton(main.seqYname);
	        seq2er.setActionCommand(main.seqYname);
	        seq2er.setText(main.seqYname);
	        seq2er.setBounds(150,315,250,25);
	        seq2er.setVisible(true);
	        seq2er.setVisible(true);
	        panel.add(seq2er);
	        //bigGroup.add(seq2er);

	        
	        JButton BLASTButton = new JButton(new AbstractAction("BLAST") {
	        	@Override
	        	public void actionPerformed(ActionEvent e) {
	        		if(Desktop.isDesktopSupported())
	        		{
	        		String seq = "";
					int startnum = Integer.parseInt(startboxinput.getText());
					int endnum = Integer.parseInt(endboxinput.getText());
					
					boolean runnable = true;//for error handling

					//error handling of start and stop for BLAST
					if(startnum >= endnum){
						runMessage("The windows specified are out of range.");
						runnable = false;
					}
					if(startnum < 0){
						runMessage("The windows specified are out of range.");
						runnable = false;
					}
					if(seq1Button.isSelected() && endnum > main.numberIterationsX){
						runMessage("The windows specified are out of range.");
						runnable = false;
					}
					if(seq2Button.isSelected() && endnum > main.numberIterationsY){
						runMessage("The windows specified are out of range.");
						runnable = false;
					}
					//
					if(runnable == true){
					if (seq1Button.isSelected() && main.byGene==false) {
						seq = main.x.substring(startnum * main.slidingSize, endnum * main.slidingSize);
						seq = seq.toUpperCase();
						JFileChooser saveFile =new JFileChooser();
						saveFile.setDialogTitle("Specify a file to save");

						File workingDirectory = new File(System.getProperty("user.dir"));    
						saveFile.setCurrentDirectory(workingDirectory);
						saveFile.addChoosableFileFilter(new FileNameExtensionFilter(".txt", "Text File"));
						int retrival = saveFile.showSaveDialog(null);
						if (retrival == JFileChooser.APPROVE_OPTION) {
						  try {
						   String pictureName = saveFile.getSelectedFile().getAbsolutePath();
						   if (pictureName.contains(".txt")==false) {
							   pictureName+=".txt";
						   }
						   try {
							runBlast(seq, pictureName);
							} catch (Exception e1) {
								e1.printStackTrace();
								}
							}
						    catch (Exception ex) {
						        ex.printStackTrace();
						    }
						}
					}
					else if (seq2Button.isSelected() && main.byGene==false) {
						seq = main.y.substring(startnum * main.slidingSize, endnum * main.slidingSize);
						seq = seq.toUpperCase();
						JFileChooser saveFile =new JFileChooser();
						saveFile.setDialogTitle("Specify a file to save");
						File workingDirectory = new File(System.getProperty("user.dir"));    
						saveFile.setCurrentDirectory(workingDirectory);
						saveFile.addChoosableFileFilter(new FileNameExtensionFilter(".txt", "Text File"));
						int retrival = saveFile.showSaveDialog(null);
						if (retrival == JFileChooser.APPROVE_OPTION) {
						  try {
						   String pictureName = saveFile.getSelectedFile().getAbsolutePath();
						   if (pictureName.contains(".txt")==false) {
							   pictureName+=".txt";
						   }
						   try {
							runBlast(seq, pictureName);
							} catch (Exception e1) {
								e1.printStackTrace();
								}
							}
						    catch (Exception ex) {
						        ex.printStackTrace();
						    }
						}
					}
					else if (seq1Button.isSelected() && main.byGene==true) {
						JFileChooser saveFile =new JFileChooser();
						saveFile.setDialogTitle("Specify a file root to save");
						File workingDirectory = new File(System.getProperty("user.dir"));    
						saveFile.setCurrentDirectory(workingDirectory);
						saveFile.addChoosableFileFilter(new FileNameExtensionFilter(".txt", "Text File"));
						int retrival = saveFile.showSaveDialog(null);
						if (retrival == JFileChooser.APPROVE_OPTION) {
						  try {
						   String pictureName = saveFile.getSelectedFile().getAbsolutePath();
						   if (pictureName.contains(".txt")==false) {
							   pictureName+=".txt";
						   }
						   try {
							   for (int i = startnum; i<endnum; i++) {
								   seq = main.sequencesX.get(i).sequence;
								   runBlast(seq, pictureName.replaceAll(".txt", i + ".txt"));
							   }
						} catch (Exception e1) {
							e1.printStackTrace();
							}
						}
						 catch (Exception ex) {
						      ex.printStackTrace();
						   }
						}	
					}
					else if (seq2Button.isSelected() && main.byGene==true) {
						JFileChooser saveFile =new JFileChooser();
						saveFile.setDialogTitle("Specify a file root to save");
						File workingDirectory = new File(System.getProperty("user.dir"));    
						saveFile.setCurrentDirectory(workingDirectory);
						saveFile.addChoosableFileFilter(new FileNameExtensionFilter(".txt", "Text File"));
						int retrival = saveFile.showSaveDialog(null);
						if (retrival == JFileChooser.APPROVE_OPTION) {
						  try {
						   String pictureName = saveFile.getSelectedFile().getAbsolutePath();
						   if (pictureName.contains(".txt")==false) {
							   pictureName+=".txt";
						   }
						   try {
							   for (int i = startnum; i<endnum; i++) {
								   seq = main.sequencesY.get(i).sequence;
								   runBlast(seq, pictureName.replaceAll(".txt", i + ".txt"));
							   }
						} catch (Exception e1) {
							e1.printStackTrace();
							}
						}
						 catch (Exception ex) {
						      ex.printStackTrace();
						   }
						}	
					}
						
					else {
						System.out.println("you did not select a sequence");
					}
	        		}
	        	}
	        	}
	        });
	        BLASTButton.setBounds(250, 190, 150, 25);
	        BLASTButton.setToolTipText("BLASTs specified windows against NCBI nr database");
	        panel.add(BLASTButton);  
	    }
	  
		public void runGUI() {
	      frame = new JFrame("Analysis");
	      frame.setSize(500, 550);
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
	      JPanel panel = new JPanel();    
	      frame.add(panel);
	      placeComponents(panel);
	      frame.setVisible(true);
		}
		
		public File getFile() {
			JFileChooser saveFile =new JFileChooser();
			String pictureName = null;
			saveFile.setDialogTitle("Specify a file to save");
			File workingDirectory = new File(System.getProperty("user.dir"));    
			saveFile.setCurrentDirectory(workingDirectory);
			saveFile.addChoosableFileFilter(new FileNameExtensionFilter(".fasta", "FASTA File"));
			int retrival = saveFile.showSaveDialog(null);
			if (retrival == JFileChooser.APPROVE_OPTION) {
			   pictureName = saveFile.getSelectedFile().getAbsolutePath();
			   if (pictureName.contains(".fasta")==false) {
				   pictureName+=".fasta";
			   }
			}
			File fileToSave = new File(pictureName);
			return fileToSave;
		}
		
		public File getFileTxt(String input) {
			JFileChooser saveFile =new JFileChooser();
			String pictureName = null;
			saveFile.setDialogTitle(input);
			File workingDirectory = new File(System.getProperty("user.dir"));    
			saveFile.setCurrentDirectory(workingDirectory);
			saveFile.addChoosableFileFilter(new FileNameExtensionFilter(".txt", "Text File"));
			int retrival = saveFile.showSaveDialog(null);
			if (retrival == JFileChooser.APPROVE_OPTION) {
			   pictureName = saveFile.getSelectedFile().getAbsolutePath();
			   if (pictureName.contains(".txt")==false) {
				   pictureName+=".txt";
			   }
			}
			File fileToSave = new File(pictureName);
			return fileToSave;
		}
		
		public File getFileFasta(String input) {
			JFileChooser saveFile =new JFileChooser();
			String pictureName = null;
			saveFile.setDialogTitle(input);
			File workingDirectory = new File(System.getProperty("user.dir"));    
			saveFile.setCurrentDirectory(workingDirectory);
			saveFile.addChoosableFileFilter(new FileNameExtensionFilter(".fasta", "FASTA File"));
			int retrival = saveFile.showSaveDialog(null);
			if (retrival == JFileChooser.APPROVE_OPTION) {
			   pictureName = saveFile.getSelectedFile().getAbsolutePath();
			   if (pictureName.contains(".fasta")==false) {
				   pictureName+=".fasta";
			   }
			}
			File fileToSave = new File(pictureName);
			return fileToSave;
		}
		
		public static void runBlast(String seq, String filename) throws FileNotFoundException {
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
	    	//logger.info("Blasting");
	    	Blast blast = new Blast(put, get);
	    	blast.run();
	    	out.close();
		}
		
		public void outToFile(File file, String output, boolean type, int start, int end) {
			String id = null;
			if (type == true) {
				id = main.seqXname;
			}
			else {
				id = main.seqYname;
			}
			PrintStream out = null;
			try {
				out = new PrintStream(file);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			out.println(">" + id + "[" + start + " , " + end + "]");
			int i = 0;
			for (i = 0; i<=(output.length()/70)-1; i++) {
				out.println(output.substring(i*70, (i+1)*70));
		     }
		   out.println(output.substring(i*70, output.length()));
		out.close();
		}
		
		public void spitOutGenes (File file, Vector<Genes> output) {
			PrintStream out = null;
			try {
				out = new PrintStream(file);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			for (Genes a : output) {
				out.println(a.ID);
				int i = 0;
				for (i = 0; i<=(a.sequence.length()/70)-1; i++) {
					out.println(a.sequence.substring(i*70, (i+1)*70));
			     }
			   out.println(a.sequence.substring(i*70, a.sequence.length()));
			}
			out.close();
		}
		
		@Override
		public void actionPerformed(ActionEvent e) {
			//TODO write out a header at the beginning of the file (fasta format)
			int startnum = Integer.parseInt(startboxinput.getText());
			int endnum = Integer.parseInt(endboxinput.getText());
			boolean runnable = true;//for error handling
			
			//error handling of start and stop for BLAST
			if(startnum >= endnum){
				runMessage("Start window cannot be larger than end window");
				runnable = false;
			}
			if(startnum < 0){
				runMessage("that window number is not viable");
				runnable = false;
			}
			if(seq1Button.isSelected() && endnum > main.numberIterationsX){
				runMessage("These windows do not exist");
				runnable = false;
			}
			if(seq2Button.isSelected() && endnum > main.numberIterationsY){
				runMessage("These windows do not exist");
				runnable = false;
			}
			//
			
			if(runnable == true){
			if (seq1Button.isSelected() && main.byGene==false) {
				outToFile(getFile(), main.x.substring(startnum * main.slidingSize, endnum * main.slidingSize), true, startnum*main.slidingSize, endnum*main.slidingSize);
			}
			else if (seq2Button.isSelected() && main.byGene==false) {
				outToFile(getFile(), main.y.substring(startnum * main.slidingSize, endnum * main.slidingSize), false, startnum*main.slidingSize, endnum*main.slidingSize);
			}
			else if (seq1Button.isSelected() && main.byGene==true) {
				Vector<Genes> toOutput = new Vector<Genes>();
				for (int i = startnum; i<endnum; i++) {
					toOutput.addElement(main.sequencesX.get(i));
				}
				spitOutGenes(getFile(), toOutput);
			}
			else if (seq2Button.isSelected() && main.byGene==true) {
				Vector<Genes> toOutput = new Vector<Genes>();
				for (int i = startnum; i<endnum; i++) {
					toOutput.addElement(main.sequencesY.get(i));
				}
				spitOutGenes(getFile(), toOutput);
			}
			else {
				System.out.println("you did not select a sequence");
			}
			}
		}
}
