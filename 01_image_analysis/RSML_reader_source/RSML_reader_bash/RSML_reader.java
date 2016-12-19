
/**
 * @author Xavier Draye - Universit� catholique de Louvain
 * @author Guillaume Lobet - Universit� de Li�ge
 *   
 * Main class for the RSML improter
 */

import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;

import ij.*;

public class RSML_reader {

	private static RSML_reader instance;
	public static String path, output, data, imageSize;
	public static boolean printData = false;
	public static File[] imgs;
	
	/**
	 * Constructor
	 */
	public RSML_reader() {
     
      // Retrieve all the rsml files
      File f = new File(path);
      System.out.println(f.getAbsolutePath());
      File[] rsml = f.listFiles(new FilenameFilter() {
    	  public boolean accept(File directory, String fileName) {
    		  return fileName.endsWith(".rsml");
    	  }
      });
      
      // Retrieve all the images files
      File f1 = new File(output);
      imgs = f1.listFiles(new FilenameFilter() {
    	  public boolean accept(File directory, String fileName) {
    		  return fileName.endsWith(".jpg");
    	  }
      });      
      
      // Open the different RSML files, retrieve their data and get their size.
      RootModel model;
      PrintWriter pw = null;
      
      if(printData){
    	  try{pw = new PrintWriter(new FileWriter(data));}
    	  catch(IOException e){System.out.println("Could not save file "+data);}
    	  pw.println("image,tot_root_length,width,depth,direction,n_primary,tot_prim_length,mean_prim_length,mean_prim_diameter,mean_lat_density,n_laterals," +
    		  " tot_lat_length,mean_lat_length,mean_lat_diameter,mean_lat_angle,tree_size,tree_magnitude,tree_altitude,tree_index");
      }
      int percent = 0;
      float progression = 0;
      for(int i = 0; i < rsml.length; i++){
    	  model = new RootModel(rsml[i].getAbsolutePath());
    	  if(model.getNRoot() > 0){
    		  // Save the image
    		  if(!imageExist(rsml[i].getName())){
    			  ImagePlus ip;
    			  if(imageSize.equals("1")){
    				  ip = new ImagePlus(rsml[i].getName(),model.createImage(false, 0, true, false, false));  
    			  }else{
    				  ip = new ImagePlus(rsml[i].getName(),model.createFixedImage(false, 0, true, false, false));  
    			  }
        		  IJ.save(ip, output+System.getProperty("file.separator")+rsml[i].getName()+".jpg");	   
    		  }
    		  progression = (i/rsml.length)*100;
    		  if(progression > percent){    			
    			  IJ.log(percent+" % of the rsml files converted. "+(rsml.length-i)+" files remaining.");
    			  percent = percent + 10;
    		  }
    		      		  
    		  // Save the data
    		  if(printData) model.sendImageData(pw, rsml[i].getName());
    	  } 	  	
      }
      if(printData) pw.flush();
	}

	
	/**
	 * Constructor
	 */
	public static void RSMLToPhenoPackets() {
     
		
      // Retrieve all the rsml files
      File f = new File(path);
      System.out.println(f.getAbsolutePath());
      File[] rsml = f.listFiles(new FilenameFilter() {
    	  public boolean accept(File directory, String fileName) {
    		  return fileName.endsWith(".rsml");
    	  }
      });
      
      // Retrieve all the images files
      File f1 = new File(output);
      imgs = f1.listFiles(new FilenameFilter() {
    	  public boolean accept(File directory, String fileName) {
    		  return fileName.endsWith(".jpg");
    	  }
      });      
      
      // Open the different RSML files, retriev their data and get their size.
      RootModel model;
      PrintWriter pw = null;
      

	  try{pw = new PrintWriter(new FileWriter(data));}
	  catch(IOException e){System.out.println("Could not save file "+data);}
	  pw.println("# Test for RSML to PhenoPacket conversion");
	  pw.println("schema: rsml_example");
	  pw.println("comment: example of export from RSML files to PhenoPackets"); 

      // Encode the different entities (the different files)
	  pw.println("entities:");
      for(int i = 0; i < rsml.length; i++){
    	  model = new RootModel(rsml[i].getAbsolutePath());
    	  if(model.getNRoot() > 0){
        	  pw.println(" - id: "+rsml[i].getName());
        	  pw.println("   type: maize");

    	  } 	  	
      }
      
      // Encode the phenotypes
	  pw.println("phenotype_profile:");
      for(int i = 0; i < rsml.length; i++){
    	  model = new RootModel(rsml[i].getAbsolutePath());
    	  if(model.getNRoot() > 0){    		  
    		  // Save the data
    		  System.out.println("TEST");

    	      pw.println(" - entity: "+rsml[i].getName());
    	      pw.println("   evidence: ");
    	      pw.println("    type: image");
    	      pw.println("    source:");
    	      pw.println("     id: "+rsml[i].getName());
    	      
    		  model.sendImageDataToPhenoPackets(pw, rsml[i].getName());
    	  } 	  	
      }
      pw.flush();
	}	
	
	
	
	
	public static boolean imageExist(String s){
		
		for(int i = 0; i < imgs.length; i++){
			if(imgs[i].getName().substring(0, imgs[i].getName().length()-9).equals(s.substring(0, s.length()-5))) return true; 
		}
		return false;
	}
	
   /**
    * Get instance
    * @return
    */
   public static RSML_reader getInstance() {return instance; }

   
   /**
    * Main class
    * @param args
    */
   @SuppressWarnings("unused")
   public static void main(String args[]) {
	   if(args.length > 0){
		   path = args[0];
		   if(args.length > 1) output = args[1];
		   else output = path;
		   if(args.length > 2){
			   printData = true;
			   data = args[2];
		   }
		   if(args[3].equals("1")){
			   imageSize = args[4];
			   RSML_reader ie = new RSML_reader();
		   }
		   else{
			   RSMLToPhenoPackets();
		   }
	   }
	   else System.out.println("No path specified");
	   System.exit(0);
   }

}




