# Image generation and analysis

The script in this folder were used to create the library of simulated root images and to analyse them using the RIA-J image analysis pipeline. 

## Contents



### image_random_generator.R

R script used to generated the image library and analyse the different images. The script will generate a user-defined number of root images, together with their ground-truth and their image descriptors. 

For each root system, the script perfumes the following steps:

1.  generates a random set of input parameters for ArchiSimple (within  limits set by the user). 
2. run ArchiSimple and creates an RSML file with all the root system data. 
3. run the RSML_reader plugin, which create the corresponding root image and the stores the ground-truth data. 
4. run the RIA_J plugin and extract the different descriptors from the root image.

### RSML_reader.jar

Executable for the RSML_reader plugin. Can be run (on Unix) using:

	java -jar RSML_reader.jar [rsml_files_directory] [image_files_directory] [ground_truth_data.csv] [degradetion_level]

The `degradation_level` indicates if you want to degrade the resulting image (0 = no degradation)

### RSML_reader_source

Source files of the RSML_reader plugin for ImageJ.

### RIAJ.jar

Executable for the RSML_reader plugin. Can be run (on Unix) using:

	java -jar RIAJ.jar [image_files_directory] [descriptors.csv] [dpi] [morphometric_analysis] [outline_analysis] [shape_directory]
	

### RIAJ_source

Source files of the RIA-J plugin for ImageJ


### ArchiSimple

For copyright reasons, ArchiSimple is not archived here. If interested, please contact Loïc Pagès [loic.pages@inra.fr]()
