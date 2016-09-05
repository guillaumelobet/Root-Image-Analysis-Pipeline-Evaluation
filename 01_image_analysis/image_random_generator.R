# Guillaume Lobet - University of Liege
# Root systems Random Generator

# The aim of this script is to run ArchiSimple in a batch mode in order to create
# any number of root systems, with any combinaison of parameters.
# The script works as follow:
#   - Looping over the parameter set (user defined)
#   - For each combinaison, create the param.txt file used by ArchiSimple
#   - Run ArchiSimple and generate the corresponding RSML file
#   - Once all the root systems are generated, create the images using RSML_reader.jar
#   - Once all the images are generated, analyse them using RIAJ.jar

options(scipen=999) # Disable scientific notation

#-----------------------------------------------------------------------------------------
#--------------------------- GENERAL OPTIONS ---------------------------------------------
#-----------------------------------------------------------------------------------------

# Main directory, where everthing is stored
dir.base <- "/Users/guillaumelobet/Dropbox/research/projects/research/0_segment/segment_scripts/"

# Where is ArchiSimple folder
setwd(paste0(dir.base, "01_image_generator/")) 


# Simulation parameters
verbatim          <- F      # Display messages
delete_old        <- T      # Delete old simulation data
create_rsml       <- T      # Create the RSML files ArchiSimple
create_images     <- T      # Using RSML_reader.jar
analyse_images    <- T      # Using MARIA_J.jar

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


# Range of parameters
P_specie_range              <- c("dicot")   # Type of species ("monocot", "dicot")
repetitions                 <- 500                      # number of repetitions
replicates                  <- 10                      # Number of simulation by parameter sets ("synthetype")
P_maxSegments               <- 8000                  # Max number of segments allowed in the simulation (this prevent the simulation to get out of hands)


# Define parameters, either manually, or from and experimental dataset
P_duree_range               <- c(20,20)               # Total length lenght of the simulation [days] 
P_penteVitDiam_range        <- c(20, 60)              # Slope between the diameter and the root growth [-]
P_angInitMoyVertPrim_range  <- c(0.8, 1.5)            # Emission angle for the principal roots. Between 0 and PI/2 [radian]
P_intensiteTropisme_range   <- c(0.01, 0.3)           # strenght of the gravitropic response
P_propDiamRamif_range       <- c(0.3, 0.9)            # Relation between the diameter of a parent and a child root
P_distRamif_range           <- c(0.5, 5)              # Distance between two successive lateral roots [mm]
P_diamMax_range             <- c(0.2, 0.4)            # Max diameter for the primary roots [mm]
P_maxLatAge_range           <- c(5, 20)               # Maximal age growing age for the laterals [days]
P_angLat_range              <- c(0.8, 1.5)            # Emission angle for the laterals [radian]

P_nbMaxPrim_dicot           <- c(1,1)                # Number of primary axes. Put c(10,70) to have a monocot, c(1,1) to have a dicot
P_coeffCroissRad_dicot      <- c(0.1,0.5)            # Coefficient of radial growth. 0 for monocots
P_nbMaxPrim_monocot         <- c(1,20)               # Number of primary axes. Put c(10,70) to have a monocot, c(1,1) to have a dicot
P_coeffCroissRad_monocot    <- c(0,0)                # Coefficient of radial growth. 0 for monocots

P_vitEmissionPrim_range     <- c(0.1, 2)             # Speed of emission of the primary roots [root/day]
P_ageMaturitePointe_range   <- c(0.1, 3)             # Maturity age for the root tip [day]
P_tertiary                  <- 0

#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------

# Root system parameters

## General
P_penteDureeVieDiamTMD  <- 2e9
P_probaMaxArret         <- 0.5 # Probability of a root to stop growing
P_TMD                   <- 0.1
P_type                  <- 1        # Type of simulation (1 = 3D, 2 = 2D, 3 = shovelomics)
## Primary roots
P_slopePrimAngle        <- -0.3 # Slope btw the age of the and the insertion angle of the primary
P_simultEmiss           <- 0 # Emmission of seminal roots, 3 days after the start (0 = NO, 1 = YES)
P_nbSeminales           <- 5 # Number of seminal roots

## Secondary roots
P_diamMin               <-  0.02# 0.0014 # Min diameter for a root [mm]
P_coeffVarDiamRamif     <- 0.5 # Variation coefficient for the diameter fof the laterals

# Simulation parameters
P_IC_meca               <- 0   # Mecanical impedence of the soil
P_shovel                <- 70 # Depth of sampling for the shovelomics



# Variable for the computation time evolution
percent     <- 0
counter     <- 0
t           <- Sys.time()
tot_time    <- 0
mean_time   <- NA

for(P_specie in P_specie_range){

    # Set up the directories
    dir.out       <- paste0(dir.base, "00_data")
    dir.results   <- paste0(dir.out, "/results")
    dir.figs      <- paste0(dir.out, "/figures")
    dir.sp        <- paste0(dir.out, "/", P_specie)
    dir.img       <- paste0(dir.sp,"/images")
    dir.rsml      <- paste0(dir.sp,"/rsml")
    dir.shape     <- paste0(dir.sp,"/shapes")
    
    if(delete_old){
      unlink(dir.img, recursive = TRUE, force = TRUE)
      unlink(dir.rsml, recursive = TRUE, force = TRUE)
      unlink(dir.shape, recursive = TRUE, force = TRUE)
    }
    dir.create(dir.out, showWarnings = F)
    dir.create(dir.results, showWarnings = F)
    dir.create(dir.sp, showWarnings = F)
    dir.create(dir.img, showWarnings = F)
    dir.create(dir.rsml, showWarnings = F)
    dir.create(dir.shape, showWarnings = F)
    dir.create(dir.figs, showWarnings = F)
    
    
    
    
    if(P_specie == "dicot"){
      P_nbMaxPrim_range       <- P_nbMaxPrim_dicot
      P_coeffCroissRad_range  <- P_coeffCroissRad_dicot
    }else{
      P_nbMaxPrim_range       <- P_nbMaxPrim_monocot
      P_coeffCroissRad_range  <- P_coeffCroissRad_monocot      
    }
    
    

######################## Create the RSML files
    if(create_rsml){
      setwd("./archisimple/") # Set the correct repository
      
      params_range <- data.frame(penteVitDiam = numeric(), 
                                nbMaxPrim = numeric(),
                                vitEmissionPrim = numeric(), 
                                angInitMoyVertPrim = numeric(),
                                intensiteTropisme = numeric(), 
                                propDiamRamif = numeric(), 
                                distRamif = numeric(), 
                                coeffCroissRad = numeric(), 
                                diamMax = numeric(), 
                                angLat = numeric(),
                                maxLatAge = numeric(),
                                duree = numeric())
        
      for(i in 1:repetitions){ # Repetitions
        
          # Change the name pre-fix based on the simualution type
          basename <- paste0(dir.rsml,"/", P_specie)
         
          # Get the different parameters, from a random continuous distribution
          P_nbMaxPrim             <- round(runif(1, P_nbMaxPrim_range[1], P_nbMaxPrim_range[2]), 0)
          P_diamMax               <- round(runif(1, P_diamMax_range[1], P_diamMax_range[2]), 3)
          P_vitEmissionPrim       <- round(runif(1, P_vitEmissionPrim_range[1], P_vitEmissionPrim_range[2]), 3)
          P_angLat                <- round(runif(1,P_angLat_range[1],P_angLat_range[2]),3)
          P_angInitMoyVertPrim    <-round(runif(1,P_angInitMoyVertPrim_range[1],P_angInitMoyVertPrim_range[2]),3)
          P_penteVitDiam          <- round(runif(1,P_penteVitDiam_range[1],P_penteVitDiam_range[2]),3)
          P_intensiteTropisme     <- round(runif(1,P_intensiteTropisme_range[1],P_intensiteTropisme_range[2]),5)
          P_ageMaturitePointe     <- round(runif(1,P_ageMaturitePointe_range[1],P_ageMaturitePointe_range[2]),3) 
          P_propDiamRamif         <- round(runif(1,P_propDiamRamif_range[1],P_propDiamRamif_range[2]),3)
          P_distRamif             <- round(runif(1,P_distRamif_range[1],P_distRamif_range[2]),3) 
          P_maxLatAge             <- round(runif(1,P_maxLatAge_range[1],P_maxLatAge_range[2]),3)
          P_duree                 <- round(runif(1,P_duree_range[1],P_duree_range[2]),0)
          P_coeffCroissRad        <- round(runif(1,P_coeffCroissRad_range[1],P_coeffCroissRad_range[2]),3)

          for(k in 1:replicates){
            # Setup the name of the file, containing the principal info about the simulation
            name <- paste(basename, 
                          "sim",
                          i, k, 
                          sep="-")
            if(verbatim) message(name)
            
            var <- c(P_duree,
                   P_simultEmiss,
                   P_vitEmissionPrim,
                   P_nbSeminales,
                   P_nbMaxPrim,
                   P_diamMin,
                   P_diamMax,
                   P_penteVitDiam,
                   P_intensiteTropisme,
                   P_distRamif,
                   P_propDiamRamif,
                   P_coeffVarDiamRamif,
                   P_probaMaxArret,
                   P_TMD,
                   P_penteDureeVieDiamTMD,
                   P_coeffCroissRad,
                   P_angLat,
                   P_tertiary,
                   name,
                   P_type,
                   P_IC_meca,
                   P_shovel,
                   P_maxLatAge,
                   P_angInitMoyVertPrim,
                   P_slopePrimAngle,
                   P_ageMaturitePointe,
                   P_maxSegments
            )   
            cat(var, file="param.txt", sep='\n') # Create the input file for Archisimple
            t <- Sys.time()
            
            system("./ArchiSimp5Maria")  # Run Archisimple
            
            t1 <- Sys.time() - t
            tot_time <- tot_time + t1
            counter <- counter+1
            
          }

          # Counter to track the evolution of the simulations
          prog <- round((counter / (repetitions * replicates * length(P_specie_range))) * 100)
          if(prog > percent){
            message(paste0(prog, "% of root systems generated"))
            message(paste0("------ Mean simulation time = ", round(((tot_time/i)),2)))
            message(paste0("------ Total simulation time = ",round((tot_time),2)))
            message("-------------------------------------------------")
            percent <- percent + 5
          }
      }
      setwd("../") # Back to the old repo
    }
  
    
####################### Create and save the images  ########################

    if(create_images){
      message("Creating root images and ground truth data")
      system(paste0('java -Xmx6000m -jar RSML_reader.jar ',dir.rsml,' ',dir.img,' "../00_data/results/root_data_',P_specie,'-',repetitions,'.csv"'))
    }
  
######################## Analyse the images  ########################
    
    if(analyse_images){
      message("Analysing root images")
      system(paste0('java -Xmx6000m -jar RIAJ.jar ',dir.img,' ../00_data/results/root_estimators_',P_specie,'-',repetitions,'.csv 300 true true ',dir.shape))  
    }
}
#------------------------------------------------------------------------

# Print the results
diff = Sys.time() - t
print(paste((repetitions * replicates * length(P_specie_range))," root systems were generated and analysed in ",diff," seconds"), quote=F)

#------------------------------------------------------------------------


# References

# Pagès et al. (2013). 
# Calibration and evaluation of ArchiSimple, a simple model of root system architecture, 
# 290, 76–84. doi:10.1016/j.ecolmodel.2013.11.014
