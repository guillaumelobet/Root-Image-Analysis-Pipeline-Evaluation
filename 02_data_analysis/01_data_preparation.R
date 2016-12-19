

# The aim of this document is to investigate the link between different root image descriptors 
# and between the root system ground truths. To do so, we generate a wide array of root systems 
# using the model ArchiSimple (Pagès et al. 2014). The model generation is done in a seprate file (image_random_generator.R). 
#Dataset were generated for monocots and dicots. All simulations are generic and were not tuned to 
# match any specific plant species. The main differences between monocots and dicots are the number of 
# principal axes and the presence of a secondary growth.

# Note about the nomenclature used here: 
# `data` will refer to the ground truth values when 
# `descriptors` refere to thedescriptors coming from the image analysis. 


#---------------------------------------------------------------
#---------------------------------------------------------------
# GENERAL PARAMETERS FOR THE ANALYSIS
#---------------------------------------------------------------
#---------------------------------------------------------------

compute.raw.data   <- T                            # If the shape analysis has already been done, put FALSE to reuse the generated datafile
species            <- c("tap-rooted", "fibrous")             # the type of species to analyse
base.dir           <- "/Users/g.lobet/OneDrive - UCL/03_research/0-segment/segment_scripts/"
n                  <- 500                              # This is the number of simulated root systems

setwd(paste0(base.dir,"/00_data/results"))



#---------------------------------------------------------------
#---------------------------------------------------------------
# LOADING LIBRARIES
#---------------------------------------------------------------
#---------------------------------------------------------------

library(ggplot2)
library(gridExtra)
library(shapes)
library(Momocs)
library(geomorph)
library(reshape2)
library(lattice)
library(latticeExtra)
library(MASS)
library(ggrepel)


#---------------------------------------------------------------
#---------------------------------------------------------------
# HOME MADE FUNCTIONS
#---------------------------------------------------------------
#---------------------------------------------------------------


getGroup <- function(name){
  n1 <- gsub("/Users/guillaumelobet/Desktop/Work/segment/outputs/monocot/rsml/", "", name)
  strsplit(n1, "-")[[1]][3]
}

getID <- function(name){
  n1 <- gsub("/Users/guillaumelobet/Desktop/Work/segment/outputs/monocot/rsml/", "", name)
  strsplit(n1, "-")[[1]][4]
}

getDegradation <- function(name){
  n1 <- gsub("/Users/guillaumelobet/Desktop/Work/segment/outputs/monocot/rsml/", "", name)
  if(!grepl("deg",n1)) return(0)
  gsub(".jpg", "", gsub("deg", "", strsplit(n1, "-")[[1]][6]))
}

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}




if(compute.raw.data){
  
  #---------------------------------------------------------------
  #---------------------------------------------------------------
  # RAW DATA LOADING
  #---------------------------------------------------------------
  #---------------------------------------------------------------
  
  # First, let's prepare the ground-truth data. We just need to add a `specie` column. 
  # We also separate the tree variables from the dataset (not sure how we will use them). 
  
  # Loading the data
  
  data <- NULL
  for(sp in species){
    if(sp == "fibrous") sp1 <- "monocot"
    if(sp == "tap-rooted") sp1 <- "dicot"
    temp <- read.csv(paste0("root_data_",sp1,"-",n,".csv"))
    temp$specie <- sp
    data <- rbind(data, temp)
  }
  
  # Create a group (synthetype) and an id
  data$group <- sapply(data$image, getGroup)
  data$id <- sapply(data$image, getID)
  data$degradation <- "0"
  
  # Remove the intermediate datasets
  remove(temp)
  da.list <- colnames(data[,2:(ncol(data)-2)])
  
  
  # Loading the Descriptors
  
  descr <- NULL
  for(sp in species){
    if(sp == "fibrous") sp1 <- "monocot"
    if(sp == "tap-rooted") sp1 <- "dicot"
    
    temp <- read.csv(paste0("root_estimators_",sp1,"-",n,".csv"))
    temp$specie <- sp
    descr <- rbind(descr, temp)
  }
  
  # Create a group (synthetype) and an id
  descr$group <- sapply(descr$image, getGroup)
  descr$id <- sapply(descr$image, getID)
  descr$degradation <- sapply(descr$image, getDegradation)
  
  
  # For the "Data" dataframe, duplicate rows, to have the nuber as the "Descr"
  
  dataInit <- data
  for(i in unique(descr$degradation)){
    if(i != "0"){
      temp <- dataInit
      temp$degradation <- i
      data <- rbind(data, temp)
    }
  }
  
  remove(temp, dataInit)
  
  
  data <- data[order(data$image, data$degradation),]
  descr <- descr[order(descr$image, descr$degradation),]
  
  
  #---------------------------------------------------------------
  #---------------------------------------------------------------
  # SHAPE DATA PREPARATION
  #---------------------------------------------------------------
  #---------------------------------------------------------------
  
  # We need to generate the PC's for the pseudo landmark analysis and we add them to the `descr` table. 
  # We perfom the PCA within each species as we decided to split them.
  
  procShape <- NULL
  for(sp in species){
    
    if(sp == "fibrous") sp1 <- "monocot"
    if(sp == "tap-rooted") sp1 <- "dicot"
    
    # Load the shape analysis
    shape <- readland.tps(paste0("root_estimators_",sp1,"-",n,"-shape.tps"), specID = "ID")
    
    # Make the procust analysis
    proc <- procGPA(shape)
    
    procShape <- rbind(procShape, proc$scores)
    
    # PLOT THE SHAPE ANALYSIS
    
    # get the mean shape
    shape <- data.frame(x = proc$mshape[,1], y = -proc$mshape[,2])
    
    # Get the deviations
    rotPC1 <-  proc$pcar[, 1]
    rotPC2 <-  proc$pcar[, 2]
    rotPC3 <-  proc$pcar[, 3]
    sdPC1 <- proc$pcasd[1]
    sdPC2 <- proc$pcasd[2]
    sdPC3 <- proc$pcasd[3]
    
    # Base plot
    plot <- ggplot(data = shape, aes(x=x, y=y)) + coord_fixed()
    for(i in 1:10){#proc$n){
      temp <- data.frame(x=proc$rotated[,1,i], y=-proc$rotated[,2,i])
      plot <- plot + geom_polygon(data=temp, aes(x=x, y=y), fill="transparent",colour="#00000010") 
    }
    plot <- plot + geom_polygon(fill="transparent", size=1, colour="#00000090")  + theme_bw()
    
    devs <- c(-3, 3)
    col <- c("#FF605D", "#00B94E")
    
    ### PC1
    plot1 <- plot + ggtitle(paste("PC1 [",round(proc$percent[1],2),"%]"))
    for(i in 1:length(devs)){
      fig <- proc$mshape + devs[i] * rotPC1 * sdPC1
      k <- length(proc$mshape[,1])
      temp <- data.frame(x= fig[1:k], y=-fig[(k + 1):(2 * k)])
      plot1 <- plot1 + geom_polygon(data=temp, aes(x=x, y=y), fill="transparent",colour=col[i], size=1) 
    }
    
    ### PC2
    plot2 <- plot + ggtitle(paste("PC2 [",round(proc$percent[2],2),"%]"))
    for(i in 1:length(devs)){
      fig <- proc$mshape + devs[i] * rotPC2 * sdPC2
      k <- length(proc$mshape)/2
      temp <- data.frame(x= fig[1:k], y=-fig[(k + 1):(2 * k)])
      plot2 <- plot2 + geom_polygon(data=temp, aes(x=x, y=y), fill="transparent",colour=col[i], size=1) 
    }
    
    
    ### PC3
    plot3 <- plot + ggtitle(paste("PC3 [",round(proc$percent[3],2),"%]"))
    for(i in 1:length(devs)){
      fig <- proc$mshape + devs[i] * rotPC3 * sdPC3
      k <- length(proc$mshape)/2
      temp <- data.frame(x= fig[1:k], y=-fig[(k + 1):(2 * k)])
      plot3 <- plot3 + geom_polygon(data=temp, aes(x=x, y=y), fill="transparent",colour=col[i], size=1) 
    }
    
    g <- arrangeGrob(plot1, plot2, plot3,  ncol=3)
    ggsave(paste0(base.dir,"/00_data/figures/PL-shape-",sp,".pdf"), g, width=10, height=5)  
    
  }
  # Get the values from the PCA
  descr$PL.PC1 <- procShape[,1]
  descr$PL.PC2 <- procShape[,2]
  descr$PL.PC3 <- procShape[,3]
  
  remove(shap, shape, procShape)
  
  #---------------------------------------------------------------
  
  
  # An other shape variable we can include is the relative cumulative width at the 50% depth mark.
  
  temp <- melt(descr, id.vars = c("image"))
  temp <- temp[grepl("cumul", temp$variable),]
  temp$variable <- factor(temp$variable)
  temp$value <- as.numeric(temp$value)
  
  for(im in unique(temp$image)){ 
    temp$value[temp$image == im] <- temp$value[temp$image == im] /
      max(temp$value[temp$image == im])
  }
  temp$variable <- as.numeric(temp$variable) / max(as.numeric(temp$variable))
  
  temp <- temp[temp$variable == 0.5, c(1,3)]
  colnames(temp) <- c("image", "width50")
  
  descr <- merge(descr, temp, by="image")
  
  
  #---------------------------------------------------------------
  
  
  # We can then make the same for the horizontal pixel count (density). 
  
  temp <- melt(descr, id.vars = c("image"))
  temp <- temp[grepl("cross_hori", temp$variable),]
  temp <- temp[grepl("mean", temp$variable),]
  temp$variable <- factor(temp$variable)
  temp$value <- as.numeric(temp$value)
  
  for(im in unique(temp$image)) {
    temp$cum[temp$image == im] <- cumsum(temp$value[temp$image == im])
    temp$rel[temp$image == im] <- temp$cum[temp$image == im] / max(temp$cum[temp$image == im])
  }
  temp$variable <- as.numeric(temp$variable) / max(as.numeric(temp$variable))
  
  temp <- temp[temp$variable == 0.5, c(1,3)]
  colnames(temp) <- c("image", "count50")
  descr <- merge(descr, temp, by="image")
  
  remove(temp)
  
  
  #---------------------------------------------------------------
  
  descr$exploration <- descr$area / descr$convexhull
  
  #---------------------------------------------------------------
  #---------------------------------------------------------------
  # CLEANING UP THE DIFERENT DATATABLE AND SAVING THEM FOR FURTHER USE
  #---------------------------------------------------------------
  #---------------------------------------------------------------
  
  
  # Now, we have a very large number of columns. We might want to get rid of some. 
  # So we get the indexes of the unwanted columns and use those to create a smaller and  more comprehensive dataset 
  
  de.names <- colnames(descr)
  da.names <- colnames(data)
  
  names.cum <- de.names[grepl("cumul", de.names)]
  names.coord <- de.names[grepl("coord", de.names)]
  names.diff <- de.names[grepl("diff", de.names)]
  names.crossh <- de.names[grepl("cross_ho", de.names)]
  names.ell <- de.names[grepl("ellips", de.names)]
  names.rect <- de.names[grepl("rect", de.names)]
  ind.remove <- match(c(names.crossh, names.diff, names.coord, names.cum, names.rect, names.ell), de.names)
  
  # ind.size <- match(c("depth","width"), de.names)
  # 
  # names.shape <- de.names[grepl("PC", de.names)]
  # ind.shape <- match(c(names.shape, 
  #                      "width50", "count50", "width_depth_ratio", "com_y", "com_x", "convexhull", "exploration"), 
  #                    de.names)
  # 
  # ind.length.de <- match(c("length", "area", "width", "depth", "tip_count"), de.names)
  # names.length <- da.names[grepl("length", da.names)]
  # ind.length.da <- match(c(names.length, "width", "depth", "n_laterals", "n_primary", "mean_lat_density"), da.names)
  # 
  # 
  # ind.diam.de <- match(c("diam_max", "diam_mode", "diam_mean"), de.names)
  # names.diam <- da.names[grepl("diam", da.names)]
  # ind.diam.da <- match(names.diam, da.names)
  # 
  # names.tree <- da.names[grepl("tree", da.names)]
  # ind.tree <- match(names.tree, da.names)
  
  
  descr <- descr[,-c(ind.remove)]


  
  # Then, finaly, we re-order the data, as the `merge` step might have messed them up. 
  # We also save them, so we do not need to load them every time we run the script. 
  
  # Re order the data to be sure to have the same tables (merge might have messed it up)
  data <- data[order(data$image),]
  descr <- descr[order(descr$image),]
  
  write.csv(descr, paste0(n,"-descriptors.csv"))
  write.csv(data, paste0(n,"-data.csv"))
  
  
  # So, now, we have the full data set and we can start to play with it
  
}


