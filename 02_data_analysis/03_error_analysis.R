

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

color.plot         <- T                             # Make the plots in color
species            <- c("tap-rooted","fibrous")             # the type of species to analyse
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
library(reshape2)
library(lattice)
library(latticeExtra)
library(ggrepel)
library(plyr)


#---------------------------------------------------------------
#---------------------------------------------------------------
# HOME MADE FUNCTIONS
#---------------------------------------------------------------
#---------------------------------------------------------------

leg.col <- c("black", "grey50")
leg.col.3 <- c("black", "grey50", "grey30")
ramp.corr <- colorRamp(c("black", "white", "black"))
ramp <- colorRamp(c("white", "black"))
if(color.plot){
  leg.col <- c("#b2df8a", "#1f78b4")
  leg.col.3 <- c("#bae4bc", "#7bccc4", "#2b8cbe")
  ramp <- colorRamp(c("white", "#91bfdb", "#fc8d59"))
  ramp.corr <- colorRamp(c("#91bfdb", "white", "#fc8d59"))
}



get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}




#---------------------------------------------------------------
#---------------------------------------------------------------
# LOAD THE PROCESSED DATASETS
#---------------------------------------------------------------
#---------------------------------------------------------------


descr <- read.csv(paste0(n,"-descriptors.csv"))[,-1]
data <- read.csv(paste0(n,"-data.csv"))[,-1]


# Re-compute the index vector for further uses
de.names <- colnames(descr)
da.names <- colnames(data)
names.shape <- de.names[grepl("PC", de.names)]
names.tree <- da.names[grepl("tree", da.names)]

ind.type.de <- match(c("image","specie", "group", "id"), de.names)
ind.shape <- match(c(names.shape, 
                     "width50", "count50", "width_depth_ratio", "com_y", "com_x", "convexhull", "exploration"), 
                   de.names)
ind.type.da <- match(c("image","specie", "group", "id"), da.names)
ind.tree <- match(names.tree, da.names)

de.list <- de.names[-ind.type.de]
da.list <- da.names[-ind.type.da]


# Change all the dicot to "tap-rooted" and the monocot to "fibrous"

data$specie <- as.character(data$specie)
data$specie[data$specie == "monocot"] <- "fibrous"
data$specie[data$specie == "dicot"] <- "tap-rooted"
data$specie <- factor(data$specie)


descr$specie <- as.character(descr$specie)
descr$specie[descr$specie == "monocot"] <- "fibrous"
descr$specie[descr$specie == "dicot"] <- "tap-rooted"
descr$specie <- factor(descr$specie)


#---------------------------------------------------------------
#---------------------------------------------------------------
# OVERLAP EVOLUTION
#---------------------------------------------------------------
#---------------------------------------------------------------



temp <- data[,c("tot_root_length", "overlap", "specie")]
temp$tot_root_length <- round(temp$tot_root_length/100)*100
temp <- ddply(temp, .(tot_root_length, specie), summarize, overlap=mean(overlap))

ggplot(temp, aes(tot_root_length, overlap, colour=specie)) + 
  geom_point(alpha=0.8, size=2) + 
  theme_classic() +
  facet_wrap(~specie)+
  theme(text=element_text(size=20), 
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  xlab("Total root system length [mm]") + 
  ylab("overlap index [%]") +
  scale_color_manual(values=leg.col)+
  ggsave(file=paste0(base.dir,"/00_data/figures/overlap_evolution.pdf"), width = 10, height = 4)+
  ggsave(file=paste0(base.dir,"/00_data/figures/overlap_evolution.png"), width = 10, height = 4)


#---------------------------------------------------------------
#---------------------------------------------------------------
# SPECIFIC ERROR ANALYSIS
#---------------------------------------------------------------
#---------------------------------------------------------------

# We can have a look at the error made for the different parameters 
# and see how this error is influenced by the root system properties. 
# For instance the idea is to look at the ` length` vs `tot_root_length` 
# and compute the error as a function of the `overlap index`

pairs <- data.frame(data=c("tot_root_length", "n_laterals", "depth"), 
                    descr=c("length", "tip_count", "depth"), 
                    stringsAsFactors = F)

names <- data.frame(data=c("Total root length [mm]", "Number of lateral roots [-]", "Root system depth [mm]"), 
                    descr=c("Estimated length [mm]", "Estimated tip count [-]", "Estimated depth [mm]"), 
                    stringsAsFactors = F)

for(p in 1:nrow(pairs)){
  temp <- data.frame(x=data[[pairs$data[p]]], y=descr[[pairs$descr[p]]], z=data[["overlap"]], 
                     specie = data$specie, deg = factor(descr$degradation))
  
  # Plot the error
  g1 <- ggplot(temp, aes(x, y, colour=factor(deg))) + 
    geom_point(alpha=0.2, size=2) + 
    geom_abline(intercept = 0, slope = 1, alpha=0.5) + 
    stat_smooth(method="lm")+
    facet_wrap(~specie) + 
    theme_classic() +
    theme(legend.position = "none", text=element_text(size=20), 
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
    xlab(names$data[p]) + 
    ylab(names$descr[p]) + 
    scale_color_manual(values=leg.col.3) + 
    #coord_fixed() + 
    ylim(c(0, max(temp$y))) + 
    xlim(c(0, max(temp$x))) + 
    theme(panel.border = element_blank(),
          legend.key = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          panel.grid = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "transparent",colour = NA))
  
  temp$z <- round(temp$z*100)/100
  temp$error <- abs((temp$x - temp$y) / temp$x)
  
  # compute the Root Square Mean Error (RSME)
  rsme <- NULL
  for(sp in species){
    for(d in unique(temp$deg)){
      zz <- unique(temp$z[temp$specie == sp & temp$deg == d])
      temp1 <- data.frame(group = zz, value=numeric(length(zz)), specie=sp, degradation = as.character(d), stringsAsFactors = F)
      for(i in 1:length(zz)){
        temp1$value[i] <- mean(temp$error[temp$z == zz[i] & temp$specie == sp & temp$deg == d])
      }
      rsme <- rbind(rsme, temp1)
    }
  }
  rsme <- rsme[!is.infinite(rsme$value),]
  
  rsme$degradation[rsme$degradation == 0] <- "0 - null"
  rsme$degradation[rsme$degradation == 1] <- "1 - medium"
  rsme$degradation[rsme$degradation == 3] <- "2 - high"
  
  pd <- position_dodge(width=0.1)
  
  g2 <- ggplot(rsme, aes(x=group, y=value, colour=degradation)) + 
    geom_point(size=2, position = pd, alpha=0.5) + 
    stat_smooth(se = F, method="loess") + 
    theme_classic() +
    facet_wrap(~specie) + 
    theme(text=element_text(size=20), 
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
    scale_y_continuous(limits = c(0, max(1, max(rsme$value)))) + 
    scale_color_manual(values=leg.col.3) +
    scale_x_continuous(limits = c(0, 1)) + 
    xlab("Overlap index [-]") + 
    ylab("RELATIVE ERROR [-]") 
  
  
  lg <- get_legend(g2)
  
  g2 <- g2 + theme(legend.position = "none")
  
  g <- arrangeGrob(g1, g2, lg, widths=c(0.55, 0.35, 0.1), ncol=3)
  ggsave(paste0(base.dir,"/00_data/figures/error-",pairs$data[p],".pdf"), g, width=18, height=5) + 
    ggsave(paste0(base.dir,"/00_data/figures/error-",pairs$data[p],".png"), g, width=18, height=5)
  
}





#---------------------------------------------------------------
#---------------------------------------------------------------
# META ERROR ANALYSIS
#---------------------------------------------------------------
#---------------------------------------------------------------

# We can have a look at the error made for the different parameters 
# and see how this error is influenced by the root system properties. 
# For instance the idea is to look at the ` length` vs `tot_root_length` 
# and compute the error as a function of the `overlap index`


# Find the best descriptor for each ground-truth
to_est <-c("tot_root_length",         # Vector of parameters to estimate with the machine learning
                "width",
                "depth",
                "mean_prim_length",
                "n_laterals",
                "tot_lat_length",
                "mean_lat_angle")

vars <- colnames(descr)[-c(1,16, 17, 18, 19)]


couples <- NULL
# Look at the quality of the estimations
for(te in to_est){
  
  # Get the best single variable (regression)
  max <- 0
  keep <- ""
  for(v in vars){
    fit <- lm(data[[te]][data$degradation == 0] ~ descr[[v]][descr$degradation == 0])
    r2 <- summary(fit)$r.squared
    if(r2 > max){
      max <- r2
      keep <- v
    }
  }
  
  couples <- rbind(couples, data.frame(data=te, descr=keep, value = max, type="r2", deg = "0"))
  
  # Get the best single variable (MRE)
  min <- 1e9
  keep <- ""
  for(v in vars){
    estimation <- descr[[v]][descr$degradation == 0]
    truth <- data[[te]][data$degradation == 0]
    
    rel_diff <- abs((truth - estimation) / truth)
    rel_diff[is.infinite(rel_diff)] <- 0
    diff <- mean(rel_diff, na.rm = T)
    
    if(diff < min){
      min <- diff
      keep <- v
    }
  }
  couples <- rbind(couples, data.frame(data=te, descr=keep, value = min, type="MRE", deg="0"))
}


# Complete the dataset with the regressions for the degraded images
for(d in unique(descr$degradation)){
  if(d != "0"){
    for(te in to_est){
      message(paste0(d, " / ", te))
      # For the linear regression
      v <- as.character(couples$descr[couples$data == te & couples$type == "r2" & couples$deg == "0"][1])
      fit <- lm(data[[te]][data$degradation == "0"] ~ descr[[v]][descr$degradation == d])
      r2 <- summary(fit)$r.squared
      couples <- rbind(couples, data.frame(data=te, descr=v, value = r2, type="r2", deg = factor(d)))
      
      # For the relative error
      v <- as.character(couples$descr[couples$data == te & couples$type == "MRE" & couples$deg == "0"][1])
      estimation <- descr[[v]][descr$degradation == d]
      truth <- data[[te]][data$degradation == "0"]
      rel_diff <- abs((truth - estimation) / truth)
      rel_diff[is.infinite(rel_diff)] <- 0
      diff <- mean(rel_diff, na.rm = T)
      couples <- rbind(couples, data.frame(data=te, descr=v, value = diff, type="MRE", deg = factor(d)))
    }
  }
}

pos <- position_dodge(0.5)
ggplot(data=couples, aes(x=data, y=value, colour=deg)) +
  geom_point(size=3, position=pos)+
  facet_wrap(~type, scale="free") + 
  theme_classic() +
  ylim(range(0,2)) + 
  #scale_color_manual(values=leg.col)+
  xlab("")+
  ylab("Accuracy evaluation") + 
  theme(axis.text.x=element_text(angle=45, hjust = 1), legend.position="none",
        axis.line.x = element_line(colour = 'black', size=0.3, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'))
+ 
  ggsave(paste0("../figures/rf_accuracy.pdf"), height=4, width=6)+
  ggsave(paste0("../figures/rf_accuracy.png"), height=4, width=6)






