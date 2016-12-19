

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

leg.col <- c("black", "grey50")
ramp.corr.g <- colorRamp(c("black", "white", "black"))
ramp.g <- colorRamp(c("white", "black"))
if(color.plot){
  leg.col <- c("#b2df8a", "#1f78b4")
  ramp <- colorRamp(c("white", "#91bfdb", "#fc8d59"))
  ramp.corr <- colorRamp(c("#91bfdb", "white", "#fc8d59"))
}



#---------------------------------------------------------------
#---------------------------------------------------------------
# LOAD THE PROCESSED DATASETS
#---------------------------------------------------------------
#---------------------------------------------------------------


descr <- read.csv(paste0(n,"-descriptors.csv"))[,-1]
data <- read.csv(paste0(n,"-data.csv"))[,-1]


# for this part keep only the non degraded images
descr <- descr[descr$degradation == 0,]
data <- data[data$degradation == 0,]

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
# GENERAL ANALYSIS OF THE DATA STRUCTURE
#---------------------------------------------------------------
#---------------------------------------------------------------


# Distribution of the ground-truth data. 
names.tree <- da.names[grepl("tree", da.names)]
tree.ind <- match(c(names.tree, "direction", "degradation", "overlap"), colnames(data))
temp <- data[, -tree.ind]

for(sp in species){
  temp1 <- data[data$specie == sp,]
  print(sp)
  for(i in 1:ncol(temp1)){
    if(is.numeric(temp1[,i])){
      message(paste0( " ", colnames(temp1)[i], " ", round(min(temp1[,i]), 2), " ", round(max(temp1[,i]), 2)))
    }
  }
}

rs <- melt(temp, id.vars = c("image", "specie", "group", "id"))

ggplot(data = rs, aes(x=value, fill=specie)) + 
  geom_density(stat="bin", colour="white", alpha=0.8) + 
  xlab(n) + 
  theme_classic() + 
  facet_wrap(~variable, scales = "free") + 
  theme(axis.text.x=element_text(angle=45, hjust = 1)) + 
  scale_fill_manual(values=leg.col) + 
  ggsave(file=paste0(base.dir,"/00_data/figures/hist-data.pdf"), width = 10, height = 6) +
  ggsave(file=paste0(base.dir,"/00_data/figures/hist-data.png"), width = 10, height = 6)


# Distribution of the descriptors. 

rs <- melt(descr, id.vars = c("image", "specie", "group", "id"))

ggplot(data = rs, aes(x=value, fill=specie)) + 
  geom_density(stat="bin", colour="white", alpha=0.8) + 
  xlab(n) + 
  theme_classic() + 
  facet_wrap(~variable, scales = "free") + 
  theme(axis.text.x=element_text(angle=45, hjust = 1)) + 
  scale_fill_manual(values=leg.col) + 
  ggsave(file=paste0(base.dir,"/00_data/figures/hist-descr.pdf"), width = 10, height = 8) +
  ggsave(file=paste0(base.dir,"/00_data/figures/hist-descr.png"), width = 10, height = 8)



# Analysis of the ground-truth data using a Principal Component Analasis. 
# This is used to look at the global structure of the dataset and see if we have differences
# between the different species (monocot / dicot)

colnames(data[,-c(ind.type.da, tree.ind, 5, 20)])

pca <- prcomp(data[,-c(ind.type.da, tree.ind, 5, 20)], retx = T, scale=T)  # Make the PCA
rs <- data.frame(pca$x)[,c(1,2)]    # Store the loadings 
rs$specie <- data$specie
rs$image <- data$image
vars <- apply(pca$x, 2, var)  
props <- round((vars / sum(vars) * 100), 1)
xl <- paste0("\nPrincipal Component 1 (",props[1],"%)")
yl <-paste0("Principal Component 2 (",props[2],"%)\n")


ggplot(data =rs, aes(x=PC1, y=PC2, colour=specie)) + 
  geom_point(alpha=0.8, size=2) + 
  stat_ellipse(level = 0.95, size=1.2, alpha=0.8) + 
  theme_classic() + 
  theme(text = element_text(size=25), 
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
  xlab(xl) +
  ylab(yl) + 
  #scale_shape_manual(values = c(19,1)) + 
  scale_colour_manual(values=leg.col) + 
  ggsave(file=paste0(base.dir,"/00_data/figures/pca-global.pdf"), width = 8, height = 7) +
  ggsave(file=paste0(base.dir,"/00_data/figures/pca-global.png"), width = 8, height = 7)




z2 <- data.frame(var_names = rownames(pca$rotation), pca$rotation[, 1:2])
z2$var_names <- gsub("_", " ", z2$var_names)

ggplot(data=z2, aes(0, 0, xend=PC1, yend=PC2)) + 
  geom_segment(col="grey", size=2, arrow = arrow(length = unit(0.5,"cm")), alpha=0.9) +
  geom_text_repel(data=z2, aes(PC1, PC2, label=var_names), col="black", size=9) +
  geom_point(aes(x=0, y=0), size=5, colour="grey") +
  #scale_y_continuous(limits = c(-1, 0.3)) +
  theme_classic() +
  xlab(xl) + ylab(yl) + 
  theme(text = element_text(size=25), axis.line = element_line(size = 1.5), 
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
  ggsave(file=paste0(base.dir,"/00_data/figures/pca-loadings.pdf"), width = 10, height = 10)


# PERFOM A MANOVA ANALYSIS TO SEE WETHER THE OVERALL DIFFERENCE BETWEEN GROUPS IS SIGNIFICANT

fit <- manova(cbind(tot_root_length, width, depth,
                    direction, n_primary, tot_prim_length,   
                    mean_prim_length, mean_prim_diameter, mean_lat_density,
                    n_laterals, tot_lat_length, mean_lat_length,   
                    mean_lat_diameter, mean_lat_angle) ~ specie,  data) # Stele was removed
summary(fit, test="Roy")
summary(fit, test="Pillai")
summary(fit, test="Wilks")
summary(fit, test="Hotelling-Lawley")
summary.aov(fit)



#---------------------------------------------------------------
#---------------------------------------------------------------
# ANALYSIS OF GROUND-TRUTH DATA
#---------------------------------------------------------------
#---------------------------------------------------------------

### Relationship between (within) the ground-truth data

temp.list <- da.list[c(1:(length(da.list)-6))]
temp.list <- temp.list[-4]

r2.mat <- matrix(0, nrow = length(temp.list), ncol=length(temp.list))
colnames(r2.mat) <- temp.list
rownames(r2.mat) <- temp.list

data.m <- data[data$specie == "fibrous",-ind.tree]
data.d <- data[descr$specie == "tap-rooted",-ind.tree]

# Building of a correlation matrix
for(da1 in temp.list){
  for(da2 in temp.list){
    if(da1 != da2){
      fit.m <- lm(data.m[[da1]] ~ data.m[[da2]])
      r2.m <- round(summary(fit.m)$r.squared, 4)
      
      fit.d <- lm(data.d[[da1]] ~ data.d[[da2]])
      r2.d <- round(summary(fit.d)$r.squared, 4)    
      
      r2.mat[da1, da2] <- r2.m
      r2.mat[da2, da1] <- r2.d
    }
  }
}


# Plotting the correlation as a heatmap

plo <- levelplot(r2.mat, 
                 col.regions=rgb(ramp(seq(0, 1, length = 1000)), max = 255), 
                 at = seq(0, 1, length.out=100),
                 scales=list(x=list(rot=45, cex=1.5), y=list(cex=1.5), xlab=list(cex=.5)), 
                 xlab = "MONOCOTS", ylab = "DICOTS")
pdf(paste0(base.dir,"/00_data/figures/heatmap-internal.pdf"), width = 8, height=8)
plo
dev.off()
plo




#---------------------------------------------------------------

### Relationship between ground-truth data and descriptors.

# get the names of the descriptors
de.list.2 <- c("length","area","diam_mean","tip_count",
               "width", "depth", "width_depth_ratio","com_x", "com_y", 
               "convexhull", "exploration")

# Get the names of the ground-truth data
da.list.2 <- c("tot_root_length","n_primary", "tot_prim_length", "mean_prim_length", "mean_prim_diameter", 
               "mean_lat_density", "n_laterals", "tot_lat_length", "mean_lat_length", "mean_lat_diameter",
               "width", "depth")

### R-SQUARED
i <- 0
for(sp in species){
  i <- i+1
  r2.mat <- matrix(0, nrow = length(da.list.2), ncol=length(de.list.2))
  colnames(r2.mat) <- de.list.2
  rownames(r2.mat) <- da.list.2
  
  # Create the correlation matrix
  for(da in da.list.2){
    for(de in de.list.2){
      fit <- lm(data[[da]][data$specie == sp] ~ descr[[de]][descr$specie == sp])
      r2 <- round(summary(fit)$r.squared, 4)
      
      r2.mat[da, de] <- r2
      temp <- data.frame(x=data[[da]][data$specie == sp], y=descr[[de]][descr$specie == sp])
    }
  }
  
  # Plot the data as a heatmap
  
  pl <- levelplot(r2.mat, 
                  col.regions=rgb(ramp.g(seq(0, 1, length = 1000)), max = 255), 
                  at = seq(0, 1, length.out=100),
                  scales=list(x=list(rot=45, cex=1.5), y=list(cex=1.5), xlab=list(cex=.5)), 
                  xlab = "Ground-truth data", ylab = "Image descriptors", main="")
  if(i == 1) pl1 <- pl
  if(i == 2) pl2 <- pl
  message(paste0(length(r2.mat[r2.mat > 0.8]), "/", length(r2.mat)))
}
pdf(paste0(base.dir,"/00_data/figures/heatmap-all.pdf"), width = 13, height=13)
c(pl1, pl2, layout = c(1, 2), merge.legends = F)
dev.off()
png(paste0(base.dir,"/00_data/figures/heatmap-all.png"))
c(pl1, pl2, layout = c(1, 2), merge.legends = F)
dev.off()



### PEARSONS CORRELATIONS
i <- 0
for(sp in species){
  cor <- cor(data[data$specie == sp, da.list.2], descr[descr$specie == sp,de.list.2])
  # Plot the data as a heatmap
  pl <- levelplot(cor, 
                  col.regions=rgb(ramp.corr(seq(0, 1, length = 1000)), max = 255), 
                  at = seq(-1, 1, length.out=100),
                  scales=list(x=list(rot=45, cex=1.5), y=list(cex=1.5), xlab=list(cex=.5)), 
                  xlab = "Ground-truth data", ylab = "Image descriptors", main="")
  i <- i+1
  if(i == 1) pl1 <- pl
  if(i == 2) pl2 <- pl
}
pdf(paste0(base.dir,"/00_data/figures/heatmap-all-corr.pdf"), width = 10, height=10)
c(pl1, pl2, layout = c(1, 2), merge.legends = F)
dev.off()
