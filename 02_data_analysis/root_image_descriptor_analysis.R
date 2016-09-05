
  
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

compute.raw.data   <- F                            # If the shape analysis has already been done, put TRUE to reuse the generated datafile
color.plot         <- F                             # Make the plots in color
species            <- c("dicot","monocot")             # the type of species to analyse
base.dir           <- "~/Dropbox/research/projects/research/0_segment/segment_scripts"
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
library(reshape2)
library(lattice)
library(latticeExtra)
library(MASS)




#---------------------------------------------------------------
#---------------------------------------------------------------
# HOME MADE FUNCTIONS
#---------------------------------------------------------------
#---------------------------------------------------------------

leg.col <- c("black", "grey50")
ramp.corr <- colorRamp(c("black", "white", "black"))
ramp <- colorRamp(c("white", "black"))
if(color.plot){
  leg.col <- c("red", "blue")
  ramp <- colorRamp(c("white", "lightblue", "red"))
  ramp.corr <- colorRamp(c("blue", "white", "red"))
}

getGroup <- function(name){
  n1 <- gsub("/Users/guillaumelobet/Desktop/Work/segment/outputs/monocot/rsml/", "", name)
  strsplit(n1, "-")[[1]][3]
}

getID <- function(name){
  n1 <- gsub("/Users/guillaumelobet/Desktop/Work/segment/outputs/monocot/rsml/", "", name)
  strsplit(n1, "-")[[1]][4]
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
      temp <- read.csv(paste0("root_data_",sp,"-",n,".csv"))
      temp$specie <- sp
      data <- rbind(data, temp)
    }
    
    # Create a group (synthetype) and an id
    data$group <- sapply(data$image, getGroup)
    data$id <- sapply(data$image, getID)
    
    # Remove the intermediate datasets
    remove(tmp)
    
    
    # Loading the Descriptors
    da.list <- colnames(data[,2:(ncol(data)-2)])
    
    descr <- NULL
    for(sp in species){
      temp <- read.csv(paste0("root_estimators_",sp,"-",n,".csv"))
      temp$specie <- sp
      descr <- rbind(descr, temp)
    }
    
    # Create a group (synthetype) and an id
    descr$group <- sapply(descr$image, getGroup)
    descr$id <- sapply(descr$image, getID)
    
    # Remove the intermediate datasets
    remove(temp)
    
    
    
    
    #---------------------------------------------------------------
    #---------------------------------------------------------------
    # SHAPE DATA PREPARATION
    #---------------------------------------------------------------
    #---------------------------------------------------------------
    
    # We need to generate the PC's for the pseudo landmark analysis and we add them to the `descr` table. 
    # We perfom the PCA within each species as we decided to split them.
    
    procShape <- NULL
    for(sp in species){
      # Load the shape analysis
      shape <- readland.tps(paste0("root_estimators_",sp,"-",n,"-shape.tps"), specID = "ID")
      
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
    
    
    # Then we need to generate the PC's for the EFD analysis using `Momocs`. Again, we add them to the `descr` table. 
    
    efd <- NULL
    for(sp in species){
      
      # Get the file list
      lf <- list.files(paste0(base.dir,"/00_data/",sp,"/shapes"), full.names=TRUE)
      
      # import the images
      roots <- import_jpg(lf)
      
      # Get the outline of the shapes
      roots.o <- Out(roots)
      
      # Calibrate the number of harmonics needed
      g <- calibrate_harmonicpower(roots.o, nb.h = 50, plot=F)
      rs <- melt(g$q[,c(1:20)])
      ggplot(rs, aes(factor(Var2), value)) + geom_boxplot() + theme_classic() + 
      xlab("Harmonics number [-]") + ylab("Cumulative sum harmonic power [%]") +
      ggsave(paste0(base.dir,"/00_data/figures/harmonic-power-",sp,".pdf"), width=7, height=5)
      
      roots.f <- efourier(roots.o, nb.h = 20)
      
      # Then, do the PCA
      roots.p <- PCA(roots.f)
      
      
      pdf(paste0(base.dir,"/00_data/figures/efd-pca-",sp,".pdf"))
      plot(roots.p, pos.shp = "full_axes", points=T, labelspoints = F)
      dev.off()
      
      temp <- data.frame(image = row.names(roots.p$x), 
      EFD.PC1 = roots.p$x[,1],
      EFD.PC2 = roots.p$x[,2],
      EFD.PC3 = roots.p$x[,3])
      efd <- rbind(efd, temp)
    }
    
    efd$image <- gsub("_shape", "", efd$image)
    descr <- merge(descr, efd, by = "image")
    
    
    
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
    
    ind.type.de <- match(c("image", "specie", "group", "id"), de.names)
    ind.type.da <- match(c("image","specie", "group", "id"), da.names)
    
    names.cum <- de.names[grepl("cumul", de.names)]
    names.coord <- de.names[grepl("coord", de.names)]
    names.diff <- de.names[grepl("diff", de.names)]
    names.crossh <- de.names[grepl("cross_ho", de.names)]
    names.ell <- de.names[grepl("ellips", de.names)]
    names.rect <- de.names[grepl("rect", de.names)]
    ind.remove <- match(c(names.crossh, names.diff, names.coord, names.cum, names.rect, names.ell), de.names)
    
    ind.size <- match(c("depth","width"), de.names)
    
    names.shape <- de.names[grepl("PC", de.names)]
    ind.shape <- match(c(names.shape, 
    "width50", "count50", "width_depth_ratio", "com_y", "com_x", "convexhull", "exploration"), 
    de.names)
    
    ind.length.de <- match(c("length", "area", "width", "depth", "tip_count"), de.names)
    names.length <- da.names[grepl("length", da.names)]
    ind.length.da <- match(c(names.length, "width", "depth", "n_laterals", "n_primary", "mean_lat_density"), da.names)
    
    
    ind.diam.de <- match(c("diam_max", "diam_mode", "diam_mean"), de.names)
    names.diam <- da.names[grepl("diam", da.names)]
    ind.diam.da <- match(names.diam, da.names)
    
    names.tree <- da.names[grepl("tree", da.names)]
    ind.tree <- match(names.tree, da.names)
    
    
    # Keep a full dataset and a reduced one
    descr.full <- descr
    shape <- descr[,c(ind.type.de, ind.shape)]
    length.de <- descr[,c(ind.type.de, ind.length.de)]
    diam.de <- descr[,c(ind.type.de, ind.diam.de)]
    descr <- descr[,-c(ind.remove)]
    
    data.full <- data
    tree <- data[, c(ind.type.da, ind.tree)]
    length.da <- data[, c(ind.type.da, ind.length.da)]
    diam.da <- data[, c(ind.type.da, ind.diam.da)]
    
    
    
    # Then, finaly, we re-order the data, as the `merge` step might have messed them up. 
    # We also save them, so we do not need to load them every time we run the script. 
    
    # Re order the data to be sure to have the same tables (merge might have messed it up)
    data <- data[order(data$image),]
    descr <- descr[order(descr$image),]
    tree <- tree[order(tree$image),]
    shape <- shape[order(shape$image),]
    length.da <- length.da[order(length.da$image),]
    length.de <- length.de[order(length.de$image),]
    diam.da <- diam.da[order(diam.da$image),]
    diam.de <- diam.de[order(diam.de$image),]
    
    write.csv(descr, paste0(n,"-descriptors.csv"))
    write.csv(descr.full, paste0(n,"-descriptors_full.csv"))
    write.csv(data, paste0(n,"-data.csv"))
    write.csv(tree, paste0(n,"-tree-.csv"))
    write.csv(data.full, paste0(n,"-data_full.csv"))
    write.csv(shape, paste0(n,"-shape.csv"))
    write.csv(length.da, paste0(n,"-length-data.csv"))
    write.csv(length.de, paste0(n,"-length-descr.csv"))
    write.csv(diam.da, paste0(n,"-diam-data.csv"))
    write.csv(diam.de, paste0(n,"-diam-descr.csv"))
    
    
    # So, now, we have the full data set and we can start to play with it
    
    }



#---------------------------------------------------------------
#---------------------------------------------------------------
# LOAD THE PROCESSED DATASETS
#---------------------------------------------------------------
#---------------------------------------------------------------


descr <- read.csv(paste0(n,"-descriptors.csv"))[,-1]
descr.full <- read.csv(paste0(n,"-descriptors_full.csv"))[,-1]
data <- read.csv(paste0(n,"-data.csv"))[,-1]
tree <- read.csv(paste0(n,"-tree-.csv"))[,-1]
data.full <- read.csv(paste0(n,"-data_full.csv"))[,-1]
shape <- read.csv(paste0(n,"-shape.csv"))[,-1]
length.da <- read.csv(paste0(n,"-length-data.csv"))[,-1]
length.de <- read.csv(paste0(n,"-length-descr.csv"))[,-1]
diam.da <- read.csv(paste0(n,"-diam-data.csv"))[,-1]
diam.de <- read.csv(paste0(n,"-diam-descr.csv"))[,-1]


# Re-compute the index vector for further uses
de.names <- colnames(descr)
da.names <- colnames(data)
tree.names <- colnames(tree)
shape.names <- colnames(shape)
length.da.names <- colnames(length.da)
length.de.names <- colnames(length.de)
diam.da.names <- colnames(diam.da)
diam.de.names <- colnames(diam.de)
names.shape <- de.names[grepl("PC", de.names)]
names.tree <- da.names[grepl("tree", da.names)]

ind.type.de <- match(c("image","specie", "group", "id"), de.names)
ind.shape <- match(c(names.shape, 
"width50", "count50", "width_depth_ratio", "com_y", "com_x", "convexhull", "exploration"), 
de.names)
ind.type.da <- match(c("image","specie", "group", "id"), da.names)
ind.type <- match(c("image","specie", "group", "id"), tree.names)
ind.tree <- match(names.tree, da.names)
ind.size.da <- match(c("width", "depth"), length.da.names)

de.list <- de.names[-ind.type.de]
da.list <- da.names[-ind.type.da]
tree.list <- tree.names[-ind.type]
shape.list <- shape.names[-ind.type]
length.da.list <- length.da.names[-c(ind.type, ind.size.da)]
length.de.list <- length.de.names[-ind.type]
diam.da.list <- diam.da.names[-ind.type]
diam.de.list <- diam.de.names[-ind.type]







#---------------------------------------------------------------
#---------------------------------------------------------------
# GENERAL ANALYSIS OF THE DATA STRUCTURE
#---------------------------------------------------------------
#---------------------------------------------------------------


# Distribution of the ground-truth data. 
names.tree <- da.names[grepl("tree", da.names)]
tree.ind <- match(c(names.tree, "direction"), colnames(data))
temp <- data[, -tree.ind]

for(sp in species){
  temp <- data[data$specie == sp,]
  print(sp)
  for(i in 1:ncol(temp)){
    if(is.numeric(temp[,i])){
      message(paste0( "& \\verb|", colnames(temp)[i], "| & ", round(min(temp[,i]), 2), " & ", round(max(temp[,i]), 2), "\\\\"))
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
  ggsave(file=paste0(base.dir,"/00_data/figures/hist-data.pdf"), width = 10, height = 6)


# Distribution of the descriptors. 

rs <- melt(descr, id.vars = c("image", "specie", "group", "id"))

ggplot(data = rs, aes(x=value, fill=specie)) + 
  geom_density(stat="bin", colour="white", alpha=0.8) + 
  xlab(n) + 
  theme_classic() + 
  facet_wrap(~variable, scales = "free") + 
  theme(axis.text.x=element_text(angle=45, hjust = 1)) + 
  scale_fill_manual(values=leg.col) + 
  ggsave(file=paste0(base.dir,"/00_data/figures/hist-descr.pdf"), width = 10, height = 8)


# Analysis of the ground-truth data using a Principal Component Analasis. 
# This is used to look at the global structure of the dataset and see if we have differences
# between the different species (monocot / dicot)

pca <- prcomp(data[,-c(ind.type.da, ind.tree)], retx = T, scale=T)  # Make the PCA
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
  ggsave(file=paste0(base.dir,"/00_data/figures/pca-global.pdf"), width = 8, height = 7)



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




#---------------------------------------------------------------
#---------------------------------------------------------------
# ANALYSIS OF GROUND-TRUTH DATA
#---------------------------------------------------------------
#---------------------------------------------------------------

### Relationship between (within) the ground-truth data

temp.list <- da.list[c(1:(length(da.list)-4))]

r2.mat <- matrix(0, nrow = length(temp.list), ncol=length(temp.list))
colnames(r2.mat) <- temp.list
rownames(r2.mat) <- temp.list

data.m <- data[data$specie == "monocot",-ind.tree]
data.d <- data[descr$specie == "dicot",-ind.tree]

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

shape.ind <- match(shape.names, de.list)
de.list.2 <- de.list[-shape.ind[!is.na(shape.ind)]][-c(9,10)]

tree.ind <- match(c(tree.names, "direction"), da.list)
da.list.2 <- da.list[-tree.ind[!is.na(tree.ind)]]

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
      
      mean(abs(descr[[de]][descr$specie == sp] - data[[da]][data$specie == sp]) / data[[da]][data$specie == sp])
      
      r2.mat[da, de] <- r2
      temp <- data.frame(x=data[[da]][data$specie == sp], y=descr[[de]][descr$specie == sp])
    }
  }

  # Plot the data as a heatmap
  
  pl <- levelplot(r2.mat, 
                  col.regions=rgb(ramp(seq(0, 1, length = 1000)), max = 255), 
                  at = seq(0, 1, length.out=100),
                  scales=list(x=list(rot=45, cex=1.5), y=list(cex=1.5), xlab=list(cex=.5)), 
                  xlab = "Ground-truth data", ylab = "Image descriptors", main="")
  if(i == 1) pl1 <- pl
  if(i == 2) pl2 <- pl
}
pdf(paste0(base.dir,"/00_data/figures/heatmap-all.pdf"), width = 10, height=10)
  c(pl1, pl2, layout = c(1, 2), merge.legends = F)
dev.off()
# print(comb_levObj)


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




#---------------------------------------------------------------
#---------------------------------------------------------------
# ERROR ANALYSIS
#---------------------------------------------------------------
#---------------------------------------------------------------

# We can have a look at the error made for the different parameters 
# and see how this error is influenced by the root system properties. 
# For instance the idea is to look at the ` length` vs `tot_root_length` 
# and compute the error as a functino of the `n_laterals`

pairs <- data.frame(data=c("tot_root_length", "n_laterals", "depth"), 
                    descr=c("length", "tip_count", "depth"), 
                    stringsAsFactors = F)

names <- data.frame(data=c("Total root length [cm]", "Number of lateral roots [-]", "Root system depth [cm]"), 
                    descr=c("Estimated length [cm]", "Estimated tip count [-]", "Estimated depth [cm]"), 
                    stringsAsFactors = F)

for(p in 1:nrow(pairs)){
  temp <- data.frame(x=data[[pairs$data[p]]], y=descr[[pairs$descr[p]]], z=data[["tot_root_length"]], specie = data$specie)

  # Plot the error
  g1 <- ggplot(temp, aes(x, y, colour=specie)) + 
    geom_point(alpha=0.8, size=2) + 
    geom_abline(intercept = 0, slope = 1, alpha=0.5) + 
    facet_wrap(~specie) + 
    theme_classic() +
    theme(legend.position = "none", text=element_text(size=20), 
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
    xlab(names$data[p]) + 
    ylab(names$descr[p]) + 
    scale_color_manual(values=leg.col) + 
    coord_fixed() + 
    ylim(c(0, max(temp$x, temp$y))) + 
    xlim(c(0, max(temp$x, temp$y)))

  if(p < 3) temp$z <- round(temp$z/50)*50
  else temp$z <- round(temp$z/50)*50
  temp$error <- (temp$x - temp$y) / temp$x
  
  # compute the Root Square Mean Error (RSME)
  rsme <- NULL
  for(sp in species){
    zz <- unique(temp$z[temp$specie == sp])
    temp1 <- data.frame(group = zz, value=numeric(length(zz)), specie=sp)
    for(i in 1:length(zz)){
      temp1$value[i] <- sqrt(sum(temp$error[temp$z == zz[i] & temp$specie == sp]^2) / 
      length(temp$error[temp$z == zz[i] & temp$specie == sp]))
    }
    rsme <- rbind(rsme, temp1)
  }
  rsme <- rsme[!is.infinite(rsme$value),]
  
  pd <- position_dodge(width=10)
  
  g2 <- ggplot(rsme, aes(group, value, colour=specie)) + 
    #geom_line(position = pd) + 
    geom_point(size=2, position = pd) + 
    stat_smooth(se = F, method="loess") + 
    theme_classic() + 
    theme(text=element_text(size=20), 
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
    scale_y_continuous(limits = c(0, 1.5)) + 
    xlab("Total root length [cm]") + 
    ylab("RRSME [-]") + 
    scale_color_manual(values=leg.col)
  
  
  lg <- get_legend(g2)
  
  g2 <- g2 + theme(legend.position = "none")
  
  g <- arrangeGrob(g1, g2, lg, widths=c(0.55, 0.35, 0.1), ncol=3)
  ggsave(paste0(base.dir,"/00_data/figures/error-",pairs$data[p],".pdf"), g, width=14, height=4)

}






#---------------------------------------------------------------
#---------------------------------------------------------------
# SHAPE ANALYSIS
#---------------------------------------------------------------
#---------------------------------------------------------------


# First, we will look at the different shaope varaibles and 
# see if they are related to any of the classical variables

da.list.2 <- c(length.da.list)

### R-SQUARED
i <- 0
for(sp in species){
  r2.mat <- matrix(0, nrow = length(da.list.2), ncol=length(shape.list))
  colnames(r2.mat) <- shape.list
  rownames(r2.mat) <- da.list.2
  
  for(da in da.list.2){
    for(sh in shape.list){
      fit <- lm(data[[da]][data$specie == sp] ~ shape[[sh]][shape$specie == sp])
      r2 <- round(summary(fit)$r.squared, 4)
      r2.mat[da, sh] <- r2
    }
  }
  

  pl <- levelplot(r2.mat, 
                  col.regions=rgb(ramp(seq(0, 1, length = 1000)), max = 255), 
                  at = seq(0, 1, length.out=100),
                  scales=list(x=list(rot=45, cex=1.5), y=list(cex=1.5), xlab=list(cex=.5)), 
                  xlab = "Ground-truth data", ylab = "Shape descriptors", main=toupper(sp))
  
  i <- i+1 
  if(i == 1) pl1 <- pl
  if(i == 2) pl2 <- pl
}
pdf(paste0(base.dir,"/00_data/figures/heatmap-shape.pdf"), width = 12, height=10)
  c(pl1, pl2, layout = c(2, 1), merge.legends = F)
dev.off()


### PEARSONS CORRELATIONS
i <- 0
for(sp in species){
  cor <- cor(data[data$specie == sp, da.list.2], shape[shape$specie == sp,shape.list])
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
pdf(paste0(base.dir,"/00_data/figures/heatmap-shape-corr.pdf"), width = 10, height=10)
  c(pl1, pl2, layout = c(2, 1), merge.legends = F)
dev.off()


# Then we look if there are relationship between the shape variables themselves

r2.mat <- matrix(0, nrow = length(shape.list), ncol=length(shape.list))
colnames(r2.mat) <- shape.list
rownames(r2.mat) <- shape.list

corr.mat <- matrix(0, nrow = length(shape.list), ncol=length(shape.list))
colnames(corr.mat) <- shape.list
rownames(corr.mat) <- shape.list

shape.m <- shape[shape$specie == "monocot",]
shape.d <- shape[shape$specie == "dicot",]
for(da1 in shape.list){
  for(da2 in shape.list){
    if(da1 != da2){
      fit.m <- lm(shape.m[[da1]] ~ shape.m[[da2]])
      r2.m <- round(summary(fit.m)$r.squared, 4)

      fit.d <- lm(shape.d[[da1]] ~ shape.d[[da2]])
      r2.d <- round(summary(fit.d)$r.squared, 4)    

      r2.mat[da1, da2] <- r2.m
      r2.mat[da2, da1] <- r2.d
      
      corr.mat[da1, da2] <- cor(shape.m[[da1]] , shape.m[[da2]])#r2.m
      corr.mat[da2, da1] <- cor(shape.d[[da1]] , shape.d[[da2]])#r2.d      
    }
  }
}



plo <- levelplot(r2.mat, 
                 col.regions=rgb(ramp(seq(0, 1, length = 1000)), max = 255), 
                 at = seq(0, 1, length.out=100),
                 scales=list(x=list(rot=45, cex=1.5), y=list(cex=1.5), xlab=list(cex=.5)), 
                 xlab = "MONOCOTS", ylab = "DICOTS")
pdf(paste0(base.dir,"/00_data/figures/heatmap-shape-internal.pdf"), width = 8, height=8)
plo
dev.off()
plo


plo <- levelplot(corr.mat, 
                 col.regions=rgb(ramp.corr(seq(0, 1, length = 1000)), max = 255), 
                 at = seq(-1, 1, length.out=100),
                 scales=list(x=list(rot=45, cex=1.5), y=list(cex=1.5), xlab=list(cex=.5)), 
                 xlab = "MONOCOTS", ylab = "DICOTS")
pdf(paste0(base.dir,"/00_data/figures/heatmap-shape-internal-corr.pdf"), width = 8, height=8)
plo
dev.off()
plo


#---------------------------------------------------------------
#---------------------------------------------------------------
# PREDICTION ANALYSES
#---------------------------------------------------------------
#---------------------------------------------------------------

# An important feature of image descriptors is to be able to discriminate the different genotypes between them. 
# Knowing how many and which descriptors are needed can indeed speed up the overall analysis by removing 
# unnecessary metrics (that could take time to acquire). 

# Here, we used a Linear Discriminant Analysis (LDA) in order to estimate the prediction power of each variables. 
# Prediction power was estimated by (1) creating a model with half of the indivudas of each synthetype and 
# (2) using the model to predict the synthetype of the remaining plants. 

# The prediction accuracy was computed multiple times, by iterativelly adding new descriptors to the model. At each
# iteration, the added descriptor was the one that would increase the most the global accuracy of the model.


descr.small <- descr[descr$length > 200 & descr$length < 300,]

ind <- ind.type.de[c(1,2,3,5)]
to_analyse <- de.names[-ind.type.de]
ind <- match(to_analyse, colnames(descr))
ind.g <- match("group", colnames(descr))
lda.data <- NULL

done <- c()

for(sp in species){
  ind <- match(to_analyse, colnames(descr))
  done <- c()
  acc.all <- c()
  for(k in 1:length(ind)){
    max <- 0
    keep <- 0
    for(i in 1:length(ind)){
      temp <- c(done, ind[i])
      fit <- lda(group ~ ., data=descr.small[descr.small$specie == sp & descr.small$id <= 5, c(temp, ind.g)])
      fit.values <- predict(fit, descr.small[descr.small$specie == sp & descr.small$id > 5, c(temp, ind.g)])
      ct <- table(fit.values$class, descr.small$group[descr.small$specie == sp & descr.small$id > 5])      
      acc <- round(sum(diag(prop.table(ct)))*100, 2)
      if(acc > max){
        max <- acc
        keep <- i
      }
    }
    done <- c(done, ind[keep])
    acc.all <- c(acc.all, max)
    ind <- ind[-keep]    
    message(done)
    message(ind)
    message(max)
    message("--------------")
  
  }
  names <- done
  for(i in 1:length(done)) names[i] <- colnames(descr)[done[i]]
  temp <- data.frame(names, acc.all, id=c(1:length(names)), species = sp)
  lda.data <- rbind(lda.data, temp)
}


temp <- ind.shape
for(i in 1:length(ind.shape)) temp[i] <- colnames(descr)[ind.shape[i]]

lda.data$type <- "morphology"
lda.data$type[lda.data$names %in% temp] <- "shape"
lda.data$diff <- 0
lda.data$diff[lda.data$species == "dicot"] <- c(lda.data$acc.all[lda.data$species == "dicot"][1],
diff(lda.data$acc.all[lda.data$species == "dicot"]))
lda.data$diff[lda.data$species == "monocot"] <- c(lda.data$acc.all[lda.data$species == "monocot"][1],
diff(lda.data$acc.all[lda.data$species == "monocot"]))

plot1 <- ggplot(lda.data[lda.data$species == "dicot",], aes(reorder(names, id), acc.all, fill=type)) + 
  geom_line(aes(group="none")) +
  geom_hline(yintercept = 90, lty=2) +
  geom_point(size=5, pch=21) +
  ylim(range(lda.data$acc.all)) + 
  # scale_x_discrete(labels=lda.data$names) +
  scale_fill_manual(values=c("white", "black")) + 
  theme_bw() +
  # scale_y_continuous(breaks = seq(60, 100, 20)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position="none", text = element_text(size=20), 
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  ggtitle("DICOTS") + 
  xlab("Variable name") + 
  ylab("Prediction accuracy [%]\n")

plot2 <- ggplot(lda.data[lda.data$species == "monocot",], aes(reorder(names, id), acc.all, shape=type, fill=type)) + 
  geom_line(aes(group="none")) +
  geom_hline(yintercept = 90, lty=2) +
  geom_point(size=5, pch=21) +
  ylim(range(lda.data$acc.all)) + 
  # scale_x_discrete(labels=lda.data$names) +
  scale_fill_manual(values=c("white", "black")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none", text = element_text(size=20)) +
  ggtitle("MONOCOTS") + 
  xlab("Variable name") + 
  ylab("  ")

grid.arrange(plot1, plot2, ncol=2)
g <- arrangeGrob(plot1, plot2, ncol=2)
ggsave(paste0(base.dir,"/00_data/figures/precision-analysis.pdf"), g, width=16, height=8)


#---------------




data.small <- data#[data$depth > 10 & data$depth < 20,]

ind <- ind.type.de[c(1,2,3,5)]
to_analyse <- da.names[-ind.type.da]
ind <- match(to_analyse, colnames(data))
ind.g <- match("group", colnames(data))
lda.data <- NULL

done <- c()

for(sp in species){
  ind <- match(to_analyse, colnames(data))
  done <- c()
  acc.all <- c()
  for(k in 1:length(ind)){
    if(length(ind) > 0){
      max <- 0
      keep <- 0
      for(i in 1:length(ind)){
        if(sd(data.small[data.small$specie == sp & data.small$group == 1, ind[i]]) > 0){
          temp <- c(done, ind[i])
          fit <- lda(group ~ ., data=data.small[data.small$specie == sp & data.small$id <= 5, c(temp, ind.g)])
          fit.values <- predict(fit, data.small[data.small$specie == sp & data.small$id > 5, c(temp, ind.g)])
          ct <- table(fit.values$class, data.small$group[data.small$specie == sp & data.small$id > 5])      
          acc <- round(sum(diag(prop.table(ct)))*100, 2)
          if(acc > max){
            max <- acc
            keep <- i
          }
        }
      }
      
      done <- c(done, ind[keep])
      if(max > 0) acc.all <- c(acc.all, max)
      ind <- ind[-keep]    
      message(done)
      message(ind)
      message(max)
      message("--------------")
    }
  }
  names <- done
  for(i in 1:length(done)) names[i] <- colnames(data)[done[i]]
  temp <- data.frame(names, acc.all, id=c(1:length(names)), species = sp)
  lda.data <- rbind(lda.data, temp)
}


temp <- ind.shape
for(i in 1:length(ind.shape)) temp[i] <- colnames(data)[ind.shape[i]]

lda.data$type <- "morphology"
lda.data$type[lda.data$names %in% temp] <- "shape"
lda.data$diff <- 0
lda.data$diff[lda.data$species == "dicot"] <- c(lda.data$acc.all[lda.data$species == "dicot"][1],
                                                diff(lda.data$acc.all[lda.data$species == "dicot"]))
lda.data$diff[lda.data$species == "monocot"] <- c(lda.data$acc.all[lda.data$species == "monocot"][1],
                                                  diff(lda.data$acc.all[lda.data$species == "monocot"]))

plot1 <- ggplot(lda.data[lda.data$species == "dicot",], aes(reorder(names, id), acc.all, fill=type)) + 
  geom_line(aes(group="none")) +
  geom_hline(yintercept = 90, lty=2) +
  geom_point(size=5, pch=21) +
  ylim(range(lda.data$acc.all)) + 
  # scale_x_discrete(labels=lda.data$names) +
  scale_fill_manual(values=c("white", "black")) + 
  theme_bw() +
  # scale_y_continuous(breaks = seq(60, 100, 20)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position="none", text = element_text(size=20), 
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + 
  ggtitle("DICOTS") + 
  xlab("Variable name") + 
  ylab("Prediction accuracy [%]\n")

plot2 <- ggplot(lda.data[lda.data$species == "monocot",], aes(reorder(names, id), acc.all, shape=type, fill=type)) + 
  geom_line(aes(group="none")) +
  geom_hline(yintercept = 90, lty=2) +
  geom_point(size=5, pch=21) +
  ylim(range(lda.data$acc.all)) + 
  # scale_x_discrete(labels=lda.data$names) +
  scale_fill_manual(values=c("white", "black")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none", text = element_text(size=20)) +
  ggtitle("MONOCOTS") + 
  xlab("Variable name") + 
  ylab("  ")

grid.arrange(plot1, plot2, ncol=2)
g <- arrangeGrob(plot1, plot2, ncol=2)
ggsave(paste0(base.dir,"/00_data/figures/precision-analysis.pdf"), g, width=16, height=8)


