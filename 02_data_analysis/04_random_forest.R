

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
library(reshape2)
library(lattice)
library(latticeExtra)
library(ggrepel)
library(randomForest)
source("../../02_data_analysis/maria_learning.r")  

#---------------------------------------------------------------
#---------------------------------------------------------------
# HOME MADE FUNCTIONS
#---------------------------------------------------------------
#---------------------------------------------------------------

leg.col <- c("black", "grey50")
ramp.corr <- colorRamp(c("black", "white", "black"))
ramp <- colorRamp(c("white", "black"))
leg.col.3 <- c("black", "grey50", "grey30")
leg.col.4 <- c("black", "grey50", "grey30", "grey10")

if(color.plot){
  leg.col <- c("#b2df8a", "#1f78b4")
  ramp <- colorRamp(c("white", "#91bfdb", "#fc8d59"))
  ramp.corr <- colorRamp(c("#91bfdb", "white", "#fc8d59"))
  leg.col.3 <- c("#bae4bc", "#7bccc4", "#2b8cbe")
  leg.col.4 <- c("#bae4bc", "#7bccc4", "#2b8cbe", "grey50")
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
ind.type <- match(c("image","specie", "group", "id"), tree.names)
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
# RANDOM FOREST EVALUATION
#---------------------------------------------------------------
#---------------------------------------------------------------




# Random Forest parameters

vec.models <- c(5)                 # Vector with the number of models to try
vec.trees <- c(5)                  # Vector with the number of tree to try in each model

to_est_base <-c("tot_root_length",         # Vector of parameters to estimate with the machine learning
                "true_width","true_depth",
                "n_primary","tot_prim_length","mean_prim_length",
                "mean_lat_density","n_laterals","tot_lat_length",
                "mean_lat_length","mean_lat_angle"
)


# Make the datatsets the same size
parameters <- data
descriptors <- descr

parameters$image <- substr(parameters$image, 0, nchar(as.character(parameters$image))-5)
descriptors$image <- substr(descriptors$image, 0, nchar(as.character(descriptors$image))-9)
parameters <-parameters[parameters$image %in% descriptors$image,]

dlength <- ncol(descriptors)
vars <- colnames(descriptors)[-c(1,16, 17, 18)]

colnames(parameters)[colnames(parameters) == "width"] <- "true_width"
colnames(parameters)[colnames(parameters) == "depth"] <- "true_depth"


descriptors$image <- gsub("deg1", "", descriptors$image)
descriptors$image <- gsub("deg3", "", descriptors$image)


rs <- merge(descriptors, parameters, by = c("image", "degradation"))
#rs <- rs[sample(1:nrow(rs), repetitions),] # If the dataset is bigger than wanted, downsize it

cols <- colnames(rs)
for(c in cols) rs <- rs[!is.nan(rs[[c]]),] # Remove the row containing  NA values
for(c in cols) rs <- rs[!is.na(rs[[c]]),] # Remove the row containing  NA values


#-------------------------------------------------------------------
# RUN THE MODEL
#-------------------------------------------------------------------

# We divide the dataset in training and test. Training is used to 
# the random forest models. Test is used to evaluate how good it
# perfoms.

metaPrec <- NULL
alltest = NULL
allresults = NULL
degs <- unique(rs$degradation)

for(sp in species){
  for(d in degs){  
    temp <- rs[rs$specie.x == sp & rs$degradation == d,]
    
    descr_ind <- c(2:(ncol(descriptors))) # Indexes of the descriptors
    to_est <- to_est_base
    if(sp == "tap-rooted") to_est <- to_est_base[-4]
    to_est_ind <- match(to_est, colnames(temp))
    descr_ind <- c(descr_ind, to_est_ind)  
    
    train.id <- sample(1:nrow(temp), 3*(nrow(temp)/4))
    test.id <- c(1:nrow(temp))[-train.id]
    
    train <- temp[train.id,]
    test <- temp[test.id,]
    
    for(i in descr_ind){
       message(i)
       if(is.numeric(train[,i])) print(range(train[,i]))
     }
    
    models <- GenerateModels(fname = NULL,
                             mat.data = train,
                             vec.models = vec.models,
                             vec.trees = vec.trees,
                             vec.f = to_est_ind,
                             vec.p = descr_ind)
    
    vec.weights <- rep(1, length(to_est))
    model <- SelectModel(models, vec.weights)
    
    
    #-------------------------------------------------------------------
    # GET THE ESTIMATIONS FROM THE MODEL AND THE EXPECTED ERROR
    #-------------------------------------------------------------------
    
    t1 <- Sys.time()
    results <- PredictRFs(model, test[,c(2:dlength)])
    t2 <- Sys.time()
    
    if(sp == "tap-rooted") results <- data.frame(results[,1:4], n_primary = 1, results[5:ncol(results)])
    allresults <- rbind(allresults, results)
    alltest <- rbind(alltest, test)
    
    #-------------------------------------------------------------------
    # COMPARE THE MODELS WITH THE BEST CORRESPONDING  SINGLE LINEAR ESTIMATION
    #-------------------------------------------------------------------
    
    
    # Look at the quality of the estimations
    for(te in to_est){
    
      # Compare it to the Random Forest
      estimation <- results[[te]] 
      truth <- test[[te]]
      
      # Relative error
      rel_diff <- abs((truth - estimation) / truth)
      rel_diff[is.infinite(rel_diff)] <- 0
      diff <- mean(rel_diff)
      
      # R-squared
      fit <- lm(truth ~ estimation)
      
      # Save the data
      
      metaPrec <- rbind(metaPrec, data.frame(var=te,
                                             value=summary(fit)$r.squared,
                                             type="R-SQUARED [-]", 
                                             specie = sp,
                                             degradation = d))
      
      metaPrec <- rbind(metaPrec, data.frame(var=te,
                                             value=diff,
                                             type="RELATIVE ERROR [-]", 
                                             specie = sp,
                                             degradation = d))
    }
  }
}


# Get the best prediction from the direct descriptors
for(sp in species){
  for(te in to_est){
    test <- rs[rs$degradation == 0,]
    
    # Get the best single variable (regression)
    max <- 0
    for(v in vars){
      fit <- lm(test[[te]] ~ test[[v]])
      r2 <- summary(fit)$r.squared
      if(r2 > max){
        max <- r2
      }
    }
    metaPrec <- rbind(metaPrec, data.frame(var=te,
                                           value=max,
                                           type="R-SQUARED [-]", 
                                           specie = sp,
                                           degradation = "-1"))
    # Get the best single variable (RSME)
    min <- 1e9
    for(v in vars){
      estimation <- test[[v]]
      truth <- test[[te]]
      
      rel_diff <- abs((truth - estimation) / truth)
      rel_diff[is.infinite(rel_diff)] <- 0
      diff <- mean(rel_diff, na.rm = T)
      
      if(diff < min){
        min <- diff
      }
    }
    metaPrec <- rbind(metaPrec, data.frame(var=te,
                                           value=min,
                                           type="RELATIVE ERROR [-]", 
                                           specie = sp,
                                           degradation = "-1"))
    
  }
}


metaPrec$degradation[metaPrec$degradation == 0] <- "0 - null"
metaPrec$degradation[metaPrec$degradation == 1] <- "1 - medium"
metaPrec$degradation[metaPrec$degradation == 3] <- "2 - high"
metaPrec$degradation[metaPrec$degradation == -1] <- "direct estimation"

metaPrec$thrs <- 0.1
metaPrec$thrs[grepl("SQU", metaPrec$type)] <- 0.9
metaPrec$max <- 0
metaPrec$max[grepl("SQU", metaPrec$type)] <- 1.5

seps <- c(1:length(unique(metaPrec$var)))+0.5

pos <- position_dodge(0.4)
ggplot(data=metaPrec, aes(x=var, y=value, colour=degradation)) +
  geom_point(size=3, position=pos)+
  #geom_point(colour="white", size=1, position=pos)+  
  facet_grid(specie~type, scale="free") + 
  geom_vline(xintercept = seps, colour="grey", lty=3)+
  geom_hline(aes(yintercept = thrs), lty=2) +
  theme_classic() +
  scale_color_manual(values=leg.col.4)+
  xlab("")+
  ylab("Accuracy evaluation") + 
  theme(axis.text.x=element_text(angle=45, hjust = 1),
        axis.line.x = element_line(colour = 'black', size=0.3, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'),
        panel.background = element_rect(fill = NA, color = "black")) + 
  ggsave(paste0("../figures/rf_accuracy.pdf"), height=6, width=8)+
  ggsave(paste0("../figures/rf_accuracy.png"), height=6, width=8)


remove(temp, fit, diff, estimation)
remove(v, keep, max, r2, te)

#-------------------------------------------------------------------
# PLOT THE REGRESSIONS BETWEEN THE GROUND-TRUTH AND THE MODEL VALUES
#-------------------------------------------------------------------

# For the total root length

x = alltest[["tot_root_length"]]
y1 = allresults[["tot_root_length"]]
y2 = alltest[["length"]]

temp1 <- data.frame(x = x, y = y1, type="RANDOM FOREST", specie = alltest$specie.x, degradation = alltest$degradation)
temp2 <- data.frame(x = x, y = y2, type="DIRECT ESTIMATION", specie = alltest$specie.x, degradation = alltest$degradation)
temp <- rbind(temp2, temp1)

temp$degradation[temp$degradation == 0] <- "0 - null"
temp$degradation[temp$degradation == 1] <- "1 - medium"
temp$degradation[temp$degradation == 3] <- "2 - high"
temp$degradation[temp$degradation == -1] <- "direct estimation"

ggplot(temp, aes(x, y, colour=degradation)) + 
  geom_point() + 
  facet_grid(specie~type) + 
  stat_smooth(method="lm", se=FALSE) + 
  geom_abline(intercept = 0, slope = 1, lty=2) + 
  xlab(paste0("\n Estimated total root length [mm]")) + 
  ylab(paste0("Total root length [mm]\n")) + 
  theme_classic() + 
  scale_colour_manual(values=leg.col.3) +
  theme(text=element_text(size=20),axis.text.x=element_text(angle=45, hjust = 1),
        panel.background = element_rect(fill = NA, color = "black") ,
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
  ggsave(paste0(base.dir,"/00_data/figures/random_vs_direct_length.pdf"), width = 10, height=8)+
  ggsave(paste0(base.dir,"/00_data/figures/random_vs_direct_length.png"), width = 10, height=8)

# For the tip counting

x = alltest[["n_laterals"]]
y1 = allresults[["n_laterals"]]
y2 = alltest[["tip_count"]]

temp1 <- data.frame(x = x, y = y1, type="RANDOM FOREST", specie = alltest$specie.x, degradation = alltest$degradation)
temp2 <- data.frame(x = x, y = y2, type="DIRECT ESTIMATION", specie = alltest$specie.x, degradation = alltest$degradation)
temp <- rbind(temp2, temp1)

temp$degradation[temp$degradation == 0] <- "0 - null"
temp$degradation[temp$degradation == 1] <- "1 - medium"
temp$degradation[temp$degradation == 3] <- "2 - high"
temp$degradation[temp$degradation == -1] <- "direct estimation"


ggplot(temp, aes(x, y, colour=degradation)) + 
  geom_point() + 
  facet_grid(specie~type, scales="free_y") + 
  stat_smooth(method="lm", se=FALSE) + 
  geom_abline(intercept = 0, slope = 1, lty=2) + 
  xlab(paste0("\n Estimated number of roots [-]")) + 
  ylab(paste0("Number of roots [-]\n")) + 
  theme_classic() + 
  scale_colour_manual(values=leg.col.3) +
  theme(text=element_text(size=20),axis.text.x=element_text(angle=45, hjust = 1),
        panel.background = element_rect(fill = NA, color = "black") ,
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
  ggsave(paste0(base.dir,"/00_data/figures/random_vs_direct_tips.pdf"), width = 10, height=8)+
  ggsave(paste0(base.dir,"/00_data/figures/random_vs_direct_tips.png"), width = 10, height=8)




# Make the plots for all the regressions



pdf(paste0("../figures/rf_regressions.pdf"))
par(mfrow=c(4,4), mar=c(4,4,4,3))
for(te in to_est){
  tryCatch({
    x = test[[te]]
    y = results[[te]]
    
    fit = lm(y~x)
    #Prediction intervals
    newx <- data.frame(x=seq(0, max(x), length.out=500))
    preds <-  predict(fit, newdata =newx, interval="prediction")
    error <- data.frame("value"=((preds[,1] - preds[,2]) / max(preds[,1]))*100, "x"= newx)
    #Confidence intervals
    fitted.values = preds[,1]
    pred.lower = preds[,2]
    pred.upper = preds[,3]
    
    title <- paste0(te, " \n r2 = ", round(summary(fit)$r.squared,2))
    plot(x,y, main=title, col="#00000050", xlab="Ground-truth", ylab="Prediction")
    polygon(c(rev(newx$x), newx$x), c(rev(pred.upper), pred.lower), col = '#dc021880', border = NA)
    abline(a = 0, b = 1, col="blue", lwd=2)
    abline(fit, lwd=2, lty=2, col="red")
  }, warning = function(w) {
    message(w)
  }, error = function(e) {
    message(e)
  })
}
remove(te, fit, title, preds, error)
dev.off()


