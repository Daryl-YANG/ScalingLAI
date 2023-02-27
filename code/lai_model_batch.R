###########################################################################################
#
#  this script extracts spectral reflectance (and quality control layer if preferred) from
#  original modis brdf-corrected reflectance hdf files downloaded from DACC, convert hdf 
#  to tif and merge spectral bands. the output is a single raster file for each hdf file
#  that contains all desired spectral bands and quality control layers
#
#    --- Last updated:  2021.06.22 By Daryl Yang <dediyang@bnl.gov>
###########################################################################################

#******************** close all devices and delete all variables *************************#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?
#*****************************************************************************************#

#****************************** load required libraries **********************************#
### install and load required R packages
list.of.packages <- c("ggplot2", "randomForest", "caTools", "ggpmisc")  
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# load libraries
invisible(lapply(list.of.packages, library, character.only = TRUE))
#*****************************************************************************************#

#************************************ user parameters ************************************#
# define an output directory to store outputs
our.dir <- "/Users/darylyang/Desktop/gdrive/github/lai_scaling/figures"
# create output directory if not exist
if (! file.exists(our.dir)) dir.create(our.dir,recursive=TRUE)
# creat an temporary to store files temporarily generated during the course of processing
temp.dir <- file.path(paste0(our.dir, "/", 'temporary'))
if (! file.exists(temp.dir)) dir.create(temp.dir,recursive=TRUE)

# define how many models you want to build?
nmodel <- 100
#*****************************************************************************************#

#*************************************** load data ***************************************#
# define the directory to hdf files
source.dir <- "/Users/darylyang/Desktop/gdrive/github/lai_scaling/data"
data.dir <- list.files(source.dir, pattern = 'teller_lai_database_v2.csv',
                       full.names = TRUE)
data.org <- read.csv(data.dir)
str(data.org)
#data.org <- na.omit(data.org)
data.org <- data.org[, c(9, 12, 13, 14, 17, 18, 19, 20:27)]
#data.org$LAI[data.org$LAI > 5] <- NA
data.org <- na.omit(data.org)
#*****************************************************************************************#

#************************************* biomass model *************************************#
### split data into training and validation
data.split <- sample.split(data.org, SplitRatio = 0.6) 
# get training and validation data
rf.data.train <- subset(data.org, data.split == "TRUE") 
rf.data.test <- subset(data.org, data.split == "FALSE") 
#rf.data.train <- rf.data.train
#rf.data.test <- rf.data.train
# extract data for random forest regression
#rf.data.train <- data.train[, c(15:20)]
# build a random forest model 
test.pred <- c()
for (i in 1:nmodel)
{
  # select a portion of training data for each model train
  train.split <- sample.split(rf.data.train, SplitRatio = 0.9) 
  data.train <- subset(rf.data.train, train.split == "TRUE") 
  # build a random forest model
  lai.model <- randomForest(LAI ~ ., data = data.train,
                                ntree = 100, mtry = 4, importance = TRUE)
  # apply model to validation data
  rf.model.test <- predict(lai.model, newdata = rf.data.test)
  # store prediction
  test.pred <- cbind(test.pred, rf.model.test)
}
test.pred <- data.frame(test.pred)
pred.mean <- apply(test.pred, 1, FUN = mean)
pred.sd <- apply(test.pred, 1, FUN = sd)
#*****************************************************************************************#
#*
#************************************** make a plot **************************************#
# make a plot of the predict result
plt.data <- data.frame('truth' = rf.data.test$LAI,
                       'predicted' = pred.mean,
                       'unc' = pred.sd)

rmse <- sqrt(sum((plt.data$predicted-plt.data$truth)^2)/nrow(plt.data))

formula <- y ~ x
ggplot(data = plt.data, aes(x = truth, y = predicted)) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=1.5) +
  geom_point(size = 2, shape = 1) + 
  geom_smooth(method = 'lm', formula = formula, col = 'black') +
  geom_errorbar(aes(ymin = predicted-unc, ymax = predicted+unc)) + 
  ylim(c(0, 5)) + xlim(c(0, 5)) +
  labs(x = 'Ground LAI', y = 'Predicted LAI') +
  theme(legend.position = 'none') +
  theme(axis.text = element_text(size=12), axis.title=element_text(size=13)) +
  stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
               label.x.npc = 0.65, label.y.npc = 0.33,
               eq.with.lhs = "italic(hat(y))~`=`~",
               eq.x.rhs = "~italic(x)",
               formula = formula, parse = TRUE, size = 4, hjust = 0) +
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               label.x.npc = 0.65, label.y.npc = 0.25,
               formula = formula, parse = TRUE, size = 4, hjust = 0) +
  annotate('text', x= 3.33, y = 0.7, label = 'RMSE = 0.85', size = 4, hjust = 0) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pngName = paste0(our.dir, "/",'model_validation_withunc.pdf')
ggsave(pngName, plot = last_plot(), width = 12, height = 12, units = 'cm')
#*****************************************************************************************#
