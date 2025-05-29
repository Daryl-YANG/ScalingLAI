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
list.of.packages <- c("ggplot2", "randomForest", "caTools", "ggpmisc", "terra", 
                      "foreach", "doParallel")  
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# load libraries
invisible(lapply(list.of.packages, library, character.only = TRUE))
#*****************************************************************************************#

#************************************ user parameters ************************************#
# define an output directory to store outputs
out.dir <- "/Volumes/Mercury/projects/lai_scaling/output/maps"
# create output directory if not exist
if (! file.exists(out.dir)) dir.create(out.dir,recursive=TRUE)
# creat an temporary to store files temporarily generated during the course of processing
temp.dir <- file.path(paste0(out.dir, "/", 'temporary'))
if (! file.exists(temp.dir)) dir.create(temp.dir,recursive=TRUE)

# define how many models you want to build?
nmodel <- 100

# define the output resolution of biomass map
reso <- 1 # m
#*****************************************************************************************#

#*************************************** load data ***************************************#
# define the directory to hdf files
source.dir <- "/Volumes/Mercury/GitHub/lai_scaling/data"
data.dir <- list.files(source.dir, pattern = 'teller_lai_database_v2.csv',
                       full.names = TRUE)
data.org <- read.csv(data.dir)
str(data.org)
#data.org <- na.omit(data.org)
data.org <- data.org[, c(9, 12, 13, 14, 17, 18, 19, 20:27)]
#data.org$LAI[data.org$LAI > 5] <- NA
data.org <- na.omit(data.org)

# load in uas datasets
# chm data
uas.dir <- "/Volumes/Mercury/projects/lai_scaling/data"
chm.dir <- list.files(uas.dir, pattern = '*teller_dji_chm_combined_v2.tif',
                      full.names = TRUE, recursive = TRUE)
chm.rst <- terra::rast(chm.dir)
# solo refl data
refl.dir <- list.files(uas.dir, pattern = '*SewPen_SOLO_20220713_Teller_MM27_Flight14_MBS',
                       full.names = TRUE, recursive = TRUE)
refl.rst <- terra::rast(refl.dir)
names(refl.rst) <- c('green', 'red', 'rededge', 'nir', 'none')
refl.rst <- c(refl.rst$green, refl.rst$red, refl.rst$rededge, refl.rst$nir)
#*****************************************************************************************#

#************************************* prepare data **************************************#
# create a mask to mask out non value regions
mask <- refl.rst$nir
mask[mask > 60000] <- NA
mask[mask < 60000] <- 1
# clip and resample chm data to solo refl region
chm.rst.clp <- terra::crop(chm.rst, refl.rst)
chm.rst.rsp <- resample(chm.rst.clp, refl.rst)
# stack refl and chm data
rst.stack <- c(refl.rst, chm.rst.rsp)
rst.input <- terra::mask(rst.stack, mask)

# calculate vegetation index
ndvi <- (rst.input$nir - rst.input$red)/(rst.input$nir + rst.input$red)
dvi <- rst.input$nir - rst.input$red
gci <- rst.input$nir/rst.input$green - 1
gdvi <- rst.input$nir - rst.input$green
gndvi <- (rst.input$nir - rst.input$green)/(rst.input$nir + rst.input$green)
NIRv <- ndvi * rst.input$nir
savi <- (1.5*(rst.input$nir-rst.input$red))/(rst.input$nir+rst.input$red + 0.5)
sr <- rst.input$nir / rst.input$red

# combine vegetation index raster with reflectance and chm
rst.input.add <- c(rst.input, ndvi, dvi, gci, gdvi, NIRv, savi, sr)

# calculate grid mean and sd
window.size <- reso/xres(rst.input.add)
grid.mean <- aggregate(rst.input.add, window.size, mean, na.rm = TRUE)
names(grid.mean) <- c('green', 'red', 'rededge', 'nir', 'chm', 'ndvi',
                      'dvi', 'gci', 'gdvi', 'NIRv', 'savi', 'sr')
grid.sd <- aggregate(rst.input.add, window.size, sd, na.rm = TRUE)
names(grid.sd) <- c('green', 'red', 'rededge', 'nir', 'chm', 'ndvi',
                      'dvi', 'gci', 'gdvi', 'NIRv', 'savi', 'sr')

rst.use <- c(grid.mean$green, grid.mean$red, grid.mean$rededge, grid.mean$nir,
             grid.sd$green, grid.sd$red, grid.sd$rededge, grid.sd$nir,
             grid.mean$chm, grid.sd$chm, grid.mean$ndvi, grid.mean$dvi,
             grid.mean$gci, grid.mean$gdvi, grid.mean$NIRv, grid.mean$savi,
             grid.mean$sr)
#*****************************************************************************************#


#************************************** apply model **************************************#
# convert raster to dataframe for apply random forest model
stack.dataframe <- as.data.frame(rst.use, xy=TRUE)
coords <- stack.dataframe[, c(1,2)]
predictors <- stack.dataframe[, c(4:6, 9:19)]
names(predictors) <- names(data.org)[2:15]
predictors[is.na(predictors)] <- NA
predictors[mapply(is.infinite, predictors)] <- NA

rm(rst.input.add, rst.use, grid.mean, grid.sd)
### apply model on application data
#Setup backend to use many processors
totalCores = detectCores()
#Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1]-1) 
registerDoParallel(cluster)

rf.data.train <- data.org
result <- foreach(i=1:100) %dopar% {
  train.split <- caTools::sample.split(rf.data.train, SplitRatio = 0.9)
  data.train <- subset(rf.data.train, train.split == "TRUE")
  lai.model <- randomForest::randomForest(LAI ~ ., data = data.train,
                                              ntree = 100, mtry = 4, importance = TRUE)
  lai.pred <- predict(lai.model, newdata = predictors)
}
pred.df <- data.frame(result)

### calculate mean biomass prediction
pred.mean <- apply(pred.df, 1, FUN = mean, na.rm = TRUE)
# convert biomass dataframe to raster
pred.mean <- cbind(coords, pred.mean)
lai.pred.rst <- terra::rast(pred.mean, type = 'xyz')
crs(lai.pred.rst) <- crs(refl.rst)
basename <- basename(refl.dir)
outname <- paste0(out.dir, '/', 'parallel_lai_', basename)
writeRaster(lai.pred.rst, outname, overwrite=TRUE)

### calculate biomass prediction uncertainty
pred.unc <- apply(pred.df, 1, FUN = sd, na.rm = TRUE)
# convert biomass dataframe to raster
pred.unc <- cbind(coords, pred.unc)
lai.unc.rst <- terra::rast(pred.unc, type = 'xyz')
crs(lai.unc.rst) <- crs(refl.rst)
basename <- basename(refl.dir)
outname <- paste0(out.dir, '/', 'parallel_lai_unc_', basename)
writeRaster(lai.unc.rst, outname, overwrite=TRUE)
#*****************************************************************************************#















