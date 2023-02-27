###########################################################################################
#
#  this script extracts spectral reflectance (and quality control layer if preferred) from
#  original modis brdf-corrected reflectance hdf files downloaded from DACC, convert hdf 
#  to tif and merge spectral bands. the output is a single raster file for each hdf file
#  that contains all desired spectral bands and quality control layers
#
#    --- Last updated:  2023.02.16 By Julia Mei Hanzl
###########################################################################################

#******************** close all devices and delete all variables *************************#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?
#*****************************************************************************************#

#****************************** load required libraries **********************************#
### install and load required R packages
list.of.packages <- c("ggplot2", "rgdal", "terra")  
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# load libraries
invisible(lapply(list.of.packages, library, character.only = TRUE))
#*****************************************************************************************#

#************************************ user parameters ************************************#
# define an output directory to store outputs
out.dir <- "/Users/darylyang/Desktop/gdrive/github/lai_scaling/data"
# create output directory if not exist
if (! file.exists(out.dir)) dir.create(out.dir,recursive=TRUE)

out.crs <- "+proj=utm +zone=3 +datum=WGS84 +units=m +no_defs"
#*****************************************************************************************#

#*************************************** load data ***************************************#
# load in leaf area index (lai) data
lai.tab.dir <- "/Users/darylyang/Desktop/gdrive/github/lai_scaling/data/lai_data"
lai.tab.lst <- list.files(lai.tab.dir, pattern = '2022_Teller_LAI.csv', full.names = TRUE
                          , recursive = FALSE)
lai.tab <- read.csv(lai.tab.lst, header = TRUE)

# load in lai polygon shp file
roi.stats.dir <- "/Users/darylyang/Desktop/gdrive/github/lai_scaling/output/roi_stats_v1_clean"

#*****************************************************************************************#
# extract roi list from lai data table
roi.ids <- lai.tab$ROI_ID
roi.stat.combn <- data.frame()
for (roi.id in roi.ids)
{
  roi.stat.dir <- paste0(roi.stats.dir, '/', 'lai_roi_', roi.id, '.csv')
  # if file exist, then start processing
  if (file.exists(roi.stat.dir))
  {
    # read in roi stat data
    roi.stat.org <- read.csv(roi.stat.dir)
    roi.stat.use <- roi.stat.org[, -1]
    # caclulate mean of all uas flighst
    roi.stat.mean <- colMeans(roi.stat.use)
    va.names <- names(roi.stat.mean)
    roi.stat.mean <- unname(roi.stat.mean)
    roi.stat.mean <- data.frame(t(roi.stat.mean))
    names(roi.stat.mean) <- va.names
    roi.stat.combn <- rbind(roi.stat.combn, roi.stat.mean)
  }
  # file does not exist, then use NA for current ROI
  if (! file.exists(roi.stat.dir))
  {
    roi.stat.mean <- data.frame(t(rep(0, 10)))
    names(roi.stat.mean) <- va.names
    roi.stat.combn <- rbind(roi.stat.combn, roi.stat.mean)
  }
}

roi.stat.combn <- data.frame(roi.stat.combn)
roi.stat.combn[roi.stat.combn == 0] <- NA

ndvi <- (roi.stat.combn$nir_ave-roi.stat.combn$red_ave)/(roi.stat.combn$nir_ave+roi.stat.combn$red_ave)
dvi <-  roi.stat.combn$nir_ave-roi.stat.combn$red_ave
gci <- roi.stat.combn$nir_ave/roi.stat.combn$green_ave - 1
gdvi <- roi.stat.combn$nir_ave-roi.stat.combn$green_ave 
gndvi <- (roi.stat.combn$nir_ave-roi.stat.combn$green_ave)/(roi.stat.combn$nir_ave+roi.stat.combn$green_ave)
NIRv <- ndvi*roi.stat.combn$nir_ave
savi <- (1.5*(roi.stat.combn$nir_ave-roi.stat.combn$red_ave))/(roi.stat.combn$nir_ave+roi.stat.combn$red_ave + 0.5)
sr <- roi.stat.combn$nir_ave/roi.stat.combn$red_ave


# combine lai data with roi stats
roi.lai.data.combn <- cbind(lai.tab, roi.stat.combn, ndvi, dvi, gci, gdvi, NIRv, savi, sr)

out.name <- paste0(out.dir, '/', 'teller_lai_database_v2', '.csv')
write.csv(roi.lai.data.combn, out.name, row.names = FALSE)









