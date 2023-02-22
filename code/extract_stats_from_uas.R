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
out.dir <- "/Users/darylyang/Desktop/gdrive/github/lai_scaling/output/roi_stats"
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
lai.roi.dir <- "/Users/darylyang/Desktop/gdrive/github/lai_scaling/data/lai_data/2022_Teller_LAI_ROI_Shapefile"
lai.roi.lst <- list.files(lai.roi.dir, pattern = '2022_Teller_LAI_ROIs.shp',
                                    full.names = TRUE, recursive = TRUE)
lai.roi <- readOGR(lai.roi.lst)
# reproject the roi file to the projection system of uas imagery
lai.roi <- spTransform(lai.roi, out.crs)


# load in UAS data list
uas.img.dir <- "/Volumes/data1/projects/lai_scaling/data/solo_data"
uas.img.lst <- list.files(uas.img.dir, pattern = '*MBS.tif', full.names = TRUE,
                           recursive = TRUE)

# load in UAS canopy height data
dji.chm.dir <- "/Volumes/data1/projects/lai_scaling/data/dji_chm"
dji.chm.lst <- list.files(dji.chm.dir, pattern = '*CHM.tif', full.names = TRUE,
                           recursive = TRUE)
#*****************************************************************************************#

#************************************* main function *************************************#
# extract roi list from lai data table
roi.ids <- lai.tab$ROI_ID
for (roi.id in roi.ids)
{
  # extract polygon for current processing ROI
  roi.shp <- subset(lai.roi, id == roi.id)
  # search for this polygon in all uas images
  refl.tab <- c()
  for (uas.img.file in uas.img.lst)
  {
    # load in current uas imagery
    uas.img.rst <- terra::rast(uas.img.file)
    # crop the polygon region from the uas imagery
    roi.img.rst <- try(terra::crop(uas.img.rst, roi.shp, mask = TRUE), silent = TRUE)
    # if the two overlapped then do the calculation
    if (class(roi.img.rst) == 'SpatRaster')
    {
      # convert roi raster data into data frame
      roi.dataframe <- as.data.frame(roi.img.rst)
      roi.dataframe <- roi.dataframe[, c(1:4)]
      roi.dataframe <- na.omit(roi.dataframe)
      names(roi.dataframe) <- c('green', 'red', 'rededge', 'nir')
      
      # calculate the mean reflectance of each band
      mean.refl <- apply(roi.dataframe, 2, FUN = mean, na.rm = TRUE)
      sd.refl <- apply(roi.dataframe, 2, FUN = sd, na.rm = TRUE)

      uas.flight.name <- basename(uas.img.file)
      plot.refl <- data.frame(uas.flight.name, t(mean.refl), t(sd.refl))
      names(plot.refl) <- c('uas_flight', 'green_ave', 'red_ave', 'rededge_ave', 'nir_ave',
                            'green_sd', 'red_sd', 'rededge_sd', 'nir_sd')
      # save the data into the tab
      refl.tab <- rbind(refl.tab, plot.refl)
    }
  }
  
  # extract canopy height info
  chm.tab <- c()
  for (dji.chm.file in dji.chm.lst)
  {
    # load in current dji chm raster
    dji.chm.rst <- terra::rast(dji.chm.file)
    roi.chm.rst <- try(terra::crop(dji.chm.rst, roi.shp, mask = TRUE), silent = TRUE)
    if (class(roi.chm.rst) == 'SpatRaster')
    {
      chm.dataframe <- as.data.frame(roi.chm.rst)
      mean.chm <- apply(chm.dataframe, 2, FUN = mean, na.rm = TRUE)
      sd.chm <- apply(chm.dataframe, 2, FUN = sd, na.rm = TRUE)
      
      plot.chm <- data.frame(mean.chm, sd.chm)
      names(plot.chm) <- c('chm_ave', 'chm_sd')
      # save the data into the tab
      chm.tab <- rbind(chm.tab, plot.chm)
    }
  }
  
  mean.plot.chm.stats <- apply(chm.tab, 2, FUN = mean, na.rm = TRUE)
  mean.plot.chm.stats <- data.frame(t(mean.plot.chm.stats))
  
  if(is.null(refl.tab) == 'FALSE')
  {
    chm.add <- mean.plot.chm.stats[rep(1, nrow(refl.tab)),]
    names(chm.add) <- c('chm_ave', 'chm_sd')
    #chm.add <- data.frame('chm_ave' = chm.rep[, 1], 'chm_sd' = chm.rep[, 2])
    plot.data.combn <- cbind(refl.tab, chm.add)
    out.name <- paste0(out.dir, '/', 'lai_roi_', roi.id, '.csv')
    write.csv(plot.data.combn, out.name, row.names = FALSE)
  }
}
#*****************************************************************************************#