# CytoPb run all script

#working_folder <- "D:/images/IMC/data1"
working_folder <- choose.dir(caption = "Select data folder")

#source("CytoPb_0 Import CODEX tif data.R")
source("CytoPb_0 Import IMC txt data.R")
#image_scale_umperpixel = 0.37744  # um per pixel

source("CytoPb_1 Determine Channel Ranges.R")
source("CytoPb_2 Channel Probabilities.R")
source("CytoPb_3 Cell Type Probabilities.R")
