# CytoPb run all script

# Specify or select a working folder
#working_folder <- "D:/images/IMC/data1"
working_folder <- choose.dir(caption = "Select data folder")

# Optionally run an import script
#source("CytoPb_0 Import CODEX tif data.R")
#source("CytoPb_0 Import IMC txt data.R")

# Specify options
#image_scale_umperpixel = 0.75488  # um per pixel, default is 1
use_global_ranges = TRUE   # can use a global range for each channel or individual for each image, default is FALSE (in 2)
probability_threshold = 0.1   # Threshold for "high probability", default is 0.5 (in 3)

# Choose which processing to do
source("CytoPb_1 Determine Channel Ranges.R")
source("CytoPb_2 Channel Probabilities.R")
source("CytoPb_3 Cell Type Probabilities.R")
