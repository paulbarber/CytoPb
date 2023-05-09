# CytoPb run all script
# Best to copy this into the working folder with your data
# So you can edit for your data set.
# Raw data in a "raw" subfolder, tiff images will be in a "img" folder etc.

# Specify or select a working folder
#working_folder <- "D:/images/IMC/data1"
working_folder <- choose.dir(caption = "Select data folder")

# Optionally run an import script
#source("CytoPb_0 Import CODEX tif data.R")
#source("CytoPb_0 Import IMC txt data.R")

# Specify options
#image_scale_umperpixel = 0.75488  # um per pixel, default is 1
use_global_ranges = TRUE   # can use a global range for each channel or individual for each image, default is TRUE (in 2)
probability_threshold = 0.1   # Threshold for "high probability", default is 0.5 (in 3)

# Choose which processing to do
source("CytoPb_1 Determine Channel Ranges.R")
source("CytoPb_2 Channel Probabilities.R")
source("CytoPb_3 Cell Type Probabilities.R")

# Neighbourhood processing
source("CytoPb_4 Cell Neighbourhoods.R")

# Cell type interaction processing
cell_type_pairs <- list(c("CD8Tcell", "Progenitor"), c("Teffcell", "Progenitor"))
source("CytoPb_4a Cell Type Colocalisations.R")


# Single channel optimisation 
#TEST_specific_channel <- "CD11b"
#TEST_specific_image <- "ROI_001"
#source("CytoPb_2 Channel Probabilities.R")
