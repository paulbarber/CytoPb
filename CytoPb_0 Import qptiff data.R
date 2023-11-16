# CytoPb: Script 0
# P R Barber, Nov 2023

# read raw data from qptiff files, make tiff and panel.csv
# This is the format of the data provided by Akoya systems

# RBioFormats requires Java JDK to be installed

# Make sure JDK >8 is installed
#if (!require("BiocManager", quietly=TRUE)) install.packages("BiocManager")
#BiocManager::install("remotes")
#BiocManager::install("aoles/RBioFormats")


cat("This CytoPb 0 script is under development.
    Your best option is to use Fiji to open the qptiff
    and save the image you want as a tiff.
    And manually create a panel.csv file.\n")
stop()

# The problem is that only the first channel is saved into all channels of the tiff.



# Give lots of memory for Java to be used by RBioFormats
#detach("package:RBioFormats", unload = TRUE) # NB after this, Java remains with old memory limit
# Restart R to change this limit
# Also, NB that cytomapper seems to interfere with RBioFormats, Java?
options( java.parameters = "-Xmx32g" ) # 32GB seems enough for a large image (16 fails on 2nd load!)
library(RBioFormats)
#checkJavaMemory(units = "g")

# Image scale
#image_scale_umperpixel = 0.4988  # um per pixel 10x vectra for resolution=1?
#image_scale_umperpixel = 0.4988*2 # um per pixel 10x vectra for resolution=2?
if(!exists("qptiff_resolution")){
  qptiff_res = 2    # default to 2nd resolution level (about 1um per pixel)
}else{
  qptiff_res = qptiff_resolution
  rm(qptiff_resolution)   # this stops in going into global_data_filename and being overwritten if changed and rerun
}
image_scale_umperpixel = 0.4988*qptiff_res # um per pixel 10x vectra for resolution?


if(!exists("working_folder")){
  working_folder <- choose.dir(caption = "Select data folder")
}

cat(paste("CytoPb 0 Working in: ", working_folder, "\n"))

# where are the raw tif files
raw_folder <- paste0(working_folder, "/raw")
if(!dir.exists(raw_folder)){
  stop("There is no 'raw' data folder")
}

# define file names
channelNames_filename <- paste0(working_folder, "/channelNames.txt")
panel_filename <- paste0(working_folder, "/panel.csv")
panel_names <- NULL

# create output folders
img_folder <- paste0(working_folder, "/img/")
dir.create(img_folder, showWarnings = F)

# Get names of raw images, recursively into all folders
raw_filenames <- list.files(raw_folder, 
                            pattern = "*.qptiff$", 
                            full.names = T,
                            recursive = T)

# Loop over all the files
pb = txtProgressBar(min = 0, max = length(raw_filenames), initial = 0)
for(i in 1:length(raw_filenames)){

  image_name <- strex::str_before_last(basename(raw_filenames[i]), "\\.")
  print(image_name)
    
  # Load raw images from the qptiff stack
  # series = 1 is usually the raw data (other series are slide images etc)
  # resolution = 1 is the highest resolution, specify qptiff_resolution above
  img <- read.image(raw_filenames[i], series = 1, resolution = qptiff_res)
  
  setTxtProgressBar(pb,i)

  # Detect the channel number and names
  g <- globalMetadata(img)
  names_list <- g[grep("^Name", names(g))]
  # order channels by the metadata name
  names_list = names_list[order(names(names_list))]
  names <- unlist(names_list)
  
  # Make a panel file, all we know is the channel names (not the marker)
  panel <- data.frame(channel = names)
  panel$name <- names
  panel$keep <- 1

  # keep the panel
  assign(paste0("panel", i), panel)
  panel_names <- c(panel_names, image_name)
  
  # The image seems to be loaded to fill 0-1
  # We can specify 16-bit tiff to save
  # but we must rescale for the bit depth
  # which we read from the metadata
  #metadataentry <- paste0("BitDepth #", qptiff_res)
  #bitdepth <- as.numeric(g[metadataentry][[1]])
  #if(is.na(bitdepth) | bitdepth<8 | bitdepth>16) bitdepth=16  # check
  # or just use a value because the metadata is not consistent!
  bitdepth <- 12
  
  # Save the channels to keep in a tiff in the img folder
  c <- coreMetadata(img)
  filename <- paste0(img_folder, image_name, ".tiff")
  #writeImage(img, filename, bits.per.sample = 16)
  write.image(img * 2^(bitdepth), filename, force = T, 
              pixelType = "uint16",
              littleEndian = TRUE)

  rm(g,c)
}
close(pb)

# Check panels for all file are identical
# assume panel1 is the one
panel <- panel1
write.csv(panel, panel_filename, row.names = F)

# Check and save any others that are not the same
if(length(panel_names) > 1){
  all_the_same = TRUE
  for(i in 2:length(panel_names)){
    p <- get(paste0("panel", i))
    if(!identical(p, panel)){
      filename <- paste0(working_folder, "/panel_", panel_names[i], ".csv")
      write.csv(p, filename, row.names = F)
      all_the_same = FALSE
    }
  }
  
  if(!all_the_same){
    warning("The raw image files do not all have the same panel.")
  }
}

# remove all the unwanted panel files
rm(list = ls(pattern = "panel+[0-9]"))

# remove this large object from the environment
rm(img)


cat("If your image names are long, now is a good time to shorten them.
    Edit the panel file names column with the names of your markers.\n")

# Example image renaming, modify it for your needs and run shorten_image_names()

shorten_image_names <- function(){
  
  filenames <- list.files(path = img_folder)
  
  for(i in filenames){

    new_name <- paste0(strex::str_before_first(i, "_"),
                       "_",
                       strex::str_after_last(i, "_"))
  
    file.rename(paste0(img_folder, i), paste0(img_folder, new_name))      
    
  }
}

#shorten_image_names()


