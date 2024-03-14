# CytoPb: Script 0
# P R Barber, C Eke, Mar 2024

# read single plane tif files from a folder that contains all planes.
# The folder name will form the new multi-plane tif image file name.
# We expect the images to be named with the channel name.
# Assume images end ".tiff"
# Also, code below saves as 32-bit integer which will probably not preserve the 
# original image scale. It was written for the output of IMCDenoise which saves 
# floating point images of unknown dynamic range.

suppressMessages(library(EBImage))
suppressMessages(library(strex))

# Image scale, can be changed in the parent script
image_scale_umperpixel = 1.0  # um per pixel, expecting IMC images

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
panel_filename <- paste0(working_folder, "/panel.csv") # to be created

# create output folders
img_folder <- paste0(working_folder, "/img/")
dir.create(img_folder, showWarnings = F)


# get a list of the raw data folders
folder_names <- list.dirs(raw_folder, recursive = F)

panel_names <- vector()
j = 1 
for(i in 1:length(folder_names)){
  
  folder <- folder_names[i]
  image_names <- list.files(path = folder, pattern = "*.tif")
  
  if(length(image_names) < 1) next    # skip empty folders
  
  # store all panels
  assign(paste0("panel", j), strex::str_before_last(image_names, ".tif"))
  j = j + 1
  
  # keep list of panel names
  panel_names <- c(panel_names, strex::str_after_last(folder, "/"))
}

# Check panels for all file are identical
# assume panel1 is the one
# make the panel.csv table
channel = 1:length(panel1)
panel <- data.frame(channel)
panel$name <- panel1
panel$keep <- 1
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

# delete all unwanted panels
rm(list = ls(pattern = "panel[0-9]"))

# read in the data and save in new format
for(i in 1:length(folder_names)){
  
  folder <- folder_names[i]
  image_paths <- list.files(path = folder, pattern = "*.tif", full.names = T)

  print(folder)
  
  if(length(image_paths) < 1) next    # skip empty folders
  
    # somewhere to store the planes for this image
  composite_image = NULL
  
  pb = txtProgressBar(min = 0, max = length(image_paths), initial = 0)
  
  for(j in 1:length(panel$name)){
    
    setTxtProgressBar(pb, j)

    image_path <- paste0(folder, "/", panel$name[j], ".tiff")
    
    # Load raw images
    raw <- loadImages(image_path)
    img <- raw@listData[[1]]     # assuming from the start there is only one plane
    
    composite_image <- EBImage::combine(composite_image, img)
  }
  
  filename <- paste0(strex::str_after_last(folder, "/"), ".tiff")
  path <- paste0(img_folder, filename)
  
  # Will preserve the image scale
  #writeImage(composite_image/2^16, path, bits.per.sample = 16)
  
  # OR Will extend the dynamic range to cope with floating point input images
  writeImage(composite_image/2^16, path, bits.per.sample = 32)
  
  rm(composite_image)
}

close(pb)
rm(img)


cat("If your image names are long, now is a good time to shorten them.")

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


