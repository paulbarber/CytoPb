# CytoPb: Script 0
# P R Barber, Jan 2023

# read raw data from IMC txt files, make 16-bit tiff and panel.csv

#BiocManager::install("imcRtools")
suppressMessages(library(imcRtools))
suppressMessages(library(EBImage))
suppressMessages(library(strex))

# Image scale
image_scale_umperpixel = 1  # um per pixel


if(!exists("working_folder")){
  working_folder <- choose.dir(caption = "Select data folder")
}

cat(paste("CytoPb 0 Working in: ", working_folder, "\n"))

# where are the raw txt files, they can be in subfolders
raw_folder <- paste0(working_folder, "/raw")
if(!dir.exists(raw_folder)){
  stop("There is no 'raw' data folder")
}

# define file names
panel_filename <- paste0(working_folder, "/panel.csv")
panel_names <- NULL

# create output folders
img_folder <- paste0(working_folder, "/img/")
dir.create(img_folder, showWarnings = F)


# Get names of raw images, recursively into all folders
raw_filenames <- list.files(raw_folder, 
                            pattern = "*.txt$", 
                            full.names = T,
                            recursive = T)

# Loop over all the files
for(i in 1:length(raw_filenames)){
  
  this_filename <- basename(raw_filenames[i])
  this_folder <- dirname(raw_filenames[i])
  
  # Read in the raw data
  raw <- NULL
  try(raw <- readImagefromTXT(this_folder, pattern = this_filename))
  if(is.null(raw)){
    cat(paste("Skipping:   ", this_filename, " (not an image file)\n"))
    next
  } else {
    cat(paste("Processing: ", this_filename, "\n"))
  }
  
  image_name <- names(raw)[1]
  img <- raw@listData[[1]]

  channel <- sub("Di", "", dimnames(img)[[3]])
  
  # There is Target information in the first line of the text file
  con <- file(paste0(this_folder, "/", image_name, ".txt"), "r")
  header <- readLines(con, n = 1)
  close(con)
  cols <- strsplit(header, "\t")
  
  # which col has the first channel? Should only be one!
  first <- grep(channel[1], cols[[1]])
  last <- first + length(channel) - 1
  # pull out relevant col names
  cols <- cols[[1]][first:last]
  
  # Some will have a '-' and these are labelled
  named <- grep("-", cols)
  names <- sub("^[0-9]{2,3}[A-Za-z]{1,2}-([-_a-zA-Z0-9]+) *\\([a-zA-Z0-9]+\\)$", "\\1", cols)
  
  panel <- data.frame(channel)
  panel$name <- channel
  panel$name[named] <- names[named]  # replace channel with the names we know
  panel$keep <- 0
  panel$keep[named] <- 1   # mark those with names to keep
  
  # keep the panel
  assign(paste0("panel", i), panel)
  panel_names <- c(panel_names, image_name)
  
  # Save the channels to keep in a tiff in the img folder
  filename <- paste0(img_folder, image_name, ".tiff")
  img_keep <- img[,,named]
  writeImage(img_keep/2^16, filename, bits.per.sample = 16)
  
  # NB these are 16bit images with original scale, Steinbock output are 32bit (unknown range).
  # If they are saved without scaling, they somehow get thresholded and become binary.
  # This seems to be a problem with EBImage.
}

rm(this_filename, this_folder)

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
rm(raw, img_keep)

cat("If your image names are long, now is a good time to shorten them.\n")

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


