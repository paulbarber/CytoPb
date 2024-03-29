# CytoPb: Script 0
# P R Barber, Mar 2023

# read raw data from CODEX tiff files, make tiff and panel.csv
# This is the format of the data provided by 
# Schürch et. al: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7479520/

suppressMessages(library(cytomapper))
suppressMessages(library(strex))

# Image scale
#image_scale_umperpixel = 0.37744  # um per pixel
image_scale_umperpixel = 0.75488  # um per pixel in the montage files

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

# create output folders
img_folder <- paste0(working_folder, "/img/")
dir.create(img_folder, showWarnings = F)


# Get names of raw images
raw_filenames <- list.files(raw_folder, pattern = "*.tif", full.names = T)
#raw_names <- str_before_last_dot(list.files(raw_folder, pattern = "*.tif", full.names = F))

# read in the channel names
cn <- read.delim(channelNames_filename, header = F)
channel_descriptions <- str_after_last(cn$V1, " - ")
channel_names <- str_before_first(cn$V1, " - ")
# keep channel names without a description as well
channel_names[is.na(channel_names)] <- cn$V1[is.na(channel_names)]
# Can also mark those to keep, not marked "blank" or "empty"
# also just keep the first HOECHST channel
keep <- rep(1, length(channel_names))
keep[grep("blank", channel_names)] <- 0
keep[grep("empty", channel_names)] <- 0
keep[grep("HOECHST", channel_names)] <- 0
keep[grep("HOCHST", channel_names)] <- 0   # sometimes misspelled
keep[1] <- 1  # the first HOECHST channel

# make the panel.csv table
channel = 1:length(channel_names)
panel <- data.frame(channel)
panel$name <- channel_names
panel$keep <- keep
panel$description <- channel_descriptions

# Loop over all the files
pb = txtProgressBar(min = 0, max = length(raw_filenames), initial = 0)
for(i in 1:length(raw_filenames)){
  
  # Load raw images
  raw <- loadImages(raw_filenames[i])
  image_name <- names(raw)[1]
  img <- raw@listData[[1]]
  
  setTxtProgressBar(pb,i)
  
  # check the channel numbers
  nChannels <- dim(img)[3]
  if(nChannels != length(keep)){
    stop("Number of channels in the images and ChannelNames do not match.")
  }

  # Save the channels to keep in a tiff in the img folder
  filename <- paste0(img_folder, image_name, ".tiff")
  img_keep <- img[,,(keep==1)]
  writeImage(img_keep, filename, bits.per.sample = 16)
  # NB these are 16bit images
  rm(raw)
}
close(pb)

# write the panel file
write.csv(panel, panel_filename, row.names = F)

rm(img)
rm(img_keep)


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


