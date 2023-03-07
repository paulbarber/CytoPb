# Overlap Signal Maps: Script 1
# P R Barber, Jan 2023

# Takes in the multi-channel tiff images and the panel.csv file.
# The tiff is expected to just have the "keep" channels
# The panel.csv is expected to have "ALL" the channels. 
#   Columns "channel" "name" "keep"
# Can define channels_needed here to reduce processing later.
# Estimates the foreground of each image.
# Save image to channel_png folder
# Save convenient tables of channels with image numbers.
# output: Channel Lists.txt
# Save table of positive and negative values for each channel of each image.
# output: pos_value_table.csv neg_value_table.csv
# Save the complete workspace for follow on scripts.


#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("cytomapper")

library(cytomapper)
library(ggplot2)

# REMEMBER TO SET WORKING DIRECTORY
#setwd("...")

# Set the scale of the blurring in pixels, usually 5 pixels (5 um)
sigma = 3

# User check of working directory.
print("Working in:")
print(getwd())
print("Enter 'y' to proceed:")
proceed = readLines(n=1)
stopifnot(proceed == "y")

# image and panel file locations
image_location <- "img"
panel_location <- "panel.csv"

# read in channel names
panel <- read.csv(panel_location)
panel_keep <- subset(panel, keep == 1)
panel_keep$image_number <- 1:dim(panel_keep)[1]

# Channels needed for identification
#channels_needed <- c("CD3", "CD4", "CD8", "CD19", "CD25", "FOXP3", "CD34", "CD31", "CD45")
channels_needed <- panel_keep$name
panel_needed <- panel_keep[panel_keep$name %in% channels_needed, c("image_number", "name")] 

# folder to save channel QC images to
folder <- "channel_png"
dir.create(folder, showWarnings = F)

sink(file = "Channel Lists.txt")
print(paste("All channel names in order are:"))
print(panel_keep)
print(paste("Channels needed for cell identification are:"))
print(panel_needed)
sink(file = NULL)

# load images
images <- loadImages(image_location)


# Account for different ways things may have been stored
# Check number of channels in all images
nChannels <- vector()
for(i in 1:length(images)){
  c <- dim(images[i]@listData[[1]])[3]
  nChannels <- c(nChannels, c)
}
nChannels <- unique(nChannels)
stopifnot(length(nChannels)==1)

# Check number of rows in panel.csv
nRows <- dim(panel)[1]
nKeep <- dim(panel_keep)[1]

# If more entries in the keep panel than the images
# or if more entries in the images than the keep panel
if(nKeep > nChannels){
  stop("There are not enough channels in the images for those in the keep column of panel.csv.")
} else if(nKeep < nChannels){
  # All image channels must be identified in the panel
  if(nRows != nChannels){
    stop("Image channels are not all or uniquely identified in the panel.csv file.")
  }
  
  # extract the image channels that we need
  for(i in 1:length(images)){
    images[i]@listData[[1]] <- images[i]@listData[[1]][,,which(panel$keep==1)]
  }

}
# Otherwise the number of image channels must be the keep channels


channelNames(images) <- panel_keep$name

# Identify negative and positive channel values
estimateNegPosValue <- function(image, channel, sigma = 10){
  
  # image is a Cytoimagelist
  # channel is a name of a channel
  
  # Get EBImage
  img <- image@listData[[1]][,,channel]
  img_blur <- gblur(img, sigma = sigma) # to get the final intensity from
  
  img2 <- medianFilter(img/max(img), 2)
  img3 <- img2 > otsu(img2, range = c(0, max(img2)))
  img4 <- opening(img3, kern = makeBrush(7, shape='disc'))  # to determine the foreground
  
  
  image_name <- names(image)[1]
  filename <- paste0(folder, "/", image_name, "_", channel, "_fgmask.png")
  writeImage(img4, filename)
  
  
  # calculate negative and positive image values in the original image units
  # suppressWarnings on max in case there is no foreground
  posv <- suppressWarnings((mean(img_blur[img4==1])))
  negv <- suppressWarnings((mean(img_blur[img4==0])))
  
  # If posv comes out 0 (max returns -Inf), there was no foreground = NA
  # If posv comes out NaN (mean returns NaN), there was no foreground = NA
  if(is.nan(posv)) posv = NA
  else if(posv <= 0) posv = NA
  if(is.nan(negv)) negv = NA
  else if(negv <= 0) negv = NA
  
  # return
  c(negv, posv)
}

# estimate the negative and positive values for every channel (in panel_keep) of every image
# write to csv file to check/store, in Excel
# If table already there then no need to do this
if(!file.exists("pos_value_table.csv")){
  image_list <- names(images)
  b <- rep(" ", dim(panel_keep)[1])
  pos_table <- data.frame(panel_keep$name, 
                          panel_keep$image_number, 
                          matrix(nrow = dim(panel_keep)[1], ncol = length(image_list)))
  names(pos_table) <- c("Channel", "No.", image_list)
  neg_table <- pos_table
  
  pb = txtProgressBar(min = 0, max = length(images), initial = 0)
  for(i in 1:length(images)){
    
    image_name <- names(images)[i]
    setTxtProgressBar(pb,i)
    
    for(channel in channelNames(images)){
      
      vals <- estimateNegPosValue(images[i], channel, sigma = sigma)
      
      neg_table[which(neg_table$Channel == channel), image_name] = vals[1] 
      pos_table[which(pos_table$Channel == channel), image_name] = vals[2] 
      
    }
  }
  close(pb)
}

# Write blank files, only if file does not already exist since good work could be overwritten!
if(!file.exists("pos_value_table.csv")) write.csv(pos_table, file = "pos_value_table.csv")
if(!file.exists("neg_value_table.csv")) write.csv(neg_table, file = "neg_value_table.csv")

# Save everything so far
#save.image(file = "Overlap Signal Maps.RData")

