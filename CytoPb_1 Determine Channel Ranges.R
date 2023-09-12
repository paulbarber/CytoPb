# CytoPb: Script 1 Determine Channel Ranges
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

suppressMessages(library(cytomapper))
suppressMessages(library(ggplot2))
suppressMessages(library(strex))

# source the helper functions
source("R/estimateNegPosValue.R")
source("R/plotPosValues.R")

# Set the scale of the blurring in pixels, in tests 3 pixels is best (3 um)
# Final results are quite insensitive to this. Appearance of cell maps is sensitive.
if(!exists("image_scale_umperpixel")){
  image_scale_umperpixel = 1
  }
sigma = ceiling(3 / image_scale_umperpixel)


if(!exists("working_folder")){
  working_folder <- choose.dir(caption = "Select data folder")
}

cat(paste("CytoPb 1 Working in: ", working_folder, "\n"))

global_data_filename <- paste0(working_folder, "/CytoPb.RData")

# image and panel file locations
image_location <- paste0(working_folder, "/img")
panel_location <- paste0(working_folder, "/panel.csv")

# File locations
channel_list_filename <- paste0(working_folder, "/Channel Lists.txt")
pos_value_filename <- paste0(working_folder, "/pos_value_table.csv")
neg_value_filename <- paste0(working_folder, "/neg_value_table.csv")
pos_value_plot_filename <- paste0(working_folder, "/Positive Value Plot.pdf")
neg_value_plot_filename <- paste0(working_folder, "/Negative Value Plot.pdf")
range_plot_filename <- paste0(working_folder, "/Range Plot.pdf")

# read in channel names
panel <- read.csv(panel_location)
panel_keep <- subset(panel, keep == 1)
panel_keep$image_number <- 1:dim(panel_keep)[1]

# Channels needed for identification 
channels_needed <- panel_keep$name
# (This can be specified to ignore other channels, otherwise all keep channels are used).
#channels_needed <- c("CD3", "CD4", "CD8", "CD19", "CD25", "FOXP3", "CD34", "CD31", "CD45")
panel_needed <- panel_keep[panel_keep$name %in% channels_needed, c("image_number", "name")] 

sink(file = channel_list_filename)
cat(paste("All channel names in order are:"))
print(panel_keep)
cat(paste("Channels needed for cell identification are:"))
print(panel_needed)
sink(file = NULL)

# Get names of images
img_filenames <- list.files(image_location, pattern = "*.tif", full.names = T)
img_names <- str_before_last_dot(list.files(image_location, pattern = "*.tif", full.names = F))

# Get number of rows in panel.csv
nRows <- dim(panel)[1]
nKeep <- dim(panel_keep)[1]



# estimate the negative and positive values for every channel (in panel_keep) of every image
# write to csv file to check/store, in Excel
b <- rep(" ", dim(panel_keep)[1])
pos_table <- data.frame(panel_keep$name, 
                        panel_keep$image_number, 
                        matrix(nrow = dim(panel_keep)[1], ncol = length(img_names)))
names(pos_table) <- c("Channel", "No.", img_names)
neg_table <- pos_table

pb = txtProgressBar(min = 0, max = length(img_filenames), initial = 0)
for(i in 1:length(img_filenames)){

  images <- loadImages(img_filenames[i])   # will be a list of one image
  
  # Check number of channels in image
  nChannels <- dim(images[1]@listData[[1]])[3]
  
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
    images[i]@listData[[1]] <- images[i]@listData[[1]][,,which(panel$keep==1)]

  }
  # Otherwise the number of image channels must be the keep channels
  
  channelNames(images) <- panel_keep$name
  
  image_name <- names(images)[1]
  setTxtProgressBar(pb,i)
  
  for(channel in channelNames(images)){
    
    vals <- estimateNegPosValue(images[1], channel, sigma = sigma)
    
    neg_table[which(neg_table$Channel == channel), image_name] = vals[1] 
    pos_table[which(pos_table$Channel == channel), image_name] = vals[2] 
    
  }
  
  rm(images)
}
close(pb)


if(length(img_filenames) > 1){
  # Calculate global levels from all the images
  # Mean level from neg values
  neg_table$global <- rowMeans(neg_table[,3:dim(neg_table)[2]], na.rm = TRUE)
  # For pos table, try to get lower of the values where staining is good
  # ie those above the mid point (t)
  # range will be invariant to the proportion of images with good staining
  # using min should temper any big outliers
  p <- pos_table[,3:dim(pos_table)[2]]
  #mn <- apply(p, 1, min, na.rm = T)
  #mx <- apply(p, 1, max, na.rm = T)
  #t <- (mn + mx)/2
  #pos_table$global <- mapply(function(x, y){min(x[x>y], na.rm = T)}, x=as.data.frame(t(p)), y=t) 
  
  # OR assume all pos values are good (since fg is strict) and take the mean or something
  pos_table$global <- rowMeans(pos_table[,3:dim(pos_table)[2]], na.rm = TRUE)
  #pos_table$global <- apply(p, 1, quantile, probs = 0.1, na.rm = TRUE)  # This was default March 2023-Sept 2023, and often caused pos<neg values.
  
  rm(p)
  #rm(mn, mx, t)
  
} else {  # only 1 image
  neg_table$global <- neg_table[,3]
  pos_table$global <- pos_table[,3]
}

plotPosValues(pos_table, pos_value_plot_filename)
plotPosValues(neg_table, neg_value_plot_filename)
plotRanges(pos_table, neg_table, range_plot_filename)

write.csv(pos_table, file = pos_value_filename, row.names = F)
write.csv(neg_table, file = neg_value_filename, row.names = F)


# Save everything so far
save.image(file = global_data_filename)

