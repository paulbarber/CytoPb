# CytoPb: Script 2 Make Channel Probabilities
# P R Barber, Jan 2023

# Following from script 1 (using same working directory).
# Plot the image positive values as a SNR indicator
# output: Positive Value Plot.pdf
# Normalise each needed channel to neg and pos values.
# Calculate marker probability map.
# output: scaled and map images to channel_png
# Save the complete workspace for follow on scripts.

library(ggplot2)



if(!exists("working_folder")){
  working_folder <- choose.dir(caption = "Select data folder")
}

print("CytoPb 2 Working in:")
print(working_folder)

global_data_filename <- paste0(working_folder, "/CytoPb.RData")

# Read previous session
load(global_data_filename)

# folder to save channel QC images to
channel_png_folder <- paste0(working_folder, "/channel_png/")
dir.create(channel_png_folder, showWarnings = F)

# folder for R objects
objects_folder <- paste0(working_folder, "/objects/")
dir.create(objects_folder, showWarnings = F)

# Read in pos and neg values, user may have tweaked them from last script
pos_table <- read.csv(pos_value_filename)
neg_table <- read.csv(neg_value_filename)

# Are we going to use the global pos and neg values, or individuals for each image?
if(!exists("use_global_ranges")){
  use_global = FALSE
}else{
  use_global = use_global_ranges
  rm(use_global_ranges)   # this stops in going into global_data_filename and being overwritten if changed and rerun
}

# Blue to red pallette
jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

# Channel images and density plots
pb = txtProgressBar(min = 0, max = length(img_filenames)*length(channels_needed), initial = 0)
for(i in 1:length(img_filenames)){
  
  images <- loadImages(img_filenames[i])   # will be a list of one image
  image_name <- names(images)[1]
  
  # set dim names for the images
  dimnames(images@listData[[1]])[[3]] <- channels_needed
  
  # somewhere to store the maps for this image
  channel_probability_maps = NULL

  for(j in 1:length(channels_needed)){
    
    setTxtProgressBar(pb,(i-1)*length(channels_needed) + j)
    
    # Channel of interest
    channel <- channels_needed[j]
    i_p1 <- images@listData[[1]][,,channel]
  
    # Get scaling parameters
    if(use_global){
      nv1 <- neg_table[which(neg_table$Channel == channel), "global"]     
      pv1 <- pos_table[which(pos_table$Channel == channel), "global"]     
    }else{
      table_name <- make.names(image_name)   # to match table col
      nv1 <- neg_table[which(neg_table$Channel == channel), table_name]
      pv1 <- pos_table[which(pos_table$Channel == channel), table_name]
      
    }
    
    # do some checks against the global positive value
    g <- pos_table[which(pos_table$Channel == channel), "global"]
    
    # Check these values, if nv is NA set to 0, if pv is NA set to global value
    if(is.na(nv1)) nv1 = 0
    if(is.na(pv1)) {
      pv1 = g
    }
    if(is.na(pv1)) {   # final check
      pv1 = 100000
    }
    
    # Check for fg>bg
    if(pv1<nv1) {
      pv1 = g
      if(pv1<nv1) {
        nv1 = 0
        if(pv1<nv1) {
          pv1 = 100000
        }
      }
    }
    
    # if pos is << global value, channel or fg segmentation must have failed
    if(pv1 < g/10) {
      pv1 = g
    }
    
    
    # Write a png of the channel for convenience
    # I cannot get EBImage to write a nice image out! This is the best I can do.
    # It uses png::writePNG
    filename <- paste0(channel_png_folder, image_name, "_", channel, ".png")
    writeImage(i_p1*256, filename)
    
    # blur to account for cell size/position uncertainties
    i_p1 <- gblur(i_p1, sigma = sigma)
    
    # re-scale 0-1 to get an estimate of the probability of being positive.
    i_p1 <- (i_p1 - nv1)/(pv1 - nv1)
    
    # Threshold/clamp to make sure strictly 0-1
    #i_p1[i_p1<0] <- 0
    #i_p1[i_p1>1] <- 1
    
    # Use Sigmoid instead of clamping, bigger f is steeper slope
    # Will ensure 0-1 without sharp cutoff
    f = 5
    i_p1 <- 1/(1+exp(-f*(i_p1-0.5)))
    
    filename <- paste0(channel_png_folder, image_name, "_", channel, "_map.png")
    y = colormap(i_p1, jet.colors(256))
    writeImage(y, filename)
    
    # Store these probability maps for later use, in order as channels_needed
    channel_probability_maps <- combine(channel_probability_maps, i_p1)

  }
  
  rm(images)
  rm(y)
  rm(i_p1)
  
  # set dim names for collections of marker maps
  dimnames(channel_probability_maps)[[3]] <- channels_needed
  
  filename <- paste0(objects_folder, 
                     image_name, 
                     "_channel_probability_maps.RData")
  
  save(channel_probability_maps, file = filename)
  rm(channel_probability_maps)
}
close(pb)
    
# Save everything so far
save.image(file = global_data_filename)



