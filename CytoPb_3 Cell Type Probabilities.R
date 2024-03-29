# CytoPb: Script 3 Cell Type Probabilities
# P R Barber, Jan 2023
# This file is specific to the panel and cell type - MODIFY FOR YOUR APPLICATION

# For each cell type that you define in the cell_type_matrix.csv,
# Calculate the probability map using fuzzy logic multiplication or
# Zadeh operators (AND = min, OR = max).
# and save as tif in celltype_maps folder.
# Plot bar charts from data about how abundant each cell type is in each image
# output: Cell Total Plots.pdf

# Cell type matrix has a column for each cell type and 
# a row for each marker. A one indicates positive channel required,
# a -1 indicates channel should be negative. Empty or NA indicates 
# we do not care about that channel.

suppressMessages(library(cytomapper))
suppressMessages(library(ggplot2))

if(!exists("working_folder")){
  working_folder <- choose.dir(caption = "Select data folder")
}

cat(paste("CytoPb 3 Working in: ", working_folder, "\n"))

global_data_filename <- paste0(working_folder, "/CytoPb.RData")

# Read previous session
load(global_data_filename)


# Output File Locations
# Check for existing data.
results_folder <- paste0(working_folder, "/CellType_Results")
r <- dir.create(results_folder, showWarnings = F)
if(!r){   # Folder already exists, data may be overwritten
  cat("WARNING: CellType_Results folder already exisits.\nData will be overwritten. You could rename the previous results folder.\n")
  cat("Enter 'y' to proceed:\n")
  proceed = readLines(n=1)
  stopifnot(proceed == "y")
}

celltype_objects_folder <- paste0(results_folder, "/objects/")
dir.create(celltype_objects_folder, showWarnings = F)

markerpercellbyimage_filename <- paste0(results_folder, "/Marker per CellType by image.pdf")
markerpercell_filename <- paste0(results_folder, "/Marker per CellType.pdf")
markerpercelltable_filename <- paste0(results_folder, "/Marker per CellType.csv")
CellTypeTotals_filename <- paste0(results_folder, "/Cell Type Totals.csv")
CellTotalPlotsfilename <- paste0(results_folder, "/Cell Total Plots.pdf")

# folder to save images to
celltype_png_folder <- paste0(results_folder, "/celltype_png/")
dir.create(celltype_png_folder, showWarnings = F)
celltype_map_folder <- paste0(results_folder, "/celltype_maps/")
dir.create(celltype_map_folder, showWarnings = F)

# Input File Locations
matrix_filename <- paste0(working_folder, "/cell_type_matrix.csv")
copy_matrix_filename <- paste0(results_folder, "/cell_type_matrix.csv")
cell_type_colours_filename <- paste0(working_folder, "/cell_type_colours.txt")


# Set the high probability threshold, but allow it to be predefined differently
if(!exists("probability_threshold")){
  prob_threshold = 0.5
}else{
  prob_threshold = probability_threshold
  rm(probability_threshold)   # this stops in going into global_data_filename and being overwritten if changed and rerun
}


# Check for cell_type_matrix
# If not there create a framework and prompt user to fill it in
if(!file.exists(matrix_filename)){

  ct_matrix <- data.frame(channels_needed)
  names(ct_matrix) <- "Marker" 
  
  # Create some standard cell type columns
  ct_matrix$Mac <- ""
  ct_matrix$Mac[which(channels_needed == "CD68")] <- 1
  
  ct_matrix$Tcell <- ""
  ct_matrix$Tcell[which(channels_needed == "CD68")] <- -1
  ct_matrix$Tcell[which(channels_needed == "CD3")] <- 1
  
  ct_matrix$Bcell <- ""
  ct_matrix$Bcell[which(channels_needed == "CD68")] <- -1
  ct_matrix$Bcell[which(channels_needed == "CD3")] <- -1
  ct_matrix$Bcell[which(channels_needed == "CD20")] <- 1
  
  write.csv(ct_matrix, matrix_filename, row.names = F)
  
  stop("cell_type_matrix.csv file did not exsist.
        A template has been created. 
        Please complete with your cell types of interest.")
}

# Set a colour pallette
# the #000000 black will be used for unclassified image, will be removed when using it for bar plots
# Useful website: https://colorbrewer2.org/#type=qualitative&scheme=Set3&n=12
#cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
#               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

if(file.exists(cell_type_colours_filename)){

  cb <- read.delim(cell_type_colours_filename, header = F)
  cbPalette <- c("black", cb$V1) 
    
}else{
  
  # A nice palette
#  cbPalette <- c("#000000", '#8dd3c7','#ffffb3','#bebada',
#                 '#fb8072','#80b1d3','#fdb462','#b3de69',
#                 '#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f')
  
  # A brighter palette
  cbPalette <- c('#000000', '#a6cee3','#1f78b4','#b2df8a','#33a02c',
                 '#fb9a99','#e31a1c','#fdbf6f','#ff7f00',
                 '#cab2d6','#6a3d9a','#ffff99','#b15928')
  
}

# Check colours are valid since it is annoying to get an error later on
areColors <- function(x) {
  sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)), 
             error = function(e) FALSE)
  })
}

colours_valid <- areColors(cbPalette)
if(!all(colours_valid)){
  cat(colours_valid)
  stop("Some provided colours are not valid R colours.")
}


ct_matrix <- read.csv(matrix_filename, row.names = 1)
# copy into the results folder for a record
file.copy(matrix_filename, copy_matrix_filename)
# make sure no spaces in the marker row names
row.names(ct_matrix) <- gsub(" ", "", row.names(ct_matrix))
# and the list of channels matches
full_channel_list <- gsub(" ", "", channels_needed)

n_cell_types <- dim(ct_matrix)[2]
n_markers <- dim(ct_matrix)[1]

if(n_cell_types > (length(cbPalette)-1)){
  cat("Not enough colours provided for the cell types requested.
       You can add entries to a cell_type_colours.txt file. Meanwhile,
        a rainbow pallette will be used.")
  cbPalette <- c("#000000", rainbow(n_cell_types+2))
}

# How to process the cell type maps
process_os <- function(os, image_name, ct_name){
  
  filename <- paste0(celltype_png_folder, "/", image_name, "_", ct_name, ".png") 
  y = colormap(os, jet.colors(256))
  writeImage(y, filename)
  rm(y)
  
  pixels <- length(os)
  
  total <- sum(os)
  density <- mean(os)
  
  area_highProb <- sum(os > prob_threshold)/pixels # threshold
  os <- ifelse(os > prob_threshold, os, 0)       # replace <prob_threshold with zero
  
  total_highProb <- sum(os)
  density_highProb <- mean(os)
  
  c(pixels, total, density, total_highProb, density_highProb, area_highProb)
}

# Storage variable arrays
ct_names <- vector()
image_names <- vector()
pixels <- vector()
total <- vector()
density <- vector()
total_highProb <- vector()
density_highProb <- vector()
area_highProb <- vector()
max_prob_area <- vector()
max_prob_area_perc <- vector()


addToGlobalArrays <- function(ct_name, image_name, scores){
  ct_names <<- c(ct_names, ct_name)
  image_names <<- c(image_names, image_name)
  pixels <<- c(pixels, scores[1])
  total <<- c(total, scores[2])
  density <<- c(density, scores[3])
  total_highProb <<- c(total_highProb, scores[4])
  density_highProb <<- c(density_highProb, scores[5])
  area_highProb <<- c(area_highProb, scores[6])
}


############# Cell types ######################

pb = txtProgressBar(min = 0, max = length(img_names), initial = 0)
for(i in 1:length(img_names)){
  
  image_name <- img_names[i]
  
  # somewhere to store the maps for this image
  celltype_probability_maps = NULL
  
  # retrieve channel probability 
  l <- loadChannelMapObject(image_name)

  setTxtProgressBar(pb,i)

  for(j in 1:n_cell_types){  # cell types
    
    ct_name <- names(ct_matrix)[j]
      
## AND = multiplication
    # start with an image of ones
##    os <- l[,,1]  # copy existing image
##    os = os/os    # divide it by itself
    os = NULL
      
    for(k in 1:n_markers){  # markers
      
      marker_name <- row.names(ct_matrix)[k]

      v <- ct_matrix[k, j]   # [marker row, cell type col]
      
      # skip if channel not needed
      if(is.na(v)) next
      if(v==0) next
      
      i_p <- l[,,which(full_channel_list == marker_name)]
      
      # invert if channel is to be -ve
      if(v < 0) i_p <- (1 - i_p)

      ## AND = multiplication
      ##      os <- os * i_p
      
      ### AND = min
      ### combine all channels into a stack and do min outside the loop
      os <- EBImage::combine(os, i_p)
      
    }
    
    ### AND = min
    os <- apply(os, c(1,2), min)
    
    scores <- process_os(os, image_name, ct_name)

    addToGlobalArrays(ct_name, image_name, scores)
    
    # Store these probability maps for later use
    celltype_probability_maps <- EBImage::combine(celltype_probability_maps, os)
    
  }
  
  ############ Total tissue area by OR of all channels ##################
  # OR by inverting all channels, AND them all, invert the result
  # or OR = max
  
  ct_name <- "Total_Tissue"
  
  ## AND = multiplication
  # start with an image of ones
  ##os <- l[,,1]  # copy existing image
  ##os = os/os    # divide it by itself
  os = NULL
  
  for(k in 1:n_markers){  # markers
    
    marker_name <- row.names(ct_matrix)[k]
    
    i_p <- l[,,which(full_channel_list == marker_name)]
    
    ## AND = multiplication
    # invert channel
    ##i_p <- (1 - i_p)
    
    ##os <- os * i_p
    
    ### OR = max
    ### combine all channels into a stack and do max outside the loop
    os <- EBImage::combine(os, i_p)
  }
  
  ## AND = multiplication
  # invert result
  ##os <- (1 - os)
  
  ### OR = max
  os <- apply(os, c(1,2), max)
  
  
  scores <- process_os(os, image_name, ct_name)
  
  addToGlobalArrays(ct_name, image_name, scores)
  
  # DON'T Store this probability map for later use
  # Only generate stats on it
  #celltype_probability_maps <- EBImage::combine(celltype_probability_maps, os)
  
  ########################################################################
  
  
  saveCellTypeMapObject(celltype_probability_maps, image_name, names(ct_matrix))
  rm(celltype_probability_maps)
  
}
close(pb)

rm(os)
rm(i_p)
rm(l)


# Make images of most likely cell type per pixel
mean_per_ct <- matrix(nrow = n_cell_types, ncol = length(full_channel_list), data = 0)
mean_per_ct_table <- data.frame()
pdf(markerpercellbyimage_filename)
pb = txtProgressBar(min = 0, max = length(img_filenames), initial = 0)
for(i in 1:length(img_filenames)){
  
  # load the image we need
  images <- loadImages(img_filenames[i])   # will be a list of one image
  image_name <- names(images)[1]
  
  # set dim names for the images
  dimnames(images@listData[[1]])[[3]] <- full_channel_list
  
  # load the cell type prob map we need
  ct <- loadCellTypeMapObject(image_name)

  setTxtProgressBar(pb,i)
  
  which_ct = apply(ct, c(1,2), which.max)  # which kind of cell is the max
  max_ct = apply(ct, c(1,2), max)   # what is it's max value
  
  # only allow high probability
  mask <- max_ct > prob_threshold
  which_ct <- which_ct * mask
  
  areas <- as.vector(table(factor(which_ct, levels = 0:n_cell_types))) # factor makes sure all cell types are represented
  areas <- c(areas, NA) # add extra element for Total_Tissue cell type
  
  total_area <- dim(which_ct)[1] * dim(which_ct)[2]
  max_prob_area <- c(max_prob_area, areas[-1])   # exclude first area, the bg
  max_prob_area_perc <- c(max_prob_area_perc, areas[-1]/total_area*100)
  
  # carefully assign 1 colour value per integer using colormap
  y <- colormap(which_ct/n_cell_types, cbPalette[1:(n_cell_types+1)])
  #display(y)
  
  filename <- paste0(celltype_map_folder, "/", image_name, "_CellMap.png") 
  writeImage(y, filename)
  rm(y)
  
  # find average marker strength per cell type
  strength_per_ct <- matrix(nrow = n_cell_types, ncol = length(full_channel_list))
  colnames(strength_per_ct) <- full_channel_list
  rownames(strength_per_ct) <- names(ct_matrix)
  img <- images@listData[[1]]   # all channels for this image
  for(j in 1:n_cell_types){
    # mask for this cell type
    mask <- which_ct == j
    # mean for each channel
    strength <- vector()
    for(k in 1:length(full_channel_list)){
      ch <- img[,,k]
      
      s <- mean(ch[mask], na.rm = TRUE)
      if(is.nan(s)){s = 0}
      strength <- c(strength, s)
    }
    strength_per_ct[j,] <- strength
  }
  
  mean_per_ct <- strength_per_ct + mean_per_ct
  
  spc <- as.data.frame(strength_per_ct)
  spc$Image <- image_name
  spc$cellType <- names(ct_matrix)
  mean_per_ct_table <- rbind(mean_per_ct_table, spc)
  
  # heatmap plot of expressions versus cell type per image
  m <- scale(strength_per_ct)
  m <- as.data.frame(m)
  m$CellType <- rownames(strength_per_ct)
  d <- tidyr::gather(m, Channel, Mean, 1:length(full_channel_list), factor_key=TRUE)
  d$Channel <- factor(d$Channel, levels=unique(d$Channel)) # keep channel order in plot
  d$CellType <- factor(d$CellType, levels=unique(d$CellType)) # keep CellType order in plot
  
  print(ggplot(d, aes(CellType, Channel, fill = Mean)) + 
          geom_tile() + 
          guides(fill=guide_colourbar(title = "mean")) +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10),
                legend.position = "right") +
          ggtitle(image_name)) 
  
}
dev.off()
close(pb)

rm(images)
rm(ct)
rm(ch)
rm(img)
rm(which_ct)
rm(mask)
rm(max_ct)

# Finish the calculation of mean marker strength per cell type
mean_per_ct <- mean_per_ct / length(img_filenames)
mean_per_ct <- as.data.frame(mean_per_ct)

# Write full table of marker per cell type
write.csv(mean_per_ct_table, file = markerpercelltable_filename, row.names = F)

cell_type_data <- data.frame(image_names, ct_names, 
                   total, density, 
                   total_highProb, density_highProb, area_highProb,
                   max_prob_area, max_prob_area_perc)

names(cell_type_data) <- c("Image", "CellType", 
                 "Total", "Density", 
                 "Total_highProb", "Density_highProb", "Area_highProb",
                "Max_probability_area", "Max_probability_area_percentage")

rm(image_names, ct_names, 
   total, density, 
   total_highProb, density_highProb, area_highProb,
   max_prob_area, max_prob_area_perc)

write.csv(cell_type_data, file = CellTypeTotals_filename, row.names = F)

# lock in a cell type order
cell_type_data$CellType <- factor(cell_type_data$CellType, levels = names(ct_matrix))


# Plots of cell content

d <- subset(cell_type_data, CellType != "Total_Tissue")

# reduce text size when number of images is large
n_images <- length(img_names)
rel_size = 1
if(n_images > 50){
  rel_size = 50/n_images
}

pdf(CellTotalPlotsfilename)

# print(ggplot(d, aes(x = Image, y = Total, fill = CellType)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values=cbPalette[-1]) +
#         coord_flip() +
#         theme(axis.text.y = element_text(size = rel(rel_size))))

 print(ggplot(d, aes(x = Image, y = Density, fill = CellType)) +
         geom_bar(stat = "identity") +
         scale_fill_manual(values=cbPalette[-1]) +
         coord_flip() +
         theme(axis.text.y = element_text(size = rel(rel_size))))
 
# print(ggplot(d, aes(x = Image, y = Total_highProb, fill = CellType)) +
#           geom_bar(stat = "identity") +
#         scale_fill_manual(values=cbPalette[-1]) +
#         coord_flip() +
#         theme(axis.text.y = element_text(size = rel(rel_size))))
 
print(ggplot(d, aes(x = Image, y = Density_highProb, fill = CellType)) +
          geom_bar(stat = "identity") +
        scale_fill_manual(values=cbPalette[-1]) +
        coord_flip() +
        theme(axis.text.y = element_text(size = rel(rel_size))))

# print(ggplot(d, aes(x = Image, y = Area_highProb, fill = CellType)) +
#         geom_bar(stat = "identity") +
#         scale_fill_manual(values=cbPalette[-1]) +
#         coord_flip() +
#         theme(axis.text.y = element_text(size = rel(rel_size))))

print(ggplot(d, aes(x = Image, y = Max_probability_area, fill = CellType)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values=cbPalette[-1]) +
        coord_flip() +
        theme(axis.text.y = element_text(size = rel(rel_size))))

# density is the same as total with percentage

 print(ggplot(d, aes(x = Image, y = Density, fill = CellType)) +
         geom_bar(position = "fill", stat = "identity") +
         scale_y_continuous(labels = scales::percent) +
         scale_fill_manual(values=cbPalette[-1]) +
         coord_flip() +
         theme(axis.text.y = element_text(size = rel(rel_size))))

print(ggplot(d, aes(x = Image, y = Density_highProb, fill = CellType)) +
        geom_bar(position = "fill", stat = "identity") +
        scale_y_continuous(labels = scales::percent) +
        scale_fill_manual(values=cbPalette[-1]) +
        coord_flip() +
        theme(axis.text.y = element_text(size = rel(rel_size))))

# print(ggplot(d, aes(x = Image, y = Area_highProb, fill = CellType)) +
#         geom_bar(position = "fill", stat = "identity") +
#         scale_y_continuous(labels = scales::percent) +
#         scale_fill_manual(values=cbPalette[-1]) +
#         coord_flip() +
#         theme(axis.text.y = element_text(size = rel(rel_size))))

print(ggplot(d, aes(x = Image, y = Max_probability_area, fill = CellType)) +
        geom_bar(position = "fill", stat = "identity") +
        scale_y_continuous(labels = scales::percent) +
        scale_fill_manual(values=cbPalette[-1]) +
        coord_flip() +
        theme(axis.text.y = element_text(size = rel(rel_size))))

dev.off()


# heatmap plot of expressions versus cell type
m <- scale(mean_per_ct)
m <- as.data.frame(m)
m$CellType <- rownames(mean_per_ct)
d <- tidyr::gather(m, Channel, Mean, 1:length(full_channel_list), factor_key=TRUE)
d$Channel <- factor(d$Channel, levels=unique(d$Channel)) # keep channel order in plot
d$CellType <- factor(d$CellType, levels=unique(d$CellType)) # keep CellType order in plot

pdf(markerpercell_filename)
print(ggplot(d, aes(CellType, Channel, fill = Mean)) + 
        geom_tile() + 
        guides(fill=guide_colourbar(title = "mean")) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10),
              legend.position = "right"))
dev.off()


# For the neighborhood regions in script 4, make sure these vars are deleted
# so they can be created for the current cell types
# save this data for future runs with the same grid size
rm(list = ls(pattern = "region_data_[0-9]+"))
rm(list = ls(pattern = "regions_per_image_[0-9]+"))
rm(list = ls(pattern = "region_d1_per_image_[0-9]+"))
rm(list = ls(pattern = "region_d2_per_image_[0-9]+"))


# Save everything so far
save.image(file = global_data_filename)
