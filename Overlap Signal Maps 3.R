# Overlap Signal Maps: Script 3
# P R Barber, Jan 2023
# This file is specific to the panel and cell type - MODIFY FOR YOUR APPLICATION

# For each cell type that you define in the cell_type_matrix.csv,
# Calculate the probability map and save as tif in celltype_tif folder
# Plot bar charts from data about how abundant each cell type is in each image
# output: Cell Total Plots.pdf

# Cell type matrix has a column for each cell type and 
# a row for each marker. A one indicates positive channel required,
# a -1 indicates channel should be negative. Empty or NA indicates 
# we do not care about that channel.

library(ggplot2)
library(scales)

# User check of working directory.
print("Working in:")
print(getwd())
print("Enter 'y' to proceed:")
proceed = readLines(n=1)
stopifnot(proceed == "y")

# Read previous session
#load("Overlap Signal Maps.RData")

# Check for cell_type_matrix
# If not there create a framework and prompt user to fill it in
matrix_filename <- "cell_type_matrix.csv"
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
# the #000000 black will be used for unclassified image, remove it when using it for bar plots
# Useful website: https://colorbrewer2.org/#type=qualitative&scheme=Set3&n=12
#cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
#               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbPalette <- c("#000000", '#8dd3c7','#ffffb3','#bebada',
               '#fb8072','#80b1d3','#fdb462','#b3de69',
               '#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f')

ct_matrix <- read.csv(matrix_filename, row.names = 1)
# make sure no spaces in the marker row names
row.names(ct_matrix) <- gsub(" ", "", row.names(ct_matrix))

if(dim(ct_matrix)[2] > (length(cbPalette)-1)){
  stop("Not enough colours provided for the cell types requested.
       Add entries to cbPalette in the script.")
}

folder <- "celltype_tif"
dir.create(folder, showWarnings = F)

# How to process the cell type maps
process_os <- function(os, image_name, ct_name){
  
  filename <- paste0(folder, "/", image_name, "_", ct_name, ".tif") 
  writeImage(os, filename)
  
  pixels <- length(os)
  
  total <- sum(os)
  density <- mean(os)
  
  area_over_0.5 <- sum(os>0.5)/pixels # threshold
  os <- ifelse(os > 0.5, os, 0)       # replace <0.5 with zero
  
  total_over_0.5 <- sum(os)
  density_over_0.5 <- mean(os)
  
  c(pixels, total, density, total_over_0.5, density_over_0.5, area_over_0.5)
}

# Storage variable arrays
ct_names <- vector()
image_names <- vector()
pixels <- vector()
total <- vector()
density <- vector()
total_over_0.5 <- vector()
density_over_0.5 <- vector()
area_over_0.5 <- vector()
max_prob_area <- vector()
max_prob_area_perc <- vector()

# Create some names in the environment for the EBImage stacks of probability maps
for(i in 1:length(images)){
  assign(paste0(names(images)[i], "_ct"), NULL)
}

addToGlobalArrays <- function(ct_name, image_name, scores){
  ct_names <<- c(ct_names, ct_name)
  image_names <<- c(image_names, image_name)
  pixels <<- c(pixels, scores[1])
  total <<- c(total, scores[2])
  density <<- c(density, scores[3])
  total_over_0.5 <<- c(total_over_0.5, scores[4])
  density_over_0.5 <<- c(density_over_0.5, scores[5])
  area_over_0.5 <<- c(area_over_0.5, scores[6])
}


############# Define cell types in this loop ######################

pb = txtProgressBar(min = 0, max = length(images), initial = 0)
for(i in 1:length(images)){
  
  image_name <- names(images)[i]
  l <- get(image_name)
  setTxtProgressBar(pb,i)
  
  for(j in 1:dim(ct_matrix)[2]){  # cell types
    
    ct_name <- names(ct_matrix)[j]
      
    # start with an image of ones
    os <- l[,,1]  # copy existing image
    os = os/os    # divide it by itself
      
    for(k in 1:dim(ct_matrix)[1]){  # markers
      
      marker_name <- row.names(ct_matrix)[k]

      v <- ct_matrix[k, j]   # [marker row, cell type col]
      
      # skip if channel not needed
      if(is.na(v)) next
      if(v==0) next
      
      i_p <- l[,,which(channels_needed == marker_name)]
      
      # invert if channel is to be -ve
      if(v < 0) i_p <- (1 - i_p)

      os <- os * i_p
    }
    
    scores <- process_os(os, image_name, ct_name)

    addToGlobalArrays(ct_name, image_name, scores)
    
    ct <- get(paste0(names(images)[i], "_ct"))
    ct <- combine(ct, os)
    assign(paste0(names(images)[i], "_ct"), ct)
    
  }
}
close(pb)

# set dim names for collections of cell type maps
for(i in 1:length(images)){
  ct <- get(paste0(names(images)[i], "_ct"))
  dimnames(ct)[[3]] <- names(ct_matrix)
  assign(paste0(names(images)[i], "_ct"), ct)
}

# Make images of most likely cell type per pixel
pb = txtProgressBar(min = 0, max = length(images), initial = 0)
for(i in 1:length(images)){
  ct <- get(paste0(names(images)[i], "_ct"))
  image_name <- names(images)[i]
  setTxtProgressBar(pb,i)
  
  which_ct = apply(ct, c(1,2), which.max)  # which kind of cell is the max
  max_ct = apply(ct, c(1,2), max)   # what is it's max value
  
  # only allow high probability
  mask <- max_ct > 0.5
  which_ct <- which_ct * mask
  
  areas <- as.vector(table(which_ct))  # ERROR HERE IF NOT ALL CELL TYPES REPRESENTED
  # NEED TO FILL IN MISSING CELL TYPES WITH ZEROS
  
  total_area <- dim(which_ct)[1] * dim(which_ct)[2]
  max_prob_area <- c(max_prob_area, areas[-1])   # exclude first area, the bg
  max_prob_area_perc <- c(max_prob_area_perc, areas[-1]/total_area*100)
  
  # carefully assign 1 colour value per integer using colormap
  y <- colormap(which_ct/dim(ct_matrix)[2], cbPalette[1:(dim(ct_matrix)[2]+1)])
  #display(y)
  
  filename <- paste0(folder, "/", image_name, "_CellMap.tif") 
  writeImage(y, filename)
}
close(pb)

#data <- data.frame(image_names, ct_names, 
#                   total, density, 
#                   total_over_0.5, density_over_0.5, area_over_0.5,
#                   max_prob_area, max_prob_area_perc)

#names(data) <- c("Image", "CellType", 
#                 "Total", "Density", 
#                 "Total_over_0.5", "Density_over_0.5", "Area_over_0.5",
#                "Max_probability_area", "Max_probability_area_percentage")

data <- data.frame(image_names, ct_names, 
                   total, density, 
                   total_over_0.5, density_over_0.5, area_over_0.5)

names(data) <- c("Image", "CellType", 
                 "Total", "Density", 
                 "Total_over_0.5", "Density_over_0.5", "Area_over_0.5")

write.csv(data, file = "CellTypeTotals.csv", row.names = F)

# lock in a cell type order
data$CellType <- factor(data$CellType, levels = names(ct_matrix))


# Plots of cell content

#d <- subset(data, CellType != "CD4Tcell")
d <- data

# FIX THIS TOO - IT DOES NOT REALLY WORK FOR NOV2022 DATA
# if image names are long, make sure there are space to break it up
if(max(nchar(d$Image > 10))){
  d$Image <- gsub("_", " ", d$Image)
}

pdf("Cell Total Plots.pdf")

print(ggplot(d, aes(x = Image, y = Total, fill = CellType)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values=cbPalette[-1]) +
        scale_x_discrete(labels = label_wrap(10)) +
        coord_flip())

print(ggplot(d, aes(x = Image, y = Density, fill = CellType)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values=cbPalette[-1]) +
        scale_x_discrete(labels = label_wrap(10)) +
        coord_flip())

print(ggplot(d, aes(x = Image, y = Total_over_0.5, fill = CellType)) +
          geom_bar(stat = "identity") +
        scale_fill_manual(values=cbPalette[-1]) +
        scale_x_discrete(labels = label_wrap(10)) +
        coord_flip())
  
print(ggplot(d, aes(x = Image, y = Density_over_0.5, fill = CellType)) +
          geom_bar(stat = "identity") +
        scale_fill_manual(values=cbPalette[-1]) +
        scale_x_discrete(labels = label_wrap(10)) +
        coord_flip())

print(ggplot(d, aes(x = Image, y = Area_over_0.5, fill = CellType)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values=cbPalette[-1]) +
        scale_x_discrete(labels = label_wrap(10)) +
        coord_flip())

#print(ggplot(d, aes(x = Image, y = Max_probability_area, fill = CellType)) +
#        geom_bar(stat = "identity") +
#        scale_fill_manual(values=cbPalette[-1]) +
#        scale_x_discrete(labels = label_wrap(10)) +
#        coord_flip())

# density is the same as total with percentage

print(ggplot(d, aes(x = Image, y = Density, fill = CellType)) +
        geom_bar(position = "fill", stat = "identity") +
        scale_y_continuous(labels = scales::percent) +
        scale_fill_manual(values=cbPalette[-1]) +
        scale_x_discrete(labels = label_wrap(10)) +
        coord_flip())

print(ggplot(d, aes(x = Image, y = Density_over_0.5, fill = CellType)) +
        geom_bar(position = "fill", stat = "identity") +
        scale_y_continuous(labels = scales::percent) +
        scale_fill_manual(values=cbPalette[-1]) +
        scale_x_discrete(labels = label_wrap(10)) +
        coord_flip())

print(ggplot(d, aes(x = Image, y = Area_over_0.5, fill = CellType)) +
        geom_bar(position = "fill", stat = "identity") +
        scale_y_continuous(labels = scales::percent) +
        scale_fill_manual(values=cbPalette[-1]) +
        scale_x_discrete(labels = label_wrap(10)) +
        coord_flip())

#print(ggplot(d, aes(x = Image, y = Max_probability_area, fill = CellType)) +
#        geom_bar(position = "fill", stat = "identity") +
#        scale_y_continuous(labels = scales::percent) +
#        scale_fill_manual(values=cbPalette[-1]) +
#        scale_x_discrete(labels = label_wrap(10)) +
#        coord_flip())

dev.off()

