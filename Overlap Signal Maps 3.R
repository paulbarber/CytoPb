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
# we do not care about htat channel.

# User check of working directory.
print("Working in:")
print(getwd())
print("Enter 'y' to proceed:")
proceed = readLines(n=1)
stopifnot(proceed == "y")

# Read previous session
#load("Overlap Signal Maps.RData")

ct_matrix <- read.csv("cell_type_matrix.csv", row.names = 1)
# make sure no spaces in the marker row names
row.names(ct_matrix) <- gsub(" ", "", row.names(ct_matrix))

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

  }
}

close(pb)

data <- data.frame(image_names, ct_names, 
                   total, density, 
                   total_over_0.5, density_over_0.5, area_over_0.5)

names(data) <- c("Image", "CellType", 
                 "Total", "Density", 
                 "Total_over_0.5", "Density_over_0.5", "Area_over_0.5")

write.csv(data, file = "CellTypeTotals.csv", row.names = F)



# Plots of cell content

#d <- subset(data, CellType != "CD4Tcell")
d <- data

pdf("Cell Total Plots.pdf")

print(ggplot(d, aes(x = Image, y = Total, fill = CellType)) +
        geom_bar(stat = "identity") +
        coord_flip())

print(ggplot(d, aes(x = Image, y = Density, fill = CellType)) +
        geom_bar(stat = "identity") +
        coord_flip())

print(ggplot(d, aes(x = Image, y = Total_over_0.5, fill = CellType)) +
          geom_bar(stat = "identity") +
          coord_flip())
  
print(ggplot(d, aes(x = Image, y = Density_over_0.5, fill = CellType)) +
          geom_bar(stat = "identity") +
          coord_flip())

print(ggplot(d, aes(x = Image, y = Area_over_0.5, fill = CellType)) +
        geom_bar(stat = "identity") +
        coord_flip())

# density is the same as total with percentage

print(ggplot(d, aes(x = Image, y = Density, fill = CellType)) +
        geom_bar(position = "fill", stat = "identity") +
        scale_y_continuous(labels = scales::percent) +
        coord_flip())

print(ggplot(d, aes(x = Image, y = Density_over_0.5, fill = CellType)) +
        geom_bar(position = "fill", stat = "identity") +
        scale_y_continuous(labels = scales::percent) +
        coord_flip())

dev.off()
