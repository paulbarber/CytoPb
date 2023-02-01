# Overlap Signal Maps: Script 3
# P R Barber, Jan 2023
# This file is specific to the panel and cell type - MODIFY FOR YOUR APPLICATION

# For each cell type that you define in the loop below,
# Calculate the probability map and save as tif in celltype_tif folder
# Plot bar charts from data about how abundant each cell type is i neach image
# output: Cell Total Plots.pdf

# Read previous session
#load("Overlap Signal Maps.RData")

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

pb = txtProgressBar(min = 0, max = length(channels_needed), initial = 0)
for(i in 1:length(images)){
  
  image_name <- names(images)[i]
  l <- get(image_name)
  setTxtProgressBar(pb,i)
  
  #####################################################
  # CD4 Tcells (CD3+ CD4+)   
  # an extra pop in case of no T regs or effector
  ct_name <- "CD4Tcell"
  
  i_p1 <- l[,,which(channels_needed == "CD3")]
  i_p2 <- l[,,which(channels_needed == "CD4")]
  
  os <-i_p1 * i_p2
  
  # This bit same for all cell types
  scores <- process_os(os, image_name, ct_name)
  addToGlobalArrays(ct_name, image_name, scores)
  
  #####################################################
  # CD8 Tcells (CD3+ CD8+)
  ct_name <- "CD8Tcell"
  
  i_p1 <- l[,,which(channels_needed == "CD3")]
  i_p2 <- l[,,which(channels_needed == "CD8")]
  
  os <-i_p1 * i_p2
  
  # This bit same for all cell types
  scores <- process_os(os, image_name, ct_name)
  addToGlobalArrays(ct_name, image_name, scores)
  
  #####################################################
  # T Effector cells (CD3+ CD4+ CD25-)
  ct_name <- "TEffcell"
  
  i_p1 <- l[,,which(channels_needed == "CD3")]
  i_p2 <- l[,,which(channels_needed == "CD4")]
  i_p3 <- l[,,which(channels_needed == "CD25")]
  
  os <-i_p1 * i_p2 * (1-i_p3)
  
  # This bit same for all cell types
  scores <- process_os(os, image_name, ct_name)
  addToGlobalArrays(ct_name, image_name, scores)

  #####################################################
  # T Regulatory cells (CD3+ CD4+ CD25+ FOXP3+)
  ct_name <- "TRegcell"
  
  i_p1 <- l[,,which(channels_needed == "CD3")]
  i_p2 <- l[,,which(channels_needed == "CD4")]
  i_p3 <- l[,,which(channels_needed == "CD25")]
  i_p4 <- l[,,which(channels_needed == "FOXP3")]
  
  os <-i_p1 * i_p2 * i_p3 * i_p4
  
  # This bit same for all cell types
  scores <- process_os(os, image_name, ct_name)
  addToGlobalArrays(ct_name, image_name, scores)
  
  #####################################################
  # B cells (CD3- CD19+)
  ct_name <- "Bcell"
  
  i_p1 <- l[,,which(channels_needed == "CD3")]
  i_p2 <- l[,,which(channels_needed == "CD19")]
  
  os <-(1-i_p1) * i_p2
  
  # This bit same for all cell types
  scores <- process_os(os, image_name, ct_name)
  addToGlobalArrays(ct_name, image_name, scores)
  
  #####################################################
  # Endothelium (CD34+ CD31+)
  ct_name <- "Endothelium"
  
  i_p1 <- l[,,which(channels_needed == "CD34")]
  i_p2 <- l[,,which(channels_needed == "CD31")]
  
  os <-i_p1 * i_p2
  
  # This bit same for all cell types
  scores <- process_os(os, image_name, ct_name)
  addToGlobalArrays(ct_name, image_name, scores)
  
  #####################################################
  # Progenitor (CD34+ CD31-)
  ct_name <- "Progenitor"
  
  i_p1 <- l[,,which(channels_needed == "CD34")]
  i_p2 <- l[,,which(channels_needed == "CD31")]
  
  os <-i_p1 * (1-i_p2)
  
  # This bit same for all cell types
  scores <- process_os(os, image_name, ct_name)
  addToGlobalArrays(ct_name, image_name, scores)
  
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
