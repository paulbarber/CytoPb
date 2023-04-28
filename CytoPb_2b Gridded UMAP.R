# CytoPb: Script 2b Gridded UMAP - EXPERIMENTAL/UNDER DEVELOPMENT
# P R Barber, Apr 2023

# Following from script 2 (using same working directory).
# Grid the images up into regions (could use "superpixels")
# Calculate a UMAP plot of regions based on the marker probability maps
# NB umap optimised for clustering
# Cluster the umap
# Plot heatmap of marker map values for the determined clusters
# The aim here is probably to identify cell types, keep the 
# region_grid_size_um quite small.

suppressMessages(library(umap))
suppressMessages(library(ggplot2))
suppressMessages(library(dbscan))
suppressMessages(library(tidyr))

source("R/binImage.R")

# OPTIONS
if(!exists("region_grid_size_um")){   # c("All", "Random", "High.Marker")
  grid_um = 10
}else{
  grid_um = region_grid_size_um
  rm(region_grid_size_um)   # this stops in going into global_data_filename and being overwritten if changed and rerun
}

# Select the images to process, doing all together may be too much
# UMAP takes a while
# hdbscan clustering takes a while and needs lots of memory
#images_to_process <- img_names   # Process all images
images_to_process <- img_names[1]




if(!exists("working_folder")){
  working_folder <- choose.dir(caption = "Select data folder")
}

cat(paste("CytoPb 2b Working in:", working_folder, "\n"))

global_data_filename <- paste0(working_folder, "/CytoPb.RData")

# Read previous session
load(global_data_filename)

grid_size = grid_um/image_scale_umperpixel

# folder to save to
clustering_folder <- paste0(working_folder, "region_clustering/")
dir.create(clustering_folder, showWarnings = F)

# data will be a vector of all regions from all images
data <- NULL
# keep track of how many regions per image
regions_per_image <- NULL
region_d1_per_image <- NULL
region_d2_per_image <- NULL

cat("Collecting region data for all images...\n")
pb = txtProgressBar(min = 0, max = length(images_to_process), initial = 0)
for(i in 1:length(images_to_process)){
  
  setTxtProgressBar(pb, i)
  image_name <- images_to_process[i]
  
  # Get some channel probability data  # retrieve channel probability 
  filename <- paste0(objects_folder, 
                     image_name, 
                     "_channel_probability_maps.RData")
  load(filename)
  image <- channel_probability_maps
  rm(channel_probability_maps)
  
  
  d <- binImage(image, max, bin_size = grid_size)
  d <- d@.Data
  nRegions <- dim(d)[1]*dim(d)[2]

  # Note how many regions for this image
  regions_per_image <- c(regions_per_image, nRegions) 
  region_d1_per_image <- c(region_d1_per_image, dim(d)[1])
  region_d2_per_image <- c(region_d2_per_image, dim(d)[2])
  
  dim(d) <- c(nRegions, dim(d)[3])    # reshape to 2d [region,channel]
  
  data <- rbind(data, d)
  rm(d)
}  
close(pb)


# Copy ch names and make sure all the names are good for the plots below
colnames(data) <- make.names(panel_needed$name)

# Remove DNA channels from this analysis
data <- subset(data, select = -grep("DNA", colnames(data)))

# Set umap settings for pre-clustering, tight clusters
custom.settings = umap.defaults
custom.settings$n_neighbors = 30
custom.settings$min_dist = 0.0001
custom.settings$n_components=2
custom.settings$random_state=42


cat(paste("Performing UMAP with", dim(data)[1], "regions...\n"))
data.umap <- umap(data, config = custom.settings, fast_sgd=T)


u <- as.data.frame(data.umap$layout)
names(u) <- c("UMAP1", "UMAP2")
data2 <- cbind(data, u)

cat("Saving UMAPs per channel...")
pb = txtProgressBar(min = 0, max = length(colnames(data)), initial = 0)
pdf(paste0(clustering_folder, "region_channel_umap_plots.pdf"))
i=1
for(ch in colnames(data)){
  setTxtProgressBar(pb, i)
  print(ggplot(data2, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes_string(colour = ch), size = 1))
  i = i + 1
}
dev.off()
close(pb)
  
  
cat("Clustering...\n")
# minPts is min cluster size BUT also a smoothing factor!
#cl <- hdbscan(u, minPts = 5)
cl <- hdbscan(data, minPts = 5)

cat("Saving plots...\n")

# NB Cluster labels start at ZERO, but 0 are "noise points"!!!!!

pdf(paste0(clustering_folder, "region_clustering.pdf"))
plot(cl, show_flat = T)   # simpler plot showing most stable clusters
dev.off()

# cluster colours
nClusters = max(cl$cluster)
cbPalette <- c("#000000", rainbow(nClusters+2))

data3 <- cbind(data, cl$cluster)
colnames(data3) <- c(colnames(data), "cluster")
data4 <- cbind(data2, cl$cluster)
colnames(data4)[length(colnames(data4))] <- "cluster"
data4$cluster <- factor(data4$cluster, levels = 1:(nClusters))


pdf(paste0(clustering_folder, "region_cluster_umap_plot.pdf"))
print(ggplot(data4, aes(x = UMAP1, y = UMAP2)) +
        geom_point(aes(colour = cluster), size = 1) +
        scale_colour_manual(values=cbPalette[-1]))
dev.off()


# cluster heatmap
data3m <- tidyr::gather(as.data.frame(data3), key = "channel", value = "value", colnames(data), factor_key=TRUE)
pdf(paste0(clustering_folder, "region_cluster_heatmap.pdf"))
print(ggplot(data3m, aes(y = cluster, x = channel, fill = value)) +
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 5),
        legend.position = "none") )
dev.off()


# Region cluster maps
pb = txtProgressBar(min = 0, max = length(images_to_process), initial = 0)
st = 1
for(i in 1:length(images_to_process)){
  
  setTxtProgressBar(pb, i)
  image_name <- images_to_process[i]
  
  # Get the region cluster labels for this image
  d <- data4[st:(st+regions_per_image[i]-1), "cluster"]
  
  # reshape cluster labels like an image
  dim(d) <- c(region_d1_per_image[i], region_d2_per_image[i])
  
  d <- Image(d)
  
  d <- colormap(d/nClusters, cbPalette[1:(nClusters+1)])
  
  filename <- paste0(clustering_folder, "/", image_name, "_region_cluster_map.png") 
  writeImage(d, filename)
  rm(d)
  
  # new start pos for next image
  st = st + regions_per_image[i]
  
}  
close(pb)

rm(data, data2, data3, data3m)

# Save everything so far
save.image(file = global_data_filename)

