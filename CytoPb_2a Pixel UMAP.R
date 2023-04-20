# CytoPb: Script 2a Pixel UMAP - EXPERIMENTAL/UNDER DEVELOPMENT
# P R Barber, Jan 2023

# Following from script 2 (using same working directory).
# Calculate a UMAP plot of pixels based on the marker probability maps
# NB umap optimised for clustering
# Cluster the umap
# Plot heatmap of marker map values for the determined clusters

library(umap)
library(ggplot2)
library(dbscan)
library(tidyr)

# OPTIONS
if(!exists("umap_pixel_selection")){   # c("All", "Random", "High.Marker")
  option = "High.Marker" 
}else{
  option = umap_pixel_selection
  rm(umap_pixel_selection)   # this stops in going into global_data_filename and being overwritten if changed and rerun
}

if(!exists("umap_max_n")){   # c("All", "Random", "High.Marker")
  max.n = 10000 
}else{
  max.n = umap_max_n
  rm(umap_max_n)   # this stops in going into global_data_filename and being overwritten if changed and rerun
}


if(!exists("working_folder")){
  working_folder <- choose.dir(caption = "Select data folder")
}

print("CytoPb 2a Working in:")
print(working_folder)

global_data_filename <- paste0(working_folder, "/CytoPb.RData")

# Read previous session
load(global_data_filename)


# folder to save to
clustering_folder <- paste0(working_folder, "pixel_clustering/")
dir.create(clustering_folder, showWarnings = F)


pb = txtProgressBar(min = 0, max = length(img_names), initial = 0)
for(i in 1:length(img_names)){
  
  setTxtProgressBar(pb, i)
  image_name <- img_names[i]
  
  # Get some channel probability data  # retrieve channel probability 
  filename <- paste0(objects_folder, 
                     image_name, 
                     "_channel_probability_maps.RData")
  load(filename)
  image <- channel_probability_maps
  rm(channel_probability_maps)
  
  if(option == "All" | option == "Random"){  # Form matrix of all data
    
    data <- matrix(nrow = length(image), ncol = length(panel_needed$name))
    
    for(i in 1:length(panel_needed$name)){
      img <- image@.Data[,,i]
      data[,i] <- as.vector(img)
    }
    
  }
  
  if(option == "Random"){
    # choose random pixels (rows) - otherwise UMAP take too long (10000 = 1 minute)
    data <- data[sample(nrow(data), max.n), ]
  }
  
  if(option == "High.Marker"){ # form a matrix of pixels high in at least one marker
    
    threshold = 0.9   # for pixel selection, one marker must be bigger than this
    high_marker = 0.1   # to suppress low marker values to zero
  
    # Get a pixel mask of where some marker is >threshold
    # apply max() to the 3rd dim of 3D channel data, preserve dims 1 and 2
    # if max > threshold? for each pixel
    pixel_mask <- apply(image@.Data, c(1,2), max) > threshold
    pixel_mask <- as.vector(pixel_mask)  # reshape to 1D vector of pixels
    
    d <- image@.Data     # get all the channel data [x,y,c]
    dim(d) <- c(dim(d)[1]*dim(d)[2], dim(d)[3])    # reshape to 2d [pixel,channel]
    
    # mask to get high pixels
    d <- d[pixel_mask,]
    
    # Now ignore all the low marker values
    d <- apply(d, c(1,2), function(x){ifelse(x>high_marker, x, 0)})
    
    data <- d
    rm(d)
    
  }
  
  if(dim(data)[1] > max.n){
    # choose random pixels (rows)
    data <- data[sample(nrow(data), max.n), ]
  }
  
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
  
  data.umap <- umap(data, config = custom.settings, fast_sgd=T)
  
  
  u <- as.data.frame(data.umap$layout)
  names(u) <- c("UMAP1", "UMAP2")
  data2 <- cbind(data, u)
  
  pdf(paste0(clustering_folder, image_name, "_umap_plots.pdf"))
  for(ch in colnames(data)){
    print(ggplot(data2, aes(x = UMAP1, y = UMAP2)) +
      geom_point(aes_string(colour = ch), size = 1))
  }
  dev.off()
  
  
  # cluster
  cl <- hdbscan(u, minPts = 100)
  
  pdf(paste0(clustering_folder, image_name, "_pixel_clustering.pdf"))
  plot(cl, show_flat = T)   # simpler plot showing most stable clusters
  dev.off()
  
  data3 <- cbind(data, cl$cluster)
  colnames(data3) <- c(colnames(data), "cluster")
  
  
  # cluster heatmap
  data3m <- tidyr::gather(as.data.frame(data3), key = "channel", value = "value", colnames(data), factor_key=TRUE)
  pdf(paste0(clustering_folder, image_name, "_pixel_cluster_heatmap.pdf"))
  print(ggplot(data3m, aes(y = cluster, x = channel, fill = value)) +
    geom_tile() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 5),
          legend.position = "none") )
  dev.off()


}
close(pb)




# Save everything so far
save.image(file = global_data_filename)

