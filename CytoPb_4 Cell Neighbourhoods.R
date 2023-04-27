# CytoPb: Script 4 Neighbourhood UMAP - EXPERIMENTAL/UNDER DEVELOPMENT
# P R Barber, Apr 2023

# Following from script 3 (using same working directory).
# Grid the images up into regions
# Calculate a UMAP plot of regions based on the marker probability maps
# NB umap optimised for clustering
# Cluster the umap
# Plot heatmap of marker map values for the determined clusters
# The aim here is to identify neighbourhood types, the 
# region_grid_size_um can be large.

suppressMessages(library(umap))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))

source("R/binImage.R")

# OPTIONS
if(!exists("neighbourhood_grid_size_um")){   # c("All", "Random", "High.Marker")
  grid_um = 50
}else{
  grid_um = neighbourhood_grid_size_um
  rm(neighbourhood_grid_size_um)   # this stops in going into global_data_filename and being overwritten if changed and rerun
}



if(!exists("working_folder")){
  working_folder <- choose.dir(caption = "Select data folder")
}

cat(paste("CytoPb 4 Working in:", working_folder, "\n"))

global_data_filename <- paste0(working_folder, "/CytoPb.RData")

# Read previous session
load(global_data_filename)



# Select the images to process, doing all together may be too much
# UMAP takes a while
# hdbscan clustering takes a while and needs lots of memory
images_to_process <- img_names   # Process all images
#images_to_process <- img_names[1]
#images_to_process <- img_names[35]


max_colours_for_legend = 50

grid_size = grid_um/image_scale_umperpixel

# folder to save to
clustering_folder <- paste0(working_folder, "neighbourhood_clustering/")
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
  
  # Get some cell type probability data
  filename <- paste0(objects_folder, 
                     image_name, 
                     "_celltype_probability_maps.RData")
  load(filename)
  image <- celltype_probability_maps
  rm(celltype_probability_maps)

  
  # Get rid of low probabilities?
  # <0.1
  
  # Different binning functions
  # We are binning the cell type probabilities
  # These are like area and seem to be the best
  fun <- sum   # Account for area of each cell type and strength
  #fun <- mean  # all regions are the same size, so should be the same as sum
  # These are like cell presence?
  #fun <- max   # The strongest probability in each area for each type
    
  d <- binImage(image, fun, bin_size = grid_size)
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


# Copy celltype names and make sure all the names are good for the plots below
colnames(data) <- make.names(names(ct_matrix))

# Set umap settings for pre-clustering, tight clusters
#custom.settings = umap.defaults
#custom.settings$n_neighbors = 30
#custom.settings$min_dist = 0.0001
#custom.settings$n_components=2
#custom.settings$random_state=42


cat(paste("Performing UMAP with", dim(data)[1], "regions...\n"))
#data.umap <- umap(data, config = custom.settings, fast_sgd=T)
data.umap <- umap(data)   # default settings


u <- as.data.frame(data.umap$layout)
names(u) <- c("UMAP1", "UMAP2")
data2 <- cbind(data, u)

cat("Saving UMAPs per cell type...\n")
pb = txtProgressBar(min = 0, max = length(colnames(data)), initial = 0)
pdf(paste0(working_folder, "/Neighbourhood Celltype umap Plots.pdf"))
i=1
for(ch in colnames(data)){
  setTxtProgressBar(pb, i)
  print(ggplot(data2, aes(x = UMAP1, y = UMAP2)) +
          geom_point(aes_string(colour = ch), size = 1, alpha = 0.05))
  i = i + 1
}
dev.off()
close(pb)


#cat("hdbscan clustering...\n")
#library(dbscan)
## minPts is min cluster size BUT also a smoothing factor!
##minPts = ceiling(max(4, dim(data)[1]/500))   # some estimate
##minPts = 30
#minPts = 10
## Cluster the UMAP plot
#cl <- hdbscan(u, minPts = minPts, verbose = T)
## OR Cluster the data
##cl <- hdbscan(data, minPts = minPts, verbose = T)
## NB Cluster labels start at ZERO, but 0 are "noise points"!!!!!
#cat(paste0("Number of noise pixels: ", sum(cl$cluster==0), " (", 
#          round(100*sum(cl$cluster==0)/length(cl$cluster)), "%)\n"))
#pdf(paste0(clustering_folder, "neighbourhood_clustering.pdf"))
#plot(cl, show_flat = T)   # simpler plot showing most stable clusters
#dev.off()


cat("kmeans clustering...\n")
# Cluster the UMAP plot
#cl <- kmeans(u, 20)
# OR Cluster the data
cl <- kmeans(data, 50)
# kmeans may not converge with default algorithm
if (cl$ifault==4) { 
  cl = kmeans(data, cl$centers, algorithm="MacQueen", iter.max=1000)
}

nClusters = max(cl$cluster)
clusters <- cl$cluster

# reorder the clusters putting similar ones together
# Similarity is based on cell type content
# so collect and aggregate this data first
# remove cluster 0 if there is one
data3 <- cbind(data, clusters)
colnames(data3) <- c(colnames(data), "cluster")

data3l <- tidyr::pivot_longer(as.data.frame(data3), !cluster, names_to = "channel", values_to = "value")
data3la <- aggregate(value ~ cluster + channel, data = data3l, mean)
data3w <- tidyr::pivot_wider(data3la, names_from = channel, values_from = value)
data3ws <- as.data.frame(sapply(subset(data3w, cluster!=0), scale))  # remove cluster 0
data3ws$cluster <- 1:nClusters  # don't want cluster label scaled

cl.cl <- hclust(dist(as.matrix(data3ws[,2:dim(data3ws)[2]])))
cl_order <- cl.cl$order

# reorder for all regions
for(i in 1:nClusters){   # make sure to exclude cluster 0 here
  clusters[cl$cluster == cl_order[i]] <- i
}


cat("Saving plots...\n")

data3 <- cbind(data, clusters)
colnames(data3) <- c(colnames(data), "cluster")
data4 <- cbind(data2, clusters)
colnames(data4)[length(colnames(data4))] <- "cluster"
data4$cluster <- factor(data4$cluster, levels = 1:(nClusters))

# cluster colours
cbPalette <- c("#000000", rainbow(nClusters+2))


pdf(paste0(working_folder, "/Neighbourhood Cluster umap Plot.pdf"))
ggp1 <- ggplot(data4, aes(x = UMAP1, y = UMAP2)) +
        geom_point(aes(colour = cluster), size = 1, alpha = 0.05) +
        scale_colour_manual(values=cbPalette[-1]) +
        guides(colour=guide_legend(override.aes=list(alpha=1, size=3)))

if(nClusters > max_colours_for_legend){
  print(ggp1 + theme(legend.position = "none"))
} else {
  print(ggp1)
}

if(min(clusters)==0){
  ggp2 <- ggplot(subset(data4, !is.na(cluster)), aes(x = UMAP1, y = UMAP2)) +
          geom_point(aes(colour = cluster), size = 1, alpha = 0.05) +
          scale_colour_manual(values=cbPalette[-1]) +
          ggtitle("Noise pixels omitted") +
          guides(colour=guide_legend(override.aes=list(alpha=1, size=3)))
  
  if(nClusters > max_colours_for_legend){
    print(ggp2 + theme(legend.position = "none"))
  } else {
    print(ggp2)
  }
}

dev.off()


# cluster heatmap
# This code seems overly complex, but is there a better way to aggregate across
# all regions and then scale within each cell type?
pdf(paste0(working_folder, "/Neighbourhood Cluster Heatmap.pdf"))

data3l <- tidyr::pivot_longer(as.data.frame(data3), !cluster, names_to = "channel", values_to = "value")
data3la <- aggregate(value ~ cluster + channel, data = data3l, mean)
data3w <- tidyr::pivot_wider(data3la, names_from = channel, values_from = value)
data3ws <- as.data.frame(sapply(data3w, scale))
data3ws$cluster <- data3w$cluster  # don't want cluster label scaled
data3ls <- tidyr::pivot_longer(as.data.frame(data3ws), !cluster, names_to = "channel", values_to = "value")

print(ggplot(data3la, aes(y = cluster, x = channel, fill = value)) +
        geom_tile() + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 5),
              legend.position = "none")  +
        ggtitle("mean"))

print(ggplot(data3ls, aes(y = cluster, x = channel, fill = value)) +
        geom_tile() + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 5),
              legend.position = "none") +
        ggtitle("mean z scaled") )

# NB data is already scaled above, get better contrast in the heatmap
# if we use the default to scale by rows as well
heatmap(as.matrix(data3ws[,2:dim(data3ws)[2]]), cexCol = 0.5, 
        scale = "row", Colv = NA, Rowv = NA)
#legend(x="right", legend=c("max", "med", "min"), 
#       fill=c(hcl.colors(3, "YlOrRd", rev = FALSE)),
#       border = FALSE, bty = "n")

data3l <- tidyr::pivot_longer(as.data.frame(data3), !cluster, names_to = "channel", values_to = "value")
data3la <- aggregate(value ~ cluster + channel, data = data3l, max)
data3w <- tidyr::pivot_wider(data3la, names_from = channel, values_from = value)
data3ws <- as.data.frame(sapply(data3w, scale))
data3ws$cluster <- data3w$cluster  # don't want cluster label scaled
data3ls <- tidyr::pivot_longer(as.data.frame(data3ws), !cluster, names_to = "channel", values_to = "value")

print(ggplot(data3la, aes(y = cluster, x = channel, fill = value)) +
        geom_tile() + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 5),
              legend.position = "none")  +
        ggtitle("max"))

print(ggplot(data3ls, aes(y = cluster, x = channel, fill = value)) +
        geom_tile() + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 5),
              legend.position = "none") +
        ggtitle("max z scaled") )

# NB data is already scaled above, get better contrast in the heatmap
# if we use the default to scale by rows as well
heatmap(as.matrix(data3ws[,2:dim(data3ws)[2]]), cexCol = 0.5, 
        scale = "row", Colv = NA, Rowv = NA)
#legend(x="right", legend=c("max", "med", "min"), 
#       fill=c(hcl.colors(3, "YlOrRd", rev = FALSE)),
#       border = FALSE, bty = "n")

if(min(clusters)==0){
  
  # rescale without cluster zero
  data3ws <- as.data.frame(sapply(subset(data3w, cluster!=0), scale))
  data3ws$cluster <- 1:nClusters  # don't want cluster label scaled
  data3ls <- tidyr::pivot_longer(as.data.frame(data3ws), !cluster, names_to = "channel", values_to = "value")
  
  print(ggplot(data3ls, aes(y = cluster, x = channel, fill = value)) +
          geom_tile() + 
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 5),
                legend.position = "none") +
          ggtitle("z scaled (no cluster zero)") )
  
  
  # NB data is already scaled above, get better contrast in the heatmap
  # if we use the default to scale by rows
  heatmap(as.matrix(data3ws[,2:dim(data3ws)[2]]), cexCol = 0.5, 
          scale = "row", Colv = NA, Rowv = NA)
  #legend(x="right", legend=c("max", "med", "min"), 
  #       fill=c(hcl.colors(3, "YlOrRd", rev = FALSE)),
  #       border = FALSE, bty = "n")
}

dev.off()
rm(data3l, data3la, data3w, data3ws, data3ls)

# make a table of areas
neib_image <- vector()
neib_cluster <- vector()
neib_area <- vector()
neib_density <- vector()


# Region cluster maps
pb = txtProgressBar(min = 0, max = length(images_to_process), initial = 0)
st = 1
for(i in 1:length(images_to_process)){
  
  setTxtProgressBar(pb, i)
  image_name <- images_to_process[i]
  
  # Get the region cluster labels for this image
  d <- data4[st:(st+regions_per_image[i]-1), "cluster"]
  
  # store areas for later
  neib_image <- c(neib_image, rep(image_name, nClusters))
  neib_cluster <- c(neib_cluster, 1:nClusters)
  areas <- as.vector(table(d)) # factor makes sure all cell types are represented
  neib_area <- c(neib_area, areas)
  densities <- areas/length(d)
  neib_density <- c(neib_density, densities)
  
  # Make NA cluster zero again so it becomes black
  d <- as.numeric(d)
  d[is.na(d)] <- 0
  
  # reshape cluster labels like an image
  #  dim(d) <- c(region_d1_per_image[i], region_d2_per_image[i])
  
  # Resize to original size and reshape
  # This took a little fiddling, repeat each element, reshape into cols,
  # repeat each col, and transpose.
  d <- rep(d, each = grid_size)
  dim(d) <- c(region_d1_per_image[i]*grid_size, region_d2_per_image[i])
  d <- t(apply(d, 1, rep, each = grid_size))

  d <- Image(d)
  
  # Resize to original size
  #d <- resize(d, region_d1_per_image[i]*grid_size) # This is no good, it interpolates
  
  d <- colormap(d/nClusters, cbPalette[1:(nClusters+1)])
  
  filename <- paste0(clustering_folder, "/", image_name, "_neighbourhood_cluster_map.png") 
  writeImage(d, filename)
  rm(d)
  
  # new start pos for next image
  st = st + regions_per_image[i]
  
}  
close(pb)

data5 <- data.frame(neib_image, neib_cluster, 
                   neib_area, neib_density)

names(data5) <- c("Image", "Cluster", 
                 "Area", "Density")

write.csv(data5, paste0(working_folder, "/Neighbourhood Cluster Totals.csv"), row.names = F)

# lock in a cell type order
data5$Cluster <- factor(data5$Cluster, levels = 1:nClusters)


# Plots of image cluster content

d <- data5

# reduce text size when number of images is large
n_images <- length(img_names)
rel_size = 1
if(n_images > 50){
  rel_size = 50/n_images
}


pdf(paste0(working_folder, "/Neighbourhood Cluster Totals Plot.pdf"))

ggp1 <- ggplot(d, aes(x = Image, y = Area, fill = Cluster)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values=cbPalette[-1]) +
        coord_flip() +
        theme(axis.text.y = element_text(size = rel(rel_size)))

ggp2 <- ggplot(d, aes(x = Image, y = Density, fill = Cluster)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values=cbPalette[-1]) +
        coord_flip() +
        theme(axis.text.y = element_text(size = rel(rel_size)))

ggp3 <- ggplot(d, aes(x = Image, y = Area, fill = Cluster)) +
        geom_bar(position = "fill",stat = "identity") +
        scale_fill_manual(values=cbPalette[-1]) +
        coord_flip() +
        theme(axis.text.y = element_text(size = rel(rel_size)))

ggp4 <- ggplot(d, aes(x = Image, y = Density, fill = Cluster)) +
        geom_bar(position = "fill",stat = "identity") +
        scale_fill_manual(values=cbPalette[-1]) +
        coord_flip() +
        theme(axis.text.y = element_text(size = rel(rel_size)))

if(nClusters > max_colours_for_legend){
  print(ggp1 + theme(legend.position = "none"))
  print(ggp2 + theme(legend.position = "none"))
  print(ggp3 + theme(legend.position = "none"))
  print(ggp4 + theme(legend.position = "none"))
} else {
  print(ggp1)
  print(ggp2)
  print(ggp3)
  print(ggp4)
}

dev.off()


rm(data2, data3)

# Save everything so far
save.image(file = global_data_filename)

