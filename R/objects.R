# CytoPb map object functions
# P R Barber, Mar 2023

# Saving and loading the channel and cell type maps
# requires objects_folder and celltype_objects_folder to be defined 

channel_probability_maps_name <- "channel_probability_maps"
celltype_probability_maps_name <- "celltype_probability_maps"


saveChannelMapObject <- function(channel_probability_maps, image_name, channel_names){
  
  # set dim names for collections of marker maps
  dimnames(channel_probability_maps)[[3]] <- channel_names
  
  filename <- paste0(objects_folder, 
                   image_name, 
                   "_",
                   channel_probability_maps_name,
                   ".RData")
  
  save(channel_probability_maps, file = filename)
}

loadChannelMapObject <- function(image_name){
  
  # load the channel prob map we need
  load(paste0(objects_folder, 
              image_name, 
              "_",
              channel_probability_maps_name,
              ".RData"))
  
  # return it
  get(channel_probability_maps_name)
}

saveCellTypeMapObject <- function(celltype_probability_maps, image_name, celltype_names){
  
  # set dim names for collections of cell types
  dimnames(celltype_probability_maps)[[3]] <- celltype_names
  
  filename <- paste0(celltype_objects_folder, 
                     image_name, 
                     "_",
                     celltype_probability_maps_name,
                     ".RData")
  
  save(celltype_probability_maps, file = filename)
}

loadCellTypeMapObject <- function(image_name){
  
  # load the cell type prob map we need
  load(paste0(celltype_objects_folder, 
              image_name, 
              "_",
              celltype_probability_maps_name,
              ".RData"))
  
  # return it
  get(celltype_probability_maps_name)
}



