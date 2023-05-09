# CytoPb: Script 4a Cell Type Colocalisations
# P R Barber, Jan 2023
# Cell type colocalisations

# Input is the cell type tif files in celltype_tif
# You define the pairs of cell type you want to colocalise
# in the first_in_pair and second_in_pair lists (vectors)
# It will calculate the scale space overlap signature
# between the two cell types, and make an image at the 
# signature peak scale.
# output: signature and scale images in colocalisation_output
# Make a table of the scale space colocalisation values
# (or characteristic distance)
# output: Colocalisation Table.csv
# Make a plot of the scale space characteristic distances
# output: Colocalisation Scale.pdf

if(!exists("working_folder")){
  working_folder <- choose.dir(caption = "Select data folder")
}

cat(paste("CytoPb 4a Working in: ", working_folder, "\n"))

global_data_filename <- paste0(working_folder, "/CytoPb.RData")

# Read previous session
load(global_data_filename)

# Define pairs of cell type to compare in the list externally
if(!exists("cell_type_pairs")){
  stop("Please define the List of cell type pairs to compare in CellTypePairs, eg:\n
        > cell_type_pairs <- list(c(\"CD8cell\", \"CD4cell\"), c(\"Epith\", \"Progenitor\"))\n")
}else{
  ct_pairs = cell_type_pairs
  rm(cell_type_pairs)   # this stops in going into global_data_filename and being overwritten if changed and rerun
}


# Where to look for all the cell type maps
in_folder <- celltype_objects_folder


out_folder <- paste0(results_folder, "/colocalisation_output/")
dir.create(out_folder, showWarnings = F)

# output files
Colocalisation_table_filename <- paste0(results_folder, "/Colocalisation Table.csv")
Colocalisation_scale_filename <- paste0(results_folder, "/Colocalisation Scale.pdf")


suppressMessages(library(EBImage))
suppressMessages(library(ggplot2))
suppressMessages(library(strex))

img_list <- img_names

# make a version of an image at a given scale
scaleSpace <- function (img, t){

  sigma = sqrt(t)
  
  # filter size cannot be bigger than smallest image dimension
  size = min(2 * ceiling(3 * sigma) + 1, min(dim(img)))
  if((size %% 2) == 0) size = size - 1   # must be odd
  
  w <- makeBrush(size = size, shape = 'gaussian', sigma = sigma)
  
  filter2(img, w, boundary = "replicate")

}

# make the scale space signature of a pair of greyscale images
scaleSpaceSignature <- function (img1, img2) {

  max_scale = 20
  
  scale_t <- vector()
  scale_sigma <- vector()
  overlap <- vector()
  
  for (i in 0:max_scale){
    t = 2^i
    sigma = sqrt(t);
    
    # Create scale space images
    ss1 <- scaleSpace(img1, t)
    ss2 <- scaleSpace(img2, t)
    mean_1 <- mean(ss1)
    mean_2 <- mean(ss2)

    # Create an image of the pixel-by-pixel minimum of scale space images, this is the 'overlap', get statistics about the image
    mean_min <- min(ss1, ss2)
    
    #take average of normalising mean overlap by the red and green means
    overlap <- c(overlap, ((mean_min/mean_1) + (mean_min/mean_2)) / 2)
    
    scale_t <- c(scale_t, t)
    scale_sigma <- c(scale_sigma, sigma)
  }

  data.frame(scale_t, scale_sigma, overlap)
}

Image <- vector()
Cell_Type_Pair <- vector()
Peak_Scale_t <- vector()
Peak_Scale_sigma <- vector()
COM_Scale_t <- vector()
COM_Scale_sigma <- vector()
Range_Error <- vector()

pb = txtProgressBar(min = 0, max = length(img_list)*length(ct_pairs), 
                    initial = 0)

for (i in 1: length(ct_pairs)){
  
  for (j in 1:length(img_list)){
    
    img <- img_list[j]
    setTxtProgressBar(pb, length(img_list)*(i-1) + j)
    
    ct1 <- ct_pairs[[i]][1] 
    ct2 <- ct_pairs[[i]][2]
    ct_names <- paste0(ct1, "_", ct2)
    
        # load the cell type prob map we need
    ct <- loadCellTypeMapObject(img)
    img1 <- ct[,,which(dim_names == ct1)]
    img2 <- ct[,,which(dim_names == ct2)]
    rm(ct)

    sig <- scaleSpaceSignature(img1, img2)
    
    sig$overlap_gradient <- c(0, sig$overlap[-1] - sig$overlap[1:(length(sig$overlap)-1)])
    com_n <- sum(sig$overlap_gradient * 1:length(sig$overlap_gradient))
    com_d <- sum(sig$overlap_gradient)
    
    # plot overlap versus scale
    filename <- paste0(out_folder, img, "_", ct_names, "_signature.png")
    png(filename)
    suppressMessages(
      print(ggplot(sig, aes(x = scale_sigma, y = overlap)) +
      scale_x_log10() +
      geom_line())
    )
    dev.off()
  
    # Get peak and COM of gradient
    ii = which.max(sig$overlap_gradient)
    t = sig$scale_t[ii];
    com_t = 2^(com_n / com_d)
    
    # try to make sure peak is well contained in the scale limits used
    critical_value = sig$overlap_gradient[ii] / 3;
    bgn <- sig$overlap_gradient[1]
    end <- sig$overlap_gradient[length(sig$overlap_gradient)]
    # make sure no NaN values reach the if statement, which will throw an error
    if(is.nan(critical_value)) critical_value = 0
    if(is.nan(bgn)) bgn = critical_value
    if(is.nan(end)) end = critical_value
    
    if (critical_value > bgn && critical_value > end) {
      range_error = 0
    } else {
      range_error = 1
    }  
    
    Image <- c(Image, img)
    Cell_Type_Pair <- c(Cell_Type_Pair, ct_names)
    Peak_Scale_t <- c(Peak_Scale_t, t)
    Peak_Scale_sigma <- c(Peak_Scale_sigma, sqrt(t))
    COM_Scale_t <- c(COM_Scale_t, com_t)
    COM_Scale_sigma <- c(COM_Scale_sigma, sqrt(com_t))
    Range_Error <- c(Range_Error, range_error)
    
    # make scale space image at peak overlap
    ss1 <- scaleSpace(img1, t)
    ss2 <- scaleSpace(img2, t)
  
    ss <- rgbImage(red = ss1/max(ss1), green = ss2/max(ss2))
    filename <- paste0(out_folder, img, "_", ct_names, ".png")
    writeImage(ss, filename)
  }
  
}
close(pb)

data <- data.frame(Image, Cell_Type_Pair,
                   Peak_Scale_t, Peak_Scale_sigma, 
                   COM_Scale_t, COM_Scale_sigma, Range_Error)
write.csv(data, file = Colocalisation_table_filename, row.names = F)

# From previous simulations COM_scale_sigma is linearly proportional to the characteristic separation distance.
# See: Redefining_Colocalisation_Jan2023
data <- read.csv(file = Colocalisation_table_filename)

pdf(Colocalisation_scale_filename)
d <- subset(data, Range_Error == 0)
print(ggplot(d, aes(x = Image, y = COM_Scale_sigma, fill = Cell_Type_Pair)) +
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip())
dev.off()


