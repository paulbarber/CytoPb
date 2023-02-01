# Overlap Signal Maps: Script 2
# P R Barber, Jan 2023

# Following from script 1 (using same working directory).
# Plot the image positive values as a SNR indicator
# output: Positive Value Plot.pdf
# Normalise each needed channel to neg and pos values.
# Calculate marker probability map.
# output: scaled and map images to channel_png
# Save the complete workspace for follow on scripts.


# Read previous session
#load("Overlap Signal Maps.RData")

# Read in pos and neg values
pos_table <- read.csv("pos_value_table.csv")
neg_table <- read.csv("neg_value_table.csv")

jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

# SNR plot from positive values
d <- tidyr::gather(pos_table, Image, positive.value, 4:dim(pos_table)[2], factor_key=TRUE)

d$Channel <- factor(d$Channel, levels=unique(d$Channel)) # keep channel order in plot

pdf("Positive Value Plot.pdf")
print(ggplot(d, aes(Channel, Image, fill = log10(positive.value))) + 
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 5),
        legend.position = "none"))
dev.off()

# Channel images and density plots

# Create some names in the environment for the EBImage stacks of probability maps
for(i in 1:length(images)){
  assign(names(images)[i], NULL)
}

pb = txtProgressBar(min = 0, max = length(channels_needed), initial = 0)
for(j in 1:length(channels_needed)){
  
  # Channel of interest
  channel <- channels_needed[j]
  setTxtProgressBar(pb,j)
  
  for(i in 1:length(images)){
    
    image_name <- names(images)[i]

    i_p1 <- images@listData[[i]][,,channel]

    # Get scaling parameters, ## log10 these to apply to log10 then gblured data
    table_name <- image_name
    table_name <- gsub(" ", ".", table_name)   # to match table col
    table_name <- gsub("-", ".", table_name)   # to match table col
    nv1 <- neg_table[which(neg_table$Channel == channel), table_name]
    pv1 <- pos_table[which(pos_table$Channel == channel), table_name]

    # Check these values, if nv is NA set to 0, if pv is NA set to large number (>65k for 16 bit IMC)
    if(is.na(nv1)) nv1 = 0
    if(is.na(pv1)) pv1 = 100000

    # scale 0-1, mainly for display
    i_p1a <- (i_p1 - nv1)/(pv1 - nv1)

    # Threshold/clamp to get an estimate of the probability of being positive.
    i_p1a[i_p1a<0] <- 0
    i_p1a[i_p1a>1] <- 1

    filename <- paste0(folder, "/", image_name, "_", channel, "_scaled.png")
    writeImage(i_p1a, filename)

    
    # blur to account for cell size/position uncertainties
    i_p1 <- gblur(i_p1, sigma = sigma)

    # re-scale 0-1
    i_p1 <- (i_p1 - nv1)/(pv1 - nv1)
    
    # Threshold/clamp to get an estimate of the probability of being positive.
    i_p1[i_p1<0] <- 0
    i_p1[i_p1>1] <- 1

    filename <- paste0(folder, "/", image_name, "_", channel, "_map.png")
    y = colormap(i_p1, jet.colors(256))
    writeImage(y, filename)
    
    # Store these probability maps for later use, in order as channels_needed
    # (would be good if we could set dimnames of the EBImage to the channel names.
    # dimnames only supports 3 names. See dimnames(l))
    l <- get(image_name)
    l <- combine(l, i_p1)
    assign(image_name, l)
    
  }
}
close(pb)

# Save everything so far
save.image(file = "Overlap Signal Maps.RData")



