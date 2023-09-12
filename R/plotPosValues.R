# CytoPb helper function
# P R Barber, Mar 2023

# SNR plot from positive values
plotPosValues <- function(pos_table, pos_value_plot_filename){
  d <- tidyr::gather(pos_table, Image, positive.value, 3:dim(pos_table)[2], factor_key=TRUE)
  
  d$Channel <- factor(d$Channel, levels=unique(d$Channel)) # keep channel order in plot
  d$Image <- as.character(d$Image)
  
  # reduce text size when number of images is large
  n_images <- dim(pos_table)[2]-3
  rel_size = 1
  if(n_images > 50){
    rel_size = 50/n_images
  }
  
  pdf(pos_value_plot_filename)
  print(ggplot(d, aes(Channel, Image, fill = log10(positive.value))) + 
          geom_tile() + 
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 5),
                legend.position = "none",
                axis.text.y = element_text(size = rel(rel_size))))
  dev.off()
}

plotRanges <- function(pos_table, neg_table, pos_value_plot_filename){
  
  p <- tidyr::gather(pos_table, Image, positive.value, 3:dim(pos_table)[2], factor_key=TRUE)
  n <- tidyr::gather(neg_table, Image, negative.value, 3:dim(neg_table)[2], factor_key=TRUE)
  
  p$Channel <- factor(p$Channel, levels=unique(p$Channel)) # keep channel order in plot
  p$Image <- as.character(p$Image)
  n$Channel <- factor(n$Channel, levels=unique(n$Channel)) # keep channel order in plot
  n$Image <- as.character(n$Image)
  
  d <- merge(p, n)
  d$diff <- d$positive.value - d$negative.value
  d$mean <- (d$positive.value + d$negative.value) / 2
  d$range <- d$diff / d$mean # range normalised by the mean
  
  # reduce text size when number of images is large
  n_images <- dim(pos_table)[2]-3
  rel_size = 1
  if(n_images > 50){
    rel_size = 50/n_images
  }
  
  mn <- -0.0001 # so that anything negative shows red
  mx <- 2       # a nice max value for normalised range
  zero <- -mn/(mx-mn)
  
  pdf(pos_value_plot_filename)
  print(ggplot(d, aes(Channel, Image, fill = range)) + 
          geom_tile()  + 
          scale_fill_gradientn(colours = c("red", "white", "green"), 
                               values = c(0,zero,1), limits = c(mn, mx), oob = scales::squish) + 
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 5),
                legend.position = "none",
                axis.text.y = element_text(size = rel(rel_size))))
  dev.off()
}
