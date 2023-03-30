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