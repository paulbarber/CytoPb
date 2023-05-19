# CytoPb helper function
# P R Barber, Mar 2023

# Identify negative and positive channel values
# image is a Cytoimagelist
# channel is a name of a channel
estimateNegPosValue <- function(image, channel, sigma = 10){
  
  # Get EBImage
  img <- image@listData[[1]][,,channel]
  img_blur <- gblur(img, sigma = sigma, boundary = 0) # to get the final intensity from
  
  img2 <- medianFilter(normalize(img), 2)
  
  #  img3 <- img2 > otsu(img2, range = c(0, max(img2)))
  img3 <- normalize(img2) > 0.5
  
  #  img4 <- opening(img3, kern = makeBrush(7, shape='disc'))  # to determine the foreground
  #  img4 <- opening(img3, kern = makeBrush(3, shape='disc'))  # to determine the foreground     **********************
  
  img4 <- bwlabel(img3)   # label 4-connected particles
  t <- table(img4)    # sizes of particles
  o <- names(t)[t<10]   # get list of small objects
  img4 <- rmObjects(img4, o)    # remove small objects, fg is left
  rm(t, o)
  
  image_name <- names(image)[1]
  filename <- paste0(channel_png_folder, image_name, "_", channel, "_fgmask.png")
  writeImage(img4, filename)
  
  
  # calculate negative and positive image values in the original image units
  # suppressWarnings on max in case there is no foreground
  #posv <- suppressWarnings((mean(img_blur[img4==1])))                              *************************
  posv <- suppressWarnings((quantile(img_blur[img4>0], 0.05)))    #  percentile, err on the low side
  negv <- suppressWarnings((mean(img_blur[img4==0])))
  
  # If posv comes out 0 (max returns -Inf), there was no foreground = NA
  # If posv comes out NaN (mean returns NaN), there was no foreground = NA
  if(is.nan(posv)) posv = NA
  else if(is.na(posv)) posv = NA
  else if(posv <= 0) posv = NA
  if(is.nan(negv)) negv = NA
  else if(is.na(negv)) negv = NA
  else if(negv <= 0) negv = NA
  
  # return
  c(negv, posv)
}
