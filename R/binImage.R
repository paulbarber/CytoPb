# Binning of an EBImage

# Reduce both linear dimensions of an image by a factor of 2
# by combining groups of 4 pixels.
# Combination methods can be any function like sum, mean, max etc.
# Sum would be used for traditional binning of image data.
# Will work with arrays but will turn them into images on return.

library(EBImage)

#x <- array(1:25, dim = c(5,5))
#x <- array(1:512, dim = c(16,12,2))

# Calculate the size of the new image if we bin by a certain amount
newBinDim <- function(x, bin_size = 2){
  
  d1 <- dim(x)[1]
  d2 <- dim(x)[2]
  
  #c(d1 - d1%%bin_size, d2 - d2%%bin_size) # the image size we can process
  c(floor(d1/bin_size), floor(d2/bin_size)) # the size of the result
}

binMatrix <- function(x, fun = sum, bin_size = 2){
  # x is a matrix (2D)
  
  # Make sure size is divisible by the bin size
  binned_dim <- newBinDim(x, bin_size = bin_size)
  new_dim <- binned_dim * bin_size
    
  # truncate
  x <- x[1:new_dim[1], 1:new_dim[2]]
  
  # reshape adding 2 extra dimensions of size bin_size
  # then apply fun along these dimensions
  # NB apply with simplify=TRUE will drop the summed dimension
  dim(x) <- c(bin_size, binned_dim[1], bin_size, binned_dim[2])
  apply(apply(x, c(2,3,4), fun), c(1,3), fun)
}

binImage <- function(x, fun = sum, bin_size = 2){
  # x is a matrix, or a set with many channels in a 3rd dimension.

  binned_dim <- newBinDim(x, bin_size = bin_size)
  
  if(length(dim(x)) == 2){
    y <- binMatrix(x, fun = fun, bin_size = bin_size)
  } else {
    y <- apply(x, 3, binMatrix, fun = fun, bin_size = bin_size)
    dim(y) <- c(binned_dim[1], binned_dim[2], dim(x)[3])
  }
  
  y <- Image(y)
  
  # maintain channel names
  dimnames(y)[3] <- dimnames(x)[3]
  
  y
}
