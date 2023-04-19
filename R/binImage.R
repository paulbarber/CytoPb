# Binning of an EBImage

# Reduce both linear dimensions of an image by a factor of 2
# by combining groups of 4 pixels.
# Combination methods can be any function like sum, mean, max etc.
# Sum would be used for traditional binning of image data.
# Will work with arrays but will turn them into images on return.

x <- array(1:25, dim = c(5,5))
x <- array(1:512, dim = c(16,16,2))

binMatrixBy2 <- function(x, fun = sum){
  # x is a matrix (2D)
  
  # If odd dimensions drop first row or col
  if((dim(x)[1] %% 2) != 0) {  # odd
    x <- x[-dim(x)[1],]
  }
  
  if((dim(x)[2] %% 2) != 0) {  # odd
    x <- x[,-dim(x)[2]]
  }
  
  
  # reshape adding 2 extra dimensions of size 2
  # then apply fun along these dimensions
  # NB apply with simplify=TRUE will drop the summed dimension
  dim(x) <- c(2, dim(x)[1]/2, 2, dim(x)[2]/2)
  apply(apply(x, c(2,3,4), fun), c(1,3), fun)
}

binImageBy2 <- function(x, fun = sum){
  # x is a matrix, or a set with many channels in a 3rd dimension.
  
  if(length(dim(x)) == 2){
    y <- binMatrixBy2(x, fun = fun)
  } else {
    y <- apply(x, 3, binMatrixBy2, fun = fun)
    dim(y) <- c(floor(dim(x)[1]/2), floor(dim(x)[2]/2), dim(x)[3])
  }
  
  y <- Image(y)
  
  # maintain channel names
  dimnames(y)[3] <- dimnames(x)[3]
  
  y
}


binImage <- function(x, fun = sum, bin_size = 2){
  
  while (bin_size > 1) {
    x <- binImageBy2(x, fun = fun)
    bin_size = bin_size / 2
  }
  
  x
}

