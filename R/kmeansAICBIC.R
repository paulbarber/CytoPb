# Calculate AIC and BIC for kmeans
# http://sherrytowers.com/2013/10/24/k-means-clustering/
# https://stackoverflow.com/questions/15839774/how-to-calculate-bic-for-k-means-clustering-in-r

kmeansAICBIC = function(fit){
  
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(data.frame(AIC = D + 2*m*k,
                  BIC = D + log(n)*m*k))
}
