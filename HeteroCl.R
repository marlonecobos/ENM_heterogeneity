distance.weight <- function(x, xy, tau) {
  s0 <- ((xy[,1] -x[1] )^2 + (xy[,2] -x[2] )^2)
  d3 <- exp(-s0/ (2 * tau^2))
  return(d3)
}


covariance <- function(y, weights) {
  # y is an m by n matrix
  # weights is length m
  # Returns the weighted covariance matrix of y (by columns).
  if (missing(weights)) return (cov(y))
  w <- zapsmall(weights / sum(weights)) # Standardize the weights
  #y.bar <- apply(y * w, 2, sum)         # Compute column means
  y.bar <- colSums(y*w)
  z <- t(y) - y.bar                     # Remove the means
  z %*% (w * t(z))  
}
##
correlation <- function(y, weights) {
  z <- covariance(y, weights)
  sigma <- sqrt(diag(z))       # Standard deviations
  z / (sigma %o% sigma)
}



gw.pca <- function(x, xy, y, tau) {
  # x is a vector denoting a location
  # xy is a set of locations as row vectors
  # y is an array of attributes, also as rows
  # tau is a bandwidth
  # Returns a `princomp` object for the geographically weighted PCA
  # ..of y relative to the point x.
  w <- distance.weight(x = x, xy = xy, tau = tau)
  princomp(covmat=covariance(y, w), scale=TRUE)
}


#######################################################
####################
pres <- suppressWarnings(stack(list.files('/Users/p.joseratauchi/Documents/Project-R/phyto-two-approach/capas/clima/', full.names=TRUE, pattern='.asc')))

ssize <- 1000
Heter <- pres[[1]]

y <- rasterToPoints(pres)
pTs <- y[sample(1:nrow(y), ssize), ]
points <- y[,1:2]

#y <- y[,-c(1,2)]
y <- pTs[,3:dim(pTs)[2]]
pTs <- pTs[,1:2]



library("doFuture")
registerDoFuture()
plan(multiprocess(workers=3))
system.time(
z <- foreach(i=1:nrow(points)) %dopar% { 
  p <- points[i, ]
gw.pca(x=p, xy=pTs, y[, c(1:3, 10:15)],tau= 0.084)$sdev
})

z1 <- do.call(rbind, z)
meeanpc <- apply(z1, 1, mean)

hetero <- Heter; hetero[!is.na(raster::values(hetero))] <- meeanpc

writeRaster(hetero, "HeteroTEMP.tif", overwrite=T)



###
library(scales)
h10000 <- Het

h10000[!is.na(raster::values(h10000))] <- rescale(na.omit(raster::values(h10000)), c(0,1))

plot(h10000)




zv <- do.call(rbind, z[[]]$sdv)










