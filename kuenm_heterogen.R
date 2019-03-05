#' Environmental heterogeneity RasterLayer calculation
#' 
#' @description kuenm_heterogen calculates a heterogeneity layer based on values of 
#' multiple environmental variables. Values of heterogeneity result from Geographically 
#' Weighted Principal Component Analyses. 
#' 
#' @param var.stack a RasterStack of environmental variables of region of interest.
#' @param dist.window (numeric) distance of to be weighted and used for the
#' moving window. If not defined, \code{ncell.window} must be set.
#' @param ncell.window (numeric) number of cells to be considered as distance 
#' and weighted to then be used for the moving window. Ignored if \code{dist.window}
#' is defined.
#' @param nlayers.mean (numeric) number of weighted principal components to be 
#' considered when calculating the mean heterogeneity value per pixel. If not defined,
#' all principal components will be used.
#' @param var.normalization (logical) whether or not to apply a noramlization procedure 
#' to variables when performing the Principal Component Analysis. Default = TRUE.
#' @param rescale.result (logical) whether or not values representing environmental
#' heterogeneity should be rescaled to values from 0 to 1. Default = TRUE.
#' @param multi.parallel (logical) whether or not increase the amount of data analyzed 
#' during parallel processes. If TRUE, the number of processes performed per processor 
#' will be defined by the \code{comp.each} argument. Recommended when the number of
#' cells of layers in \code{var.stack} with non-NA values is bigger than 50000. 
#' @param n.cores (numeric) number of processors to be used if multi.parallel = TRUE.
#' If not defined all processors will be used.
#' @param comp.each (numeric) number of processes performed in each processor. Valid 
#' only if multi.parallel = TRUE. Default = 5000.
#' 
#' @return 
#' A RasterLayer representing levels of evironmental heterogeneity in each pixel.
#' 
#' @details 
#' By default calculations are parallelized using the doFuture package. However, defining 
#' the \code{multi.parallel} argument will spead up analyses as more processes are 
#' performed per processor. Using this option represent higher demands of RAM, especially
#' if number in \code{comp.each} is bigger than 5000.  
#' 
#' @examples 
#' #

kuenm_heterogen <- function(var.stack, dist.window, ncell.window, nlayers.mean, 
                            var.normalization = TRUE, rescale.result = TRUE, multi.parallel = FALSE, 
                            n.cores, comp.each = 5000) {
  
  # test for potential problems
  if (missing(var.stack)) {
    stop("var.stack is needed. Check function's help for details.")
  }
  if (missing(dist.window) & missing(ncell.window)) {
    stop("Either dist.window or ncell.window are needed. Check function's help for details.")
  }
  if (nlayers.mean > dim(var.stack)[3]) {
    warning("nlayers.mean is higher than the number of layers in var.stack.\nAll values will be considered for mean calculations.")
    nlayers.mean <- dim(var.stack)[3]
  }
  
  # preparing data
  cat("\nPreparing data for analyses...\n")
  
  heter <- var.stack[[1]] # layer for refilling values
  
  xy_values <- raster::rasterToPoints(var.stack) # raster to coordinates and values
  
  xy_coordinates <- xy_values[, 1:2] # only coordinates
  xy_values <- xy_values[, 3:dim(xy_values)[2]] # only values
  
  # preparing other arguments if needed
  if (missing(dist.window)) { # defining distance if ncel.window was used instead
    dist.window <- ncell.window * raster::xres(heter)
  }
  
  if (missing(nlayers.mean)) {
    nlayers.mean <- dim(var.stack)[3]
  }
  
  if (var.normalization == TRUE){
    #xy_values <- scale(xy_values) 
    xy_values <- sapply(1:dim(xy_values)[2],function(x) (xy_values[, x] - min(xy_values[, x])) / 
                          (max(xy_values[, x]) - min(xy_values[, x])))
  }
  
  if (missing(n.cores)) {
    n.cores <- parallel::detectCores()
  }
  
  # performing claculations
  if (multi.parallel == FALSE) {
    suppressPackageStartupMessages(library(doFuture))
    
    cat("\nRunning calculations, please wait...\n")
    
    doFuture::registerDoFuture()
    future::plan(future::multiprocess(workers = n.cores))
    
    gw_pcas <- foreach(i = 1:nrow(xy_coordinates)) %dopar% { 
      kuenm_gwpca(xy.point = xy_coordinates[i, ], xy.coordinates = xy_coordinates, 
                  xy.values = xy_values, dist.window = dist.window, var.normalization = FALSE)$sdev
    }
    
  } else {
    suppressPackageStartupMessages(library(future))
    
    cat("\nRunning calculations, please wait...\n")
    
    future::plan(future::multiprocess(workers = n.cores))
    gw_pcas <- new.env()
    
    steps <- seq(1, dim(xy_coordinates)[1], comp.each)
    kkk <- c(steps,  dim(xy_coordinates)[1] + 1)
    long_k <- length(kkk)
    
    pasos <- 1:(length(kkk) - 1)
    pasosChar <- paste0(pasos)
    
    # need to improve the way we store data, to avoid wasting time
    for (paso in pasosChar) {
      x <- as.numeric(paso)
      gw_pcas[[paso]] %<-% { # try with lists?
        seq_gwpca <- kkk[x]:(kkk[x + 1] - 1)
        
        gwpca <- lapply(seq_gwpca, function(y){
          gwp <- kuenm_gwpca(xy.point = xy_coordinates[y, ], xy.coordinates = xy_coordinates, 
                             xy.values = xy_values, dist.window = dist.window, var.normalization = FALSE)$sdev
          return(gwp)
        })
        
        return(gwpca)
      }
      avance <- (x / long_k) * 100
      cat("Computation progress: ", avance,"%" ,"\n")
    }
    
    gw_pcas <- as.list(gw_pcas)
    gw_names <- as.character(sort(as.numeric(names(gw_pcas))))
    gw_pcas <- gw_pcas[gw_names]
    gw_pcas <- do.call(c, gw_pcas)
  }
  
  future::plan(future::sequential)
  
  gw_pcas <- do.call(rbind, gw_pcas)
  gw_pcas <- apply(gw_pcas[, 1:nlayers.mean], 1, mean)
  
  if (rescale.result == TRUE){
    heter[!is.na(raster::values(heter))] <- scales::rescale(gw_pcas, to = c(0, 1)) 
  } else {
    heter[!is.na(raster::values(heter))] <- gw_pcas
  }
  
  return(heter)
}


#' Geographically weighted Principal Component Analyses
#' 
#' @description kuenm_gwpca performs a Geographically Weighted Principal 
#' Component Analysis for a given point
#' 
#' @param xy.point (numeric) a vector of a coordinate (x, y).
#' @param xy.coordinates (matrix) is a set of locations as row vectors
#' @param xy.values (matrix) an array of attributes, also as rows
#' @param dist.window (numeric) distance of to be weighted to be used for the
#' moving window.
#' @param var.normalization (logical) whether or not to apply a noramlization procedure 
#' to variables when performing the Principal Component Analysis. Default = TRUE.
#' 
#' @return 
#' A vector containing standar deviations for all principal components.
#' 
#' @details 
#' Geographically Weighted Principal Component Analyses is performed following
#' Author et al. (2010). 
#' 
#' @examples 
#' 

kuenm_gwpca <- function(xy.point, xy.coordinates, xy.values, dist.window, var.normalization = TRUE) {
  
  # testing for potential problems
  if (missing(xy.point)) {
    stop("xy.point is needed. Check function's help for details.")
  }
  if (missing(xy.coordinates)) {
    stop("xy.coordinates is needed. Check function's help for details.")
  }
  if (missing(xy.values)) {
    stop("xy.values is needed. Check function's help for details.")
  }
  if (missing(dist.window)) {
    stop("dist.window is needed. Check function's help for details.")
  }
  
  # weighted distance calculations
  g_dist <- ((xy.coordinates[, 1] - xy.point[1])^2 + (xy.coordinates[, 2] - xy.point[2] )^2)
  w_distance <- exp(-g_dist / (2 * dist.window^2))
  
  # weighted covariance matrix
  if (var.normalization == TRUE){
    #xy.values <- scale(xy.values) 
    xy.values <- sapply(1:dim(xy.values)[2],function(x) (xy.values[, x] - min(xy.values[, x])) / 
                          (max(xy.values[, x]) - min(xy.values[, x])))
  }
  
  st_weights <- zapsmall(w_distance / sum(w_distance)) # Standardize the weights
  y_bar <- colSums(xy.values * st_weights)
  xy_substracted <- t(xy.values) - y_bar                     # Substarct the sums
  cov_matrix <- xy_substracted %*% (st_weights * t(xy_substracted))
  
  # geographically weighted pca
  gwpca_result <- princomp(covmat = cov_matrix)
  
  return(gwpca_result)
}
