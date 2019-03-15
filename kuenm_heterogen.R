#' Environmental heterogeneity RasterLayer calculation
#' 
#' @description kuenm_heterogen calculates a heterogeneity layer based on values of 
#' multiple environmental variables. Values of heterogeneity result from Geographically 
#' Weighted Principal Component Analyses. 
#' 
#' @param var.stack a RasterStack of environmental variables of region of interest.
#' @param dist.window (numeric) distance of to be weighted and used for the
#' moving window. If not defined, \code{ncell.window} must be set. The value must be greater
#' than 0.01...
#' @param ncell.window (numeric) number of cells to be considered as distance 
#' and weighted to then be used for the moving window. Ignored if \code{dist.window}
#' is defined.
#' @param nlayers.mean (numeric) number of weighted principal components to be 
#' considered when calculating the mean heterogeneity value per pixel. If not defined,
#' all principal components will be used.
#' @param var.scale (logical) whether or not to scale variables in \code{var.stack} 
#' when performing the Principal Component Analysis. Default = FALSE. See details.
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
#' The scaling procedure does not center the variables only scale them from zero to one. 
#' 
#' @examples 
#' 
#' 

kuenm_heterogen <- function(var.stack, dist.window, ncell.window, nlayers.mean, 
                            var.scale = TRUE, ex.fun = 2, rescale.result = TRUE, 
                            multi.parallel = FALSE, n.cores, comp.each = 10000) {
  
  # 1. Add method
  # Geographically weighted PCA
  #gwpca
  
  # Geographically weighted Regression
  #gwr
  
  # universal kriging
  #univ.kr
  
  # test for potential problems
  if (missingArg(var.stack)) {
    stop("var.stack is needed. Check function's help for details.")
  }
  if (missingArg(dist.window) & missingArg(ncell.window)) {
    stop("Either dist.window or ncell.window are needed. Check function's help for details.")
  }
  if (missingArg(nlayers.mean)) {
    nlayers.mean <- dim(var.stack)[3]
  } else {
    if (nlayers.mean > dim(var.stack)[3]) {
      warning("nlayers.mean is higher than the number of layers in var.stack.\nAll PCs will be considered for mean calculations.")
      nlayers.mean <- dim(var.stack)[3]
    }
  }
  
  
  # preparing data
  cat("\nPreparing data for analyses...\n")
  
  heter <- var.stack[[1]] # layer for refilling values
  
  xy_values <- raster::rasterToPoints(var.stack) # raster to coordinates and values
  
  xy_coordinates <- xy_values[, 1:2] # only coordinates
  xy_values <- xy_values[, 3:dim(xy_values)[2]] # only values
  
  # preparing other arguments if needed
  if(missingArg(dist.window))  { # defining distance if ncel.window was used instead
    dist.window <- ncell.window * raster::xres(heter)
  }
  
  if(missingArg(nlayers.mean)) {
    nlayers.mean <- dim(var.stack)[3]
  }
  
  if (var.scale == TRUE){
    xy_values <- sapply(1:dim(xy_values)[2], function(x) (xy_values[, x] - min(xy_values[, x])) / 
                          (max(xy_values[, x]) - min(xy_values[, x])))
  }
  
  if(missingArg(n.cores)) {
    n.cores <- future::availableCores()
  }
  
  # performing claculations
  if (multi.parallel == FALSE) {
    suppressPackageStartupMessages(library(future))
    suppressPackageStartupMessages(library(doFuture))
    
    cat("\nRunning calculations, please wait...\n")
    
    doFuture::registerDoFuture()
    
    if(.Platform$OS.type == "unix") {
      future::plan(future::multiprocess(workers = n.cores))
    } else {
      future::plan(future::multiprocess)
    }
    
    gw_pcas <- foreach(i = 1:nrow(xy_coordinates)) %dopar% { 
      kuenm_gwpca(xy.point = xy_coordinates[i, ], xy.coordinates = xy_coordinates, 
                  xy.values = xy_values, dist.window = dist.window, 
                  var.scale = FALSE, ex.fun = ex.fun)$sdev
    }
    
    gw_pcas <- do.call(rbind, gw_pcas)
    gw_pcas <- apply(gw_pcas[, 1:nlayers.mean], 1, mean)
    
    future::plan(future::sequential)
    
  } else {
    suppressPackageStartupMessages(library(future))
    suppressPackageStartupMessages(library(magrittr))
    
    cat("\nRunning calculations, please wait...\n")
    
    if(.Platform$OS.type == "unix") {
      future::plan(future::multiprocess(workers = n.cores))
    } else {
      future::plan(future::multiprocess)
    }
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
        
        gwpca <- seq_gwpca %>%
          purrr::map_df(~data.frame(as.list(kuenm_gwpca(xy.point = xy_coordinates[.x, ], 
                                                        xy.coordinates = xy_coordinates, 
                                                        xy.values = data.frame(xy_values), 
                                                        dist.window = dist.window, 
                                                        var.scale = FALSE,
                                                        ex.fun = ex.fun)$sdev[1:nlayers.mean]))) %>% rowMeans(.)
        
        cat(paste0(gwpca, "\n"), file = paste0("heterogeneity_temp_", x, ".csv"))
        
        return(NULL)
      }
      avance <- (x / long_k) * 100
      cat("Computation progress: ", avance,"%" ,"\n")
    }
    
    future::plan(future::sequential)
    
    nfiles <- length(list.files(pattern = "^heterogeneity_temp_"))
    
    while(nfiles < x){
      Sys.sleep(.1)
      nfiles <- length(list.files(pattern = "^heterogeneity_temp_"))
    }
    
    gw_pcas <- paste0("heterogeneity_temp_", pasos, ".csv") %>% 
      purrr::map_df(~read.csv(.x, header = F)) %>% .[,1]
    
    unlink(paste0("heterogeneity_temp_", pasos, ".csv"))
  }
  
  
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
#' @param var.scale (logical) whether or not to scale variables in \code{xy.values} 
#' when performing the Principal Component Analysis. Default = FALSE. See details.
#' @param xy.point (numeric) 
#' 
#' @return 
#' A vector containing standar deviations for all principal components.
#' 
#' @details 
#' Geographically Weighted Principal Component Analyses is performed following
#' Author et al. (2010). 
#' 
#' The scaling procedure does not center the variables only scale them from zero to one.
#' 
#' @examples 
#' 

kuenm_gwpca <- function(xy.point, xy.coordinates, xy.values, dist.window, ex.fun = 2, 
                        var.scale = FALSE) {
  
  # testing for potential problems
  if (missingArg(xy.point)) {
    stop("xy.point is needed. Check function's help for details.")
  }
  if (missingArg(xy.coordinates)) {
    stop("xy.coordinates is needed. Check function's help for details.")
  }
  if (missingArg(xy.values)) {
    stop("xy.values is needed. Check function's help for details.")
  }
  if (missingArg(dist.window)) {
    stop("dist.window is needed. Check function's help for details.")
  }
  
  # weighted distance calculations
  g_dist <- ((xy.coordinates[, 1] - xy.point[1])^2 + (xy.coordinates[, 2] - xy.point[2] )^2)
  w_distance <- exp(-g_dist / (2 * dist.window^ex.fun))
  
  # weighted covariance matrix
  if (var.scale == TRUE){
    xy.values <- sapply(1:dim(xy.values)[2], function(x) (xy.values[, x] - min(xy.values[, x])) / 
                          (max(xy.values[, x]) - min(xy.values[, x])))
  }
  
  st_weights <- zapsmall(w_distance / sum(w_distance)) # Standardize the weights
  y_bar <- colSums(xy.values * st_weights)
  xy_substracted <- t(xy.values) - y_bar               # Substarct the sums
  cov_matrix <- xy_substracted %*% (st_weights * t(xy_substracted))
  
  # geographically weighted pca
  gwpca_result <- princomp(covmat = cov_matrix)
  
  return(gwpca_result)
}
