setwd("C:/Users/Marlon/Documents/R/Heterogeneity/")
peru <- raster::stack(list.files("Peru/clima/", pattern = ".asc$", full.names = T)) 
raster::plot(peru[[1]])
ext <- raster::extent(-81, -78, -7, 0)

#cuba <- raster::stack(list.files("Bios/", pattern = ".asc$", full.names = T))  
#var.stack <- cuba

var.stack <- raster::crop(peru, ext)


dist.window <- raster::xres(var.stack) * 2
#ncell.window <- 2 
nlayers.mean <- dim(var.stack)[3] - 1
var.normalization <- TRUE
rescale.result <- FALSE
multi.parallel <- TRUE
comp.each <- 5000
n.cores <- 6

# simple parallel
system.time({hetero <- kuenm_heterogen(var.stack = var.stack, dist.window = dist.window)})

# multi parallel
multi.parallel <- TRUE
system.time({hetero <- kuenm_heterogen(var.stack = var.stack, dist.window = dist.window, nlayers.mean = nlayers.mean, 
                                       var.normalization = var.normalization, multi.parallel = multi.parallel, 
                                       n.cores = n.cores, comp.each = comp.each)})


x11()
raster::plot(hetero)

x11()
raster::plot(raster::raster("cuba_heterogeneity.tif"))

raster::writeRaster(hetero, filename = "cuba_heterogeneity.tif", format = "GTiff")
