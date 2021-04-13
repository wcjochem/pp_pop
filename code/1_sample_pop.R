# 
# Preparing the Oshikoto synthetic data
# - generate samples for simulation study
#
# April 2021
# Chris Jochem (w.c.jochem@soton.ac.uk)
#


options("rgdal_show_exportToProj4_warnings"="none")
library(rgdal)
library(rgeos)
library(ggplot2)
library(inlabru)

set.seed(1126)
utm_km <- CRS("+proj=tmerc +ellps=WGS84 +lat_0=0 +lon_0=21 +k=0.9996 +x_0=800000 +y_0=2450000 +units=km")

# path to prepared data
dir <- "/home/jochem/Documents/GitHub/pp_pop"

# load population (building points)
pop <- readRDS(file.path(dir, "data", "pop.rds"))
# load boundary
bnd <- readRDS(file.path(dir, "data", "bnd.rds"))

# generate random samples
reps <- 50
nsamp <- 15

# pre-allocate storage
designA <- vector("list", length = reps)
designB <- vector("list", length = reps)

for(i in 1:reps){
  print(i)
  
  pts <- spsample(bnd, n = nsamp, type = "random")
  
  # create sample sites (Design A)
  # buffer each site a random size between 10 and 12 km
  pts$buff <- sample(10:15, size = nsamp, replace = TRUE)
  enumbuff <- gBuffer(pts, byid = TRUE, width = pts$buff)
  # dissolve
  enum <- gUnaryUnion(enumbuff)
  # clip to boundary
  enum <- gIntersection(enum, bnd, byid = TRUE)
  # remove holes
  enum <- SpatialPolygons(lapply(1:length(enum), function(i){ 
    outerRings <- Filter(function(f){f@ringDir == 1}, enum@polygons[[i]]@Polygons)
    Polygons(outerRings, ID = i)
  }), proj4string = utm_km)
  
  # create opposite areas of fill (Design B)
  negenum <- gSymdifference(enum, bnd)
  
  # get samples
  pop.samp <- pop[enum,]
  pop.fill <- pop[negenum,]
  
  # store sample size info
  designA[[i]] <- c(sum(pop.samp$pop), sum(pop.samp$hh), nrow(pop.samp))
  designB[[i]] <- c(sum(pop.fill$pop), sum(pop.fill$hh), nrow(pop.fill))
  
  # create sample files
  saveRDS(pop.samp, file = file.path(dir, "data", "sim", paste0("A_", i, ".rds")))
  saveRDS(pop.fill, file = file.path(dir, "data", "sim", paste0("B_", i, ".rds")))
}

ggplot() + gg(bnd) + gg(enum) + coord_equal() + gg(pop.samp, alpha=0.2)
ggplot() + gg(bnd) + gg(enum) + coord_equal() + gg(pop.fill, alpha=0.2)

designA <- do.call(rbind.data.frame, designA)
names(designA) <- c("samp_pop", "samp_hh", "samp_buildings")
designB <- do.call(rbind.data.frame, designB)
names(designB) <- c("samp_pop", "samp_hh", "samp_buildings")

designA$samp_pop / sum(pop$pop)
designB$samp_pop / sum(pop$pop)

# output sample stats
saveRDS(designA, file = file.path(dir, "data", "designA.rds"))
saveRDS(designB, file = file.path(dir, "data", "designB.rds"))
