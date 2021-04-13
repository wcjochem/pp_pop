# 
# Preparing the Oshikoto synthetic data
# - aggregate to "building-level" populations
# - generate separate samples for simulation study
#
# April 2021
# Chris Jochem (w.c.jochem@soton.ac.uk)
#

options("rgdal_show_exportToProj4_warnings"="none")
library(rgdal)
library(rgeos)

# path to raw data
# synthetic pop: https://doi.org/10.3390/data3030030
dir <- "/home/jochem/Dropbox/swap/proj/pp_pop"

# load population data
# attribute [b_ID] is the nearest building feature id number
hh <- readOGR(paste0(dir, "/dat"), layer = "osh_hh_pts")
hh <- rebuild_CRS(hh)
  dim(hh) # n = 37928  27
  head(hh) # household-level info
  
# load boundary map
# source: GADM
bnd <- readOGR(dsn=paste0(dir, "/dat"),
               layer="osh_bound_POP") # layer="osh_bound_CLIP")
bnd <- rebuild_CRS(bnd)
# re-define non-standard WGS84 code
proj4string(bnd) <- CRS("+init=epsg:4326") # warning
  summary(bnd)
  
# clip HH to the GADM boundary
hh <- hh[bnd,] # drop 342 hh from outside
  dim(hh) # n = 36956
  
# get OSH buildings to update household x,y positions
bldg <- readOGR(paste0(dir, "/dat"), layer = "osh_buildings")
bldg_pts <- gCentroid(bldg, byid = TRUE)
coords <- coordinates(bldg_pts)
bldg <- cbind(bldg@data, coords)
  
# update HH data to create building-level
hh_bldg <- aggregate(list("pop"=hh$hhsize, "hh"=hh$n), by=list("b_ID"=hh$b_ID), sum)
# matching HH to nearest spatial building
hh <- merge(hh_bldg, bldg[,c("uid", "x", "y")], 
            by.x = "b_ID", by.y = "uid",  all.x = TRUE)
  
coordinates(hh) <- c('x','y')
proj4string(hh) <- CRS("+init=epsg:4326")
  
# project to UTM 34S with KM units
# utm_km <- CRS("+proj=utm +south +zone=34 ellps=WGS84 +units=km")
# see: https://spatialreference.org/ref/epsg/32734/prettywkt/
utm_km <- CRS("+proj=tmerc +ellps=WGS84 +lat_0=0 +lon_0=21 +k=0.9996 +x_0=800000 +y_0=2450000 +units=km")
hh <- spTransform(hh, utm_km)
bnd <- spTransform(bnd, utm_km)

  head(hh)

# write out population data
saveRDS(hh, file.path("/home/jochem/Documents/GitHub/pp_pop", "data", "pop.rds"))
# write out boundary file
saveRDS(bnd, file.path("/home/jochem/Documents/GitHub/pp_pop", "data", "bnd.rds"))
