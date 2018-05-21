# load pacakges
library(sf)
library(tidyverse)
library(rgdal)
library(maptools)
library(ggplot2)

# set directory
setwd("~/data/Geodata") 

# read in zones data
zone_1 <- st_read("Zone1.shp")
zone_2 <- st_read("Zone2.shp")
zone_3 <- st_read("Zone3.shp")

# reformat problematic Munesa data (approach from Jeff Evans)
m <- readOGR(getwd(), "Munesa_repair")

## function to remove holes in sp polygon feature class
remove.holes <- function(x) {
  xp <- slot(x, "polygons")
  holes <- lapply(xp, function(x) sapply(slot(x, "Polygons"), slot, "hole"))
  res <- lapply(1:length(xp), function(i) slot(xp[[i]], "Polygons")[!holes[[i]]])
  IDs <- row.names(x)
  x.fill <- SpatialPolygons(lapply(1:length(res), function(i)
    Polygons(res[[i]], ID=IDs[i])),
    proj4string=CRS(proj4string(x)))
  slot(x.fill, "polygons") <- lapply(slot(x.fill, "polygons"),
                                     maptools::checkPolygonsHoles)  
  slot(x.fill, "polygons") <- lapply(slot(x.fill, "polygons"),
                                     "comment<-", NULL)
  return( x.fill )     
}

## remove holes and create SpatialPolygonsDataFrame
m <- remove.holes(m)
m <- SpatialPolygonsDataFrame(m, data.frame(row.names="0", ID="1"))

## write shapefile
writeOGR(m, getwd(), "Munesa_bdy", driver="ESRI Shapefile",
         check_exists=TRUE, overwrite_layer=TRUE)

## read working munesa shapefile
munesa <- st_read("Munesa_bdy.shp")
new_munesa <- st_read("New_Munesa.shp")
munesa <- st_transform(munesa,st_crs(zone_1))

## merge zone polygons together into one
zone_polygons <- rbind(zone_1,zone_2) %>%
                  rbind(zone_3)

## define custom bounding box
box <- st_polygon(list(rbind(c(38.71,7.22),c(38.87,7.22),c(38.87,7.34),c(38.71,7.34),
                            c(38.71,7.22))))
box = st_sfc(box) %>% st_set_crs(4326)

## subset munesa data set
new_munesa <- st_intersection(munesa,box)

# make plots
plot(st_geometry(zone_polygons),axes=T)
plot(st_geometry(new_munesa),add=T,col="grey")

