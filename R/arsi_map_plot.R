#################################################
# Script to map sampling points in Arsi Negele  #
# Stephen Wood                                  #
# Last updated: May 21, 2018                    #
#################################################


# LOAD PACKAGES
library(sf)
library(tidyverse)
library(rgdal)
library(maptools)
library(ggplot2)

# READ IN DATA
## Farm points
farms <- readxl::read_excel("data/Geodata/GPS coordinate sampling points.xlsx")
names(farms)[1] <- 'ID'

## Convert points to sf object
farms_sf <- st_as_sf(farms, coords = c("Long", "Lat"), 
                     crs = 4326)

## Zones data
zone_1 <- st_read("data/Geodata/Zone1.shp")
zone_2 <- st_read("data/Geodata/Zone2.shp")
zone_3 <- st_read("data/Geodata/Zone3.shp")

# ## reformat problematic Munesa data (approach from Jeff Evans)
# m <- readOGR(getwd(), "Munesa_repair")
# 
# ## function to remove holes in sp polygon feature class
# remove.holes <- function(x) {
#   xp <- slot(x, "polygons")
#   holes <- lapply(xp, function(x) sapply(slot(x, "Polygons"), slot, "hole"))
#   res <- lapply(1:length(xp), function(i) slot(xp[[i]], "Polygons")[!holes[[i]]])
#   IDs <- row.names(x)
#   x.fill <- SpatialPolygons(lapply(1:length(res), function(i)
#     Polygons(res[[i]], ID=IDs[i])),
#     proj4string=CRS(proj4string(x)))
#   slot(x.fill, "polygons") <- lapply(slot(x.fill, "polygons"),
#                                      maptools::checkPolygonsHoles)  
#   slot(x.fill, "polygons") <- lapply(slot(x.fill, "polygons"),
#                                      "comment<-", NULL)
#   return( x.fill )     
# }
# 
# ## remove holes and create SpatialPolygonsDataFrame
# m <- remove.holes(m)
# m <- SpatialPolygonsDataFrame(m, data.frame(row.names="0", ID="1"))
# 
# ## write shapefile
# writeOGR(m, getwd(), "Munesa_bdy", driver="ESRI Shapefile",
#          check_exists=TRUE, overwrite_layer=TRUE)
# 
## read working munesa shapefile
munesa <- st_read("data/Geodata/Munesa_bdy.shp")
#new_munesa <- st_read("data/Geodata/New_Munesa.shp")
munesa <- st_transform(munesa,st_crs(zone_1))

## merge zone polygons together into one
zone_polygons <- rbind(zone_1,zone_2) %>%
                  rbind(zone_3)
zone_polygons <- zone_polygons[,!c(zone_polygons$Name,zone_polygons$FolderPath)]

## define custom bounding box
box <- st_polygon(list(rbind(c(38.71,7.22),c(38.87,7.22),c(38.87,7.34),c(38.71,7.34),
                            c(38.71,7.22))))
box = st_sfc(box) %>% st_set_crs(4326)

## subset munesa data set
new_munesa <- st_intersection(munesa,box)

# MAKE PLOTS
## Base plotting
plot(farms_sf["Type"],
     xlim=c(38.72,38.81),
     ylim=c(7.28,7.35),
     axes=TRUE,
     cex.axis=.6,
     key.pos=3,
     pch=20,
     main=NULL,
     reset=FALSE
)
plot(st_geometry(new_munesa),
     col="grey90",
     border="grey90",
     xlim=c(38.72,38.81),
     ylim=c(7.28,7.35),
     add=T
)
plot(st_geometry(zone_polygons),
     xlim=c(38.72,38.81),
     ylim=c(7.28,7.35),
     border='grey40',
     add=T
)
plot(farms_sf["Type"],
     xlim=c(38.72,38.81),
     ylim=c(7.28,7.35),
     pch=20,
     axes=TRUE,
     cex.axis=.6,
     main=NULL,
     add=T
)

## ggplot
ggplot() + 
  geom_sf(data=new_munesa) +
  geom_sf(data=zone_polygons) +
  geom_sf(aes(color=Type),data=farms_sf)
