#################################################
# Script to map sampling points in Arsi Negele  #
# Stephen Wood                                  #
# Last updated: May 31, 2018                    #
#################################################


# LOAD PACKAGES
library(sf)
library(tidyverse)
library(rgdal)
library(maptools)
library(raster)     # For scalebar

# READ IN DATA
## Farm points
farms <- readxl::read_excel("data/Geodata/GPS coordinate sampling points.xlsx")
names(farms)[1] <- 'ID'
farms$ID <- as.factor(farms$ID)

## Convert points to sf object
farms_sp = st_as_sf(farms, coords = c("Long", "Lat"), crs = 4326)

## Zones
zones <- st_read("data/Geodata/Villages_merged.kml")

## Arsi Negele
an_town <- st_read("data/Geodata/arsi-town.kml")

## Roads
rd <- st_read("data/Geodata/road.kml")
hwy <- st_read("data/Geodata/natl-hwy.kml")

# ### merge roads into one
# roads <- rbind(rd,hwy) 

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
munesa <- st_transform(munesa,st_crs(zones))


# MAKE PLOTS
## Base plotting
plot(st_geometry(hwy),
     xlim=c(38.645,38.81),
     ylim=c(7.285,7.34),
     lwd=2
)
plot(st_geometry(munesa),
     col="grey90",
     border="grey90",
     add=T
)
plot(st_geometry(rd),
     add=T
)
plot(st_geometry(an_town),
     col="white",
     lwd=0.4,
     add=T
)
plot(st_geometry(zones),
     border='grey40',
     lwd=0.5,
     add=T
)
plot(farms_sp["Type"],
     pch=20,
     key.pos=3,
     axes=TRUE,
     cex.axis=.6,
     cex=0.75,
     pal=c("#66c2a5","#fc8d62","#8da0cb"),
     main=NULL,
     add=T
)
scalebar(5, xy=click(), type='bar', lonlat=T, below="km", divs=4, lwd=0.25)


# ## ggplot
# ggplot() + 
#   geom_sf(data=new_munesa) +
#   geom_sf(data=zone_polygons) +
#   geom_sf(data=an_town) +
#   geom_sf(aes(color=Type),size=0.5,data=farms_sf) +
#   theme(
#     legend.title = element_blank(),
#     legend.position = "bottom",
#     panel.grid = element_blank(),
#     panel.grid.major = element_line(color = 'white'),
#     panel.background = element_rect(fill="white"),
#     plot.background = element_rect(fill="white")
# )
