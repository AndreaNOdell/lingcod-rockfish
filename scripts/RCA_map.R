library("ggplot2")
library(tidyverse)
theme_set(theme_classic())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library(rgeos)
library("maps")
library(rgdal)
RCA_coords = read.csv("raw_data/RCA_coords.csv")
world <- ne_countries(scale = "medium", returnclass = "sf")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
states <- cbind(states, st_coordinates(st_centroid(states)))
shp <- readOGR(dsn = "raw_data/stanford-yr534kg9160-shapefile/", layer = "yr534kg9160", stringsAsFactors = F)

sites = data.frame(
  longitude = RCA_coords$longitude, 
  latitude = RCA_coords$latitude)

jpeg('plots/RCA_boundary.jpg')
ggplot(data = world) +
  geom_sf() +
  geom_sf(data = states, fill = NA) + 
  #geom_polygon(data = shp, aes(x = long, y = lat)) +
  geom_path(data = sites, aes(x = longitude, y = latitude), size = 0.5) +
  coord_sf(xlim = c(-116.75, -125.75), ylim = c(32, 42.25), expand = FALSE)
dev.off()
