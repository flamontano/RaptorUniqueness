rm(list=ls())

# required packages and functions #

library(dplyr)
library(ape)
# devtools::install_github("daijiang/dli55")
library(Matrix)
library(taxize)
# install.packages("wikitaxa")
library(wikitaxa)
library(ggplot2)

if(!require(viridis)) install.packages("viridis")

## function to plot maps

# Eckert IV  equal-area projection 
# https://upload.wikimedia.org/wikipedia/commons/c/c5/Ecker_IV_projection_SW.jpg
plot_world_eqaul_area = function(rdata = "Data/raw/world_equal_area_data.RData", print = FALSE,
                                 proj = "eck4",
                                 color_polygon = "gray70", fill_polygon = "gray90",
                                 color_lat_long = "gray75", color_label = "gray50",
                                 size_label = 2){
  # ~~~~~~~~~~~ Project from long-lat to Eckert IV projection ~~~~~~~~~~~ #
  # spTransform() is used for shapefiles and project() in the case of data frame
  # for more PROJ.4 strings check the followings
  #   http://proj4.org/projections/index.html
  #   https://epsg.io/
  load(rdata)
  # world = map_data("world")
  lakes = map_data("lakes")
  # __ give the PORJ.4 string for Eckert IV projection
  PROJ <- paste0("+proj=", proj, " +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") 
  # or use the short form "+proj=eck4"
  # lakes = map_data("lakes")
  lakes[, 1:2] = rgdal::project(cbind(lakes$long, lakes$lat), PROJ)
  # __ project the shapefiles
  NE_countries.prj  <- sp::spTransform(NE_countries, CRSobj = PROJ)
  NE_graticules.prj <- sp::spTransform(NE_graticules, CRSobj = PROJ)
  NE_box.prj        <- sp::spTransform(NE_box, CRSobj = PROJ)
  
  # __ project long-lat coordinates columns for data frames 
  # (two extra columns with projected XY are created)
  prj.coord <- rgdal::project(cbind(lbl.Y$lon, lbl.Y$lat), proj = PROJ)
  lbl.Y.prj <- cbind(prj.coord, lbl.Y)
  names(lbl.Y.prj)[1:2] <- c("X.prj","Y.prj")
  
  prj.coord <- rgdal::project(cbind(lbl.X$lon, lbl.X$lat), proj = PROJ)
  lbl.X.prj <- cbind(prj.coord, lbl.X)
  names(lbl.X.prj)[1:2] <- c("X.prj","Y.prj")
  
  
  # ~~~~~~~~~~~ Plot, edit layers and add legend ~~~~~~~~~~~ #
  plt = ggplot() +
    geom_polygon(data = NE_countries.prj, 
                 aes(long, lat, group = group), 
                 colour = color_polygon, fill = fill_polygon, size = .25) +
    geom_polygon(data = lakes, 
                 aes(long, lat, group = group), 
                 colour = color_polygon, fill = fill_polygon, size = .25) +
    # Note: "Regions defined for each Polygons" warning has to do with fortify transformation. 
    # fortify might get deprecated in future!
    # alternatively, use use map_data(NE_countries) to transform to data frame and then use project() to change to desired projection.
    # add projected bounding box
    # geom_polygon(data = NE_box.prj, 
    #              aes(x = long, y = lat), 
    #              colour = "black", fill = "transparent", size = .25) +
    labs(x = NULL, y = NULL) +
    # # add graticules
    # geom_path(data = NE_graticules.prj, 
    #           aes(long, lat, group = group), 
    #           linetype = "dotted", colour = color_lat_long, size = .25) +
    # # add graticule labels - latitude and longitude
    # geom_text(data = lbl.Y.prj, # latitude
    #           aes(x = X.prj, y = Y.prj, label = lbl), 
    #           colour = color_label, size = size_label) +
    # geom_text(data = lbl.X.prj, # longitude
    #           aes(x = X.prj, y = Y.prj, label = lbl), 
    #           colour = color_label, size = size_label) +
  # __ Set aspect ratio
  # the default, ratio = 1 in coord_fixed ensures that one unit on the x-axis is the same length as one unit on the y-axis
  coord_fixed(ratio = 1) +
    theme_void() 
  plt
}

p = plot_world_eqaul_area(color_polygon = "white")


