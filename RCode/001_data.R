
rm(list=ls())

### Loading trait data ###

raptor_traits<-read.csv("Data/raptor_traits.csv", row.names = 1)

### coordinates of cells #

coord_continents = readRDS("Data/coord_continent.rds") %>% as_tibble()

# read raptor community data #

comm_raptors = readRDS("Data/comm_raptors.rds")

## read raptor phylogeny

n_phy = readRDS("Data/raptors_trees.rds")#552 spp in Jetz phylogeny


