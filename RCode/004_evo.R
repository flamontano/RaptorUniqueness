
## EVOLUTIONARY DISTINCTIVENESS (ED) ## -----------
library(ape)
n_phy = readRDS("Data/raptors_trees.rds")

# Accipitriformes

species = raptor_traits[raptor_traits$Order == "Accipitriformes", ]$sp

phy = list()
for(i in 1:length(n_phy)){
  phy[[i]]<-drop.tip(n_phy[[i]],n_phy[[i]]$tip.label[-match(species, n_phy[[i]]$tip.label)])
}

#calculate evolutionary distinctiveness per species (ED) with fair.proportion option (Isaac et al. )
  pded_list = purrr::map(phy, .f = function(x){
  picante::evol.distinct(tree = x, type = "fair.proportion")
  }) # 
  names(pded_list) = paste0("tree_", 1:200)
  phy_ed = plyr::ldply(pded_list)
  phy_ed = group_by(phy_ed, Species) %>% 
    summarise_if(is.numeric, mean, na.rm = T)
  
  #unique especies 
  spp_unique<-phy_ed[phy_ed$w >= quantile(phy_ed$w, .9),]$Species#distinctive species
  
# calculating distinctiveness per community and number of distinctive species

  comm = comm_raptors[, colnames(comm_raptors) %in% phy_ed$Species]
  
  # convert to presence/absence data -- 
  comm[comm < 1] = 0
  comm[comm >= 1] = 1
  comm = comm[rowSums(comm) > 0, ] # remove empty cells
  comm= comm[, colSums(comm) > 0] # remove empty species
  comm = comm[rowSums(comm) > 0, ] # remove empty cells
  dim(comm)#250 spp
  
  # (0) overall species richness of lineage
  
  all_rich = data.frame(sites = rownames(comm), SR = rowSums(comm))
  all_rich = all_rich %>%
    tidyr::separate(sites, c("x", "y"), sep = "_")

## (1) distinct species richness 

  comm_spe = comm[, colnames(comm) %in% spp_unique]
  comm_spe = comm_spe[rowSums(comm_spe) > 0, ] # remove empty cells
  comm_spe = comm_spe[, colSums(comm_spe) > 0] # remove empty species

unique_rich = data.frame(sites = rownames(comm_spe), sr = rowSums(comm_spe))
unique_rich = unique_rich %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

ed_rich<-unique_rich %>% replace(is.na(.), 0)


# (2) distinctiveness per community

comm_ed<-comm[, colnames(comm) %in% phy_ed$Species]
setdiff(colnames(comm), phy_ed$Species)

for (i in 1:nrow(phy_ed)){
  comm_ed[,i]<-replace(comm_ed[,i], comm_ed[,i] == 1, phy_ed$w[i] )
}

## dataframe with all metrics
ed_df = data.frame(sites = rownames(comm_ed), ed_mean = rowMeans(comm_ed, na.rm=T))
ed_df = ed_df %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

ed_df<-dplyr::right_join(dplyr::select(coord_continents, x, y, country, continent) %>% 
                           mutate(x = as.character(x), y = as.character(y)),
                         ed_df,
                         by = c("x", "y"))

ed_df<-dplyr::left_join(ed_df, ed_rich,
                         by = c("x", "y"))

ed_df<-dplyr::left_join(ed_df, all_rich,
                         by = c("x", "y"))
ed_df<-ed_df %>% mutate(sr = ifelse(is.na(sr), 0, sr))

ed_df$residuals<-residuals(glm(ed_df$sr~ed_df$SR, na.action="na.exclude"))

## saving all for lineage
accip_evo<-list(phy_ed, spp_unique, ed_df)


# Falconiformes

species = raptor_traits[raptor_traits$Order == "Falconiformes", ]$sp

phy = list()
for(i in 1:length(n_phy)){
  phy[[i]]<-drop.tip(n_phy[[i]],n_phy[[i]]$tip.label[-match(species, n_phy[[i]]$tip.label)])
}

#calculate evolutionary distinctiveness per species (ED) with fair.proportion option (Isaac et al. )
pded_list = purrr::map(phy, .f = function(x){
  picante::evol.distinct(tree = x, type = "fair.proportion")
}) # 
names(pded_list) = paste0("tree_", 1:200)
phy_ed = plyr::ldply(pded_list)
phy_ed = group_by(phy_ed, Species) %>% 
  summarise_if(is.numeric, mean, na.rm = T)

#unique especies 
spp_unique<-phy_ed[phy_ed$w >= quantile(phy_ed$w, .9),]$Species#distinctive species

# calculating distinctiveness per community and number of distinctive species

comm = comm_raptors[, colnames(comm_raptors) %in% phy_ed$Species]

# convert to presence/absence data -- 
comm[comm < 1] = 0
comm[comm >= 1] = 1
comm = comm[rowSums(comm) > 0, ] # remove empty cells
comm= comm[, colSums(comm) > 0] # remove empty species
comm = comm[rowSums(comm) > 0, ] # remove empty cells
dim(comm)# 64 spp

# (0) overall species richness of lineage

all_rich = data.frame(sites = rownames(comm), SR = rowSums(comm))
all_rich = all_rich %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

## (1) distinct species richness 

comm_spe = comm[, colnames(comm) %in% spp_unique]
comm_spe = comm_spe[rowSums(comm_spe) > 0, ] # remove empty cells
comm_spe = comm_spe[, colSums(comm_spe) > 0] # remove empty species

unique_rich = data.frame(sites = rownames(comm_spe), sr = rowSums(comm_spe))
unique_rich = unique_rich %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

ed_rich<-unique_rich %>% replace(is.na(.), 0)


# (2) distinctiveness per community

comm_ed<-comm[, colnames(comm) %in% phy_ed$Species]
setdiff(colnames(comm), phy_ed$Species)

for (i in 1:nrow(phy_ed)){
  comm_ed[,i]<-replace(comm_ed[,i], comm_ed[,i] == 1, phy_ed$w[i] )
}

## dataframe with all metrics
ed_df = data.frame(sites = rownames(comm_ed), ed_mean = rowMeans(comm_ed, na.rm=T))
ed_df = ed_df %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

ed_df<-dplyr::right_join(dplyr::select(coord_continents, x, y, country, continent) %>% 
                           mutate(x = as.character(x), y = as.character(y)),
                         ed_df,
                         by = c("x", "y"))

ed_df<-dplyr::left_join(ed_df, ed_rich,
                        by = c("x", "y"))

ed_df<-dplyr::left_join(ed_df, all_rich,
                        by = c("x", "y"))
ed_df<-ed_df %>% mutate(sr = ifelse(is.na(sr), 0, sr))

ed_df$residuals<-residuals(glm(ed_df$sr~ed_df$SR, na.action="na.exclude"))

## saving all for lineage
falcon_evo<-list(phy_ed, spp_unique, ed_df)



# Strigiformes

species = raptor_traits[raptor_traits$Order == "Strigiformes", ]$sp
species<-species[!(spp_morpho$spp == "Glaucidium_mooreorum")& 
                           !(spp_morpho$spp =="Otus_feae")&
                           !(spp_morpho$spp =="Otus_moheliensis")& 
                           !(spp_morpho$spp =="Strix_hadorami")]#species with traits but no distribution data nor in phylogeny

phy = list()
for(i in 1:length(n_phy)){
  phy[[i]]<-drop.tip(n_phy[[i]],n_phy[[i]]$tip.label[-match(species, n_phy[[i]]$tip.label)])
}

#calculate evolutionary distinctiveness per species (ED) with fair.proportion option (Isaac et al. )
pded_list = purrr::map(phy, .f = function(x){
  picante::evol.distinct(tree = x, type = "fair.proportion")
}) # 
names(pded_list) = paste0("tree_", 1:200)
phy_ed = plyr::ldply(pded_list)
phy_ed = group_by(phy_ed, Species) %>% 
  summarise_if(is.numeric, mean, na.rm = T)

#unique especies 
spp_unique<-phy_ed[phy_ed$w >= quantile(phy_ed$w, .9),]$Species#distinctive species

# calculating distinctiveness per community and number of distinctive species

comm = comm_raptors[, colnames(comm_raptors) %in% phy_ed$Species]

# convert to presence/absence data -- 
comm[comm < 1] = 0
comm[comm >= 1] = 1
comm = comm[rowSums(comm) > 0, ] # remove empty cells
comm= comm[, colSums(comm) > 0] # remove empty species
comm = comm[rowSums(comm) > 0, ] # remove empty cells
dim(comm)# 64 spp

# (0) overall species richness of lineage

all_rich = data.frame(sites = rownames(comm), SR = rowSums(comm))
all_rich = all_rich %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

## (1) distinct species richness 

comm_spe = comm[, colnames(comm) %in% spp_unique]
comm_spe = comm_spe[rowSums(comm_spe) > 0, ] # remove empty cells
comm_spe = comm_spe[, colSums(comm_spe) > 0] # remove empty species

unique_rich = data.frame(sites = rownames(comm_spe), sr = rowSums(comm_spe))
unique_rich = unique_rich %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

ed_rich<-unique_rich %>% replace(is.na(.), 0)


# (2) distinctiveness per community

comm_ed<-comm[, colnames(comm) %in% phy_ed$Species]
setdiff(colnames(comm), phy_ed$Species)

for (i in 1:nrow(phy_ed)){
  comm_ed[,i]<-replace(comm_ed[,i], comm_ed[,i] == 1, phy_ed$w[i] )
}

## dataframe with all metrics
ed_df = data.frame(sites = rownames(comm_ed), ed_mean = rowMeans(comm_ed, na.rm=T))
ed_df = ed_df %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

ed_df<-dplyr::right_join(dplyr::select(coord_continents, x, y, country, continent) %>% 
                           mutate(x = as.character(x), y = as.character(y)),
                         ed_df,
                         by = c("x", "y"))

ed_df<-dplyr::left_join(ed_df, ed_rich,
                        by = c("x", "y"))

ed_df<-dplyr::left_join(ed_df, all_rich,
                        by = c("x", "y"))
ed_df<-ed_df %>% mutate(sr = ifelse(is.na(sr), 0, sr))

ed_df$residuals<-residuals(glm(ed_df$sr~ed_df$SR, na.action="na.exclude"))

## saving all for lineage
strig_evo<-list(phy_ed, spp_unique, ed_df)


## plotting Evolutionary Distinctiveness ##

library(ggplot2)
library("RColorBrewer")

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))


p.0 = plot_world_eqaul_area(color_polygon = "grey90") +
  viridis::scale_fill_viridis(direction = -1) +
  theme(legend.position =  c(0.55,0),
        legend.direction = "horizontal", 
        legend.key.height = unit(0.5, 'lines'),
        legend.key.width = unit(1.5, 'lines'),
        plot.margin = margin(-0.5, -0.1, -0.5, -0.1, "cm"))

p.1 = plot_world_eqaul_area(color_polygon = "grey90") +
  viridis::scale_fill_viridis() +
  theme(legend.position =  c(0.55,0),
        legend.direction = "horizontal", 
        legend.key.height = unit(0.5, 'lines'),
        legend.key.width = unit(1.5, 'lines'),
        plot.margin = margin(-0.5, -0.1, -0.5, -0.1, "cm"))

#Example for Accipitriformes

unique_plot<-accip_evo[[3]] #read same element for each lineage
unique_plot<-unique_plot[unique_plot$SR>0,] #delete all sites with no species
unique_plot$x<-as.numeric(unique_plot$x)
unique_plot$y<-as.numeric(unique_plot$y)

# unique richness
p_ed = p.0 + geom_tile(data = filter(unique_plot, !is.na(continent)),
                           aes(x = x, y = y, fill = sr), inherit.aes = T) +
  labs(fill = "SR")+ ggtitle("Evolutionary Distinct species")+theme(plot.title = element_text(hjust = 0.5, size =15), plot.margin = unit(c(-1, -1, -0.5, -1), "cm"))

# unique richness -  residuals
p_ed_res = p.0 + geom_tile(data = filter(unique_plot, !is.na(continent)),
                               aes(x = x, y = y, fill = residuals), inherit.aes = T) +
  labs(fill = expression(SR[Res]))+ ggtitle("Evolutionary Distinct species")+theme(plot.title = element_text(hjust = 0.5, size =15), plot.margin = unit(c(-1, -1, -0.5, -1), "cm"))

# community uniquenss (CED)

p_morphot_comm = p.1 + geom_tile(data = filter(unique_plot, !is.na(continent)),
                                 aes(x = x, y = y, fill = ed_mean), inherit.aes = F) +
  labs(fill = 'CED ')+ ggtitle("Evolutionary Distinctiveness") +theme(plot.title = element_text(hjust = 0.5, size =15), plot.margin = unit(c(-1, -1, -0.5, -1), "cm"))

# -----------------------------------------------------------


