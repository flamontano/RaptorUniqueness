## DIET UNIQUENESS ##

library(ecodist)
library(vegan)
library(hypervolume)
library(compositions)
library(factoextra)
library(ggplot2)
library(plyr)
library(dplyr)

# Accipitriformes ## -----------------------------------------------------------------------------------

diet_traits<-raptor_traits[raptor_traits$Order == "Accipitriformes", c(13:20)]
diet_traits<-diet_traits[complete.cases(diet_traits), ]#238 spp
names(diet_traits)<-c("Invertebrates", "Vert Endotherms", "Vert Ectotherms", "Fish", "Vert unk", "Scavenger", "Fruit", "Plant")

#deleting species with 100% of diet unkwnon
diet_traits<-diet_traits[diet_traits$`Vert unk`!= 100,]# 237 spp 

# pca 
diet.pca <-princomp(clr(acomp(diet_traits)))
summary(diet.pca)

# visualizing dietary space
funspace_diet<-funspace::funspace(diet.pca, PCs = c(1, 2))
plot(funspace_diet, 
     arrows =T, arrows.length = 1.6, arrows.label.cex = 0.7, arrows.label.pos = 1.1, arrows.label.col = "black", arrows.head = 0.05,
     pnt=T, pnt.col = "black", pnt.cex=0.2)

abline(v = 0, lty = "dashed", col = "grey70")
abline(h = 0, lty = "dashed", col = "grey70")


## species probability density in dietary space (DU PrDn)

spp_scores<-as.data.frame(diet.pca$scores[,1:2])#only PCA1 and PCA2 
hv<-hypervolume_gaussian(spp_scores, samples.per.point = 2000, quantile.requested = 1)

spp_diet<-as.data.frame(hypervolume_estimate_probability(hv, spp_scores));names(spp_diet)<-"uniqueness"
spp_diet$spp<-rownames(diet.pca$scores)

#unique especies 
spp_unique<-spp_diet[spp_diet$uniqueness <= quantile(spp_diet$uniqueness, .1),]$spp

# uniqueness per community and number of unique species

comm = comm_raptors[, colnames(comm_raptors) %in% row.names(diet_traits)]

# convert to presence/absence data -- 
comm[comm < 1] = 0
comm[comm >= 1] = 1
comm = comm[rowSums(comm) > 0, ] # remove empty cells
comm= comm[, colSums(comm) > 0] # remove empty species
comm = comm[rowSums(comm) > 0, ] # remove empty cells
dim(comm)#237 spp


# (0) overall species richness of lineage

all_rich = data.frame(sites = rownames(comm), SR = rowSums(comm))
all_rich = all_rich %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

# (1) unique species richness

# -- PCA! and PCA2
comm_spe = comm[, colnames(comm) %in% spp_unique]
comm_spe = comm_spe[rowSums(comm_spe) > 0, ] # remove empty cells
comm_spe = comm_spe[, colSums(comm_spe) > 0] # remove empty species

diet_rich = data.frame(sites = rownames(comm_spe), sr = rowSums(comm_spe))
diet_rich = diet_rich %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

diet_rich<-diet_rich %>% replace(is.na(.), 0)

## (2) uniqueness per community

comm_uniq<-comm[, colnames(comm) %in% spp_diet$spp]#PCA1 and PCA2
setdiff(colnames(comm_uniq), spp_diet$spp)

#replacing cells with uniqueness scores
for (i in 1:nrow(spp_diet)){
  comm_uniq[,i]<-replace(comm_uniq[,i], comm_uniq[,i] == 1, spp_diet$uniqueness[i] )
}

## dataframe with all metrics
diet_df = data.frame(sites = rownames(comm_uniq), diet_mean = rowMeans(comm_uniq))
diet_df = diet_df %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

diet_df<-dplyr::right_join(dplyr::select(coord_continents, x, y, country, continent) %>% 
                              mutate(x = as.character(x), y = as.character(y)),
                            select(diet_df, x, y, diet_mean),
                            by = c("x", "y"))
diet_df<-dplyr::full_join(diet_df, diet_rich,
                          by = c("x", "y"))
diet_df<-dplyr::left_join(diet_df, all_rich,
                          by = c("x", "y"))

diet_df[,6:7]<-diet_df[,6:7] %>% replace(is.na(.), 0)#for NA in richness

diet_df$residuals<-residuals(glm(data = diet_df, sr~SR ))

#saving all results for clade 

accip_diet<-list(diet.pca, spp_diet, hv, spp_unique, diet_df)


## Falconiformes ## ----------------------------------------------------------------------------------------

diet_traits<-raptor_traits[raptor_traits$Order == "Falconiformes", c(13:20)]
diet_traits<-diet_traits[complete.cases(diet_traits), ]#63 spp
names(diet_traits)<-c("Invertebrates", "Vert Endotherms", "Vert Ectotherms", "Fish", "Vert unk", "Scavenger", "Fruit", "Plant")

#deleting species with 100% of diet unkwnon
diet_traits<-diet_traits[diet_traits$`Vert unk`!= 100,]#  61 spp 

# pca 
diet.pca <-princomp(clr(acomp(diet_traits)))
summary(diet.pca)

# visualizing dietary space
funspace_diet<-funspace::funspace(diet.pca, PCs = c(1, 2))
plot(funspace_diet, 
     arrows =T, arrows.length = 1.6, arrows.label.cex = 0.7, arrows.label.pos = 1.1, arrows.label.col = "black", arrows.head = 0.05,
     pnt=T, pnt.col = "black", pnt.cex=0.2, ylim = c(-3, 3))

abline(v = 0, lty = "dashed", col = "grey70")
abline(h = 0, lty = "dashed", col = "grey70")


## species probability density in dietary space (DU PrDn)

spp_scores<-as.data.frame(diet.pca$scores[,1:2])#only PCA1 and PCA2 
hv<-hypervolume_gaussian(spp_scores, samples.per.point = 2000, quantile.requested = 1)

spp_diet<-as.data.frame(hypervolume_estimate_probability(hv, spp_scores));names(spp_diet)<-"uniqueness"
spp_diet$spp<-rownames(diet.pca$scores)

#unique especies 
spp_unique<-spp_diet[spp_diet$uniqueness <= quantile(spp_diet$uniqueness, .1),]$spp

# uniqueness per community and number of unique species

comm = comm_raptors[, colnames(comm_raptors) %in% row.names(diet_traits)]

# convert to presence/absence data -- 
comm[comm < 1] = 0
comm[comm >= 1] = 1
comm = comm[rowSums(comm) > 0, ] # remove empty cells
comm= comm[, colSums(comm) > 0] # remove empty species
comm = comm[rowSums(comm) > 0, ] # remove empty cells
dim(comm)#61 spp


# (0) overall species richness of lineage

all_rich = data.frame(sites = rownames(comm), SR = rowSums(comm))
all_rich = all_rich %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

# (1) unique species richness

comm_spe = comm[, colnames(comm) %in% spp_unique]
comm_spe = comm_spe[rowSums(comm_spe) > 0, ] # remove empty cells
comm_spe = comm_spe[, colSums(comm_spe) > 0] # remove empty species

diet_rich = data.frame(sites = rownames(comm_spe), sr = rowSums(comm_spe))
diet_rich = diet_rich %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

diet_rich<-diet_rich %>% replace(is.na(.), 0)

## (2) uniqueness per community

comm_uniq<-comm[, colnames(comm) %in% spp_diet$spp]#PCA1 and PCA2
setdiff(colnames(comm_uniq), spp_diet$spp)

#replacing cells with uniqueness scores
for (i in 1:nrow(spp_diet)){
  comm_uniq[,i]<-replace(comm_uniq[,i], comm_uniq[,i] == 1, spp_diet$uniqueness[i] )
}

## dataframe with all metrics
diet_df = data.frame(sites = rownames(comm_uniq), diet_mean = rowMeans(comm_uniq))
diet_df = diet_df %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

diet_df<-dplyr::right_join(dplyr::select(coord_continents, x, y, country, continent) %>% 
                             mutate(x = as.character(x), y = as.character(y)),
                           select(diet_df, x, y, diet_mean),
                           by = c("x", "y"))
diet_df<-dplyr::full_join(diet_df, diet_rich,
                          by = c("x", "y"))
diet_df<-dplyr::left_join(diet_df, all_rich,
                          by = c("x", "y"))

diet_df[,6:7]<-diet_df[,6:7] %>% replace(is.na(.), 0)#for NA in richness

diet_df$residuals<-residuals(glm(data = diet_df, sr~SR ))

#saving all results for clade 

falcon_diet<-list(diet.pca, spp_diet, hv, spp_unique, diet_df)



## Strigiformes -------------------------------------------------------------------------------------------------

diet_traits<-raptor_traits[raptor_traits$Order == "Strigiformes", c(13:20)]
diet_traits<-diet_traits[complete.cases(diet_traits), ]#202 spp
names(diet_traits)<-c("Invertebrates", "Vert Endotherms", "Vert Ectotherms", "Fish", "Vert unk", "Scavenger", "Fruit", "Plant")

#deleting species with 100% of diet unkwnon
diet_traits<-diet_traits[diet_traits$`Vert unk`!= 100,]#  202 spp 
colSums(diet_traits)
diet_traits<-diet_traits[,-c(6:7)]#no scavengers, no fruits

# pca 
diet.pca <-princomp(clr(acomp(diet_traits)))
summary(diet.pca)

# visualizing dietary space
funspace_diet<-funspace::funspace(diet.pca, PCs = c(1, 2))
plot(funspace_diet, 
     arrows =T, arrows.length = 1.6, arrows.label.cex = 0.7, arrows.label.pos = 1.1, arrows.label.col = "black", arrows.head = 0.05,
     pnt=T, pnt.col = "black", pnt.cex=0.2, ylim = c(-3, 3))

abline(v = 0, lty = "dashed", col = "grey70")
abline(h = 0, lty = "dashed", col = "grey70")


## species probability density in dietary space (DU PrDn)

spp_scores<-as.data.frame(diet.pca$scores[,1:2])#only PCA1 and PCA2 
hv<-hypervolume_gaussian(spp_scores, samples.per.point = 2000, quantile.requested = 1)

spp_diet<-as.data.frame(hypervolume_estimate_probability(hv, spp_scores));names(spp_diet)<-"uniqueness"
spp_diet$spp<-rownames(diet.pca$scores)

#unique especies 
spp_unique<-spp_diet[spp_diet$uniqueness <= quantile(spp_diet$uniqueness, .1),]$spp

# uniqueness per community and number of unique species
comm = comm_raptors[, colnames(comm_raptors) %in% row.names(diet_traits)]

# convert to presence/absence data -- 
comm[comm < 1] = 0
comm[comm >= 1] = 1
comm = comm[rowSums(comm) > 0, ] # remove empty cells
comm= comm[, colSums(comm) > 0] # remove empty species
comm = comm[rowSums(comm) > 0, ] # remove empty cells
dim(comm)#200 spp with distributional data


# (0) overall species richness of lineage

all_rich = data.frame(sites = rownames(comm), SR = rowSums(comm))
all_rich = all_rich %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

# (1) unique species richness

comm_spe = comm[, colnames(comm) %in% spp_unique]
comm_spe = comm_spe[rowSums(comm_spe) > 0, ] # remove empty cells
comm_spe = comm_spe[, colSums(comm_spe) > 0] # remove empty species

diet_rich = data.frame(sites = rownames(comm_spe), sr = rowSums(comm_spe))
diet_rich = diet_rich %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

diet_rich<-diet_rich %>% replace(is.na(.), 0)

## (2) uniqueness per community

comm_uniq<-comm[, colnames(comm) %in% spp_diet$spp]#PCA1 and PCA2
setdiff(spp_diet$spp,colnames(comm_uniq))#"Glaucidium_mooreorum" "Otus_moheliensis" without distribution data
spp_diet<-spp_diet[(!spp_diet$spp == "Glaucidium_mooreorum" &!spp_diet$spp == "Otus_moheliensis"), ]

#replacing cells with uniqueness scores
for (i in 1:nrow(spp_diet)){
  comm_uniq[,i]<-replace(comm_uniq[,i], comm_uniq[,i] == 1, spp_diet$uniqueness[i] )
}

## dataframe with all metrics
diet_df = data.frame(sites = rownames(comm_uniq), diet_mean = rowMeans(comm_uniq))
diet_df = diet_df %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

diet_df<-dplyr::right_join(dplyr::select(coord_continents, x, y, country, continent) %>% 
                             mutate(x = as.character(x), y = as.character(y)),
                           select(diet_df, x, y, diet_mean),
                           by = c("x", "y"))
diet_df<-dplyr::full_join(diet_df, diet_rich,
                          by = c("x", "y"))
diet_df<-dplyr::left_join(diet_df, all_rich,
                          by = c("x", "y"))

diet_df[,6:7]<-diet_df[,6:7] %>% replace(is.na(.), 0)#for NA in richness

diet_df$residuals<-residuals(glm(data = diet_df, sr~SR ))

#saving all results for clade 

strig_diet<-list(diet.pca, spp_diet, hv, spp_unique, diet_df)



## plotting diet uniqueness per lineage


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

unique_plot<-accip_diet[[5]] #read same element from falcon_diet or strig_diet lists

unique_plot<-unique_plot[unique_plot$SR>0,] #delete all sites with no species
unique_plot$x<-as.numeric(unique_plot$x)
unique_plot$y<-as.numeric(unique_plot$y)

# unique richness
p_diet = p.0 + geom_tile(data = filter(unique_plot, !is.na(continent)),
                                  aes(x = x, y = y, fill = sr), inherit.aes = T) +
  labs(fill = "SR")+ ggtitle("Dietary unique species")+theme(plot.title = element_text(hjust = 0.5, size =15), plot.margin = unit(c(-1, -1, -0.5, -1), "cm"))

# unique richness -  residuals
p_diet_res = p.0 + geom_tile(data = filter(unique_plot, !is.na(continent)),
                                      aes(x = x, y = y, fill = residuals), inherit.aes = T) +
  labs(fill = expression(SR[Res]))+ ggtitle("Dietary unique species")+theme(plot.title = element_text(hjust = 0.5, size =15), plot.margin = unit(c(-1, -1, -0.5, -1), "cm"))


# community uniquenss

p_diet_comm = p.1 + geom_tile(data = filter(unique_plot, !is.na(continent)),
                                    aes(x = x, y = y, fill = diet_mean*10000), inherit.aes = T) +
  labs(fill = 'PrDn ')+ ggtitle("Diet Uniqueness") +theme(plot.title = element_text(hjust = 0.5, size =15), plot.margin = unit(c(-1, -1, -0.5, -1), "cm"))

