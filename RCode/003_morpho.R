## MORPHOLOGICAL UNIQUENESS ##

library(ecodist)
library(vegan)
library(hypervolume)
library(compositions)
library(factoextra)
library(ggplot2)
library(plyr)
library(dplyr)


# Accipitriformes ## -----------------------------------------------------------------------------------

morpho_traits<-raptor_traits[raptor_traits$Order == "Accipitriformes", c(6:12)]
morpho_traits<-morpho_traits[complete.cases(morpho_traits), ]#250 spp
morpho_traits<-scale(log(morpho_traits))
names(morpho_traits)<-c("Bill length", "Bill width", "Bill depth", "Tarsus", "Wing", "Secondary length", "Tail")

# pca 
morpho.pca <- princomp(morpho_traits, cor=T)#not including body mass
summary(morpho.pca)

# visualizing dietary space
funspace_morpho<-funspace::funspace(morpho.pca, PCs = c(1, 2))
plot(funspace_morpho, 
     arrows =T, arrows.length = 1.6, arrows.label.cex = 0.7, arrows.label.pos = 1.1, arrows.label.col = "black", arrows.head = 0.05,
     pnt=T, pnt.col = "black", pnt.cex=0.2)

abline(v = 0, lty = "dashed", col = "grey70")
abline(h = 0, lty = "dashed", col = "grey70")


## species probability density in morphological space (MU PrDn)

spp_scores<-as.data.frame(morpho.pca$scores[,1:2])#only PCA1 and PCA2 
hv<-hypervolume_gaussian(spp_scores, samples.per.point = 2000, quantile.requested = 1)

spp_morpho<-as.data.frame(hypervolume_estimate_probability(hv, spp_scores));names(spp_morpho)<-"uniqueness"
spp_morpho$spp<-rownames(morpho.pca$scores)

#unique especies 
spp_unique<-spp_morpho[spp_morpho$uniqueness <= quantile(spp_morpho$uniqueness, .1),]$spp

# uniqueness per community and number of unique species

comm = comm_raptors[, colnames(comm_raptors) %in% row.names(morpho_traits)]

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

# (1) unique species richness

# -- PCA! and PCA2
comm_spe = comm[, colnames(comm) %in% spp_unique]
comm_spe = comm_spe[rowSums(comm_spe) > 0, ] # remove empty cells
comm_spe = comm_spe[, colSums(comm_spe) > 0] # remove empty species

morpho_rich = data.frame(sites = rownames(comm_spe), sr = rowSums(comm_spe))
morpho_rich = morpho_rich %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

morpho_rich<-morpho_rich %>% replace(is.na(.), 0)

## (2) uniqueness per community

comm_uniq<-comm[, colnames(comm) %in% spp_morpho$spp]#PCA1 and PCA2
setdiff(colnames(comm_uniq), spp_morpho$spp)

#replacing cells with uniqueness scores
for (i in 1:nrow(spp_morpho)){
  comm_uniq[,i]<-replace(comm_uniq[,i], comm_uniq[,i] == 1, spp_morpho$uniqueness[i] )
}

## dataframe with all metrics
morpho_df = data.frame(sites = rownames(comm_uniq), morpho_mean = rowMeans(comm_uniq))
morpho_df = morpho_df %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

morpho_df<-dplyr::right_join(dplyr::select(coord_continents, x, y, country, continent) %>% 
                             mutate(x = as.character(x), y = as.character(y)),
                           select(morpho_df, x, y, morpho_mean),
                           by = c("x", "y"))
morpho_df<-dplyr::full_join(morpho_df, morpho_rich,
                          by = c("x", "y"))
morpho_df<-dplyr::left_join(morpho_df, all_rich,
                          by = c("x", "y"))

morpho_df[,6:7]<-morpho_df[,6:7] %>% replace(is.na(.), 0)#for NA in richness

morpho_df$residuals<-residuals(glm(data = morpho_df, sr~SR ))

#saving all results for clade 

accip_morpho<-list(morpho.pca, spp_morpho, hv, spp_unique, morpho_df)


# Falconiformes ## -----------------------------------------------------------------------------------

morpho_traits<-raptor_traits[raptor_traits$Order == "Falconiformes", c(6:12)]
morpho_traits<-morpho_traits[complete.cases(morpho_traits), ]#64 spp
morpho_traits<-scale(log(morpho_traits))
names(morpho_traits)<-c("Bill length", "Bill width", "Bill depth", "Tarsus", "Wing", "Secondary length", "Tail")

# pca 
morpho.pca <- princomp(morpho_traits, cor=T)#not including body mass
summary(morpho.pca)

# visualizing dietary space
funspace_morpho<-funspace::funspace(morpho.pca, PCs = c(1, 2))
plot(funspace_morpho, 
     arrows =T, arrows.length = 1.6, arrows.label.cex = 0.7, arrows.label.pos = 1.1, arrows.label.col = "black", arrows.head = 0.05,
     pnt=T, pnt.col = "black", pnt.cex=0.2)

abline(v = 0, lty = "dashed", col = "grey70")
abline(h = 0, lty = "dashed", col = "grey70")


## species probability density in morphological space (MU PrDn)

spp_scores<-as.data.frame(morpho.pca$scores[,1:2])#only PCA1 and PCA2 
hv<-hypervolume_gaussian(spp_scores, samples.per.point = 2000, quantile.requested = 1)

spp_morpho<-as.data.frame(hypervolume_estimate_probability(hv, spp_scores));names(spp_morpho)<-"uniqueness"
spp_morpho$spp<-rownames(morpho.pca$scores)

#unique especies 
spp_unique<-spp_morpho[spp_morpho$uniqueness <= quantile(spp_morpho$uniqueness, .1),]$spp

# uniqueness per community and number of unique species

comm = comm_raptors[, colnames(comm_raptors) %in% row.names(morpho_traits)]

# convert to presence/absence data -- 
comm[comm < 1] = 0
comm[comm >= 1] = 1
comm = comm[rowSums(comm) > 0, ] # remove empty cells
comm= comm[, colSums(comm) > 0] # remove empty species
comm = comm[rowSums(comm) > 0, ] # remove empty cells
dim(comm)#64 spp

# (0) overall species richness of lineage

all_rich = data.frame(sites = rownames(comm), SR = rowSums(comm))
all_rich = all_rich %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

# (1) unique species richness

# -- PCA! and PCA2
comm_spe = comm[, colnames(comm) %in% spp_unique]
comm_spe = comm_spe[rowSums(comm_spe) > 0, ] # remove empty cells
comm_spe = comm_spe[, colSums(comm_spe) > 0] # remove empty species

morpho_rich = data.frame(sites = rownames(comm_spe), sr = rowSums(comm_spe))
morpho_rich = morpho_rich %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

morpho_rich<-morpho_rich %>% replace(is.na(.), 0)

## (2) uniqueness per community

comm_uniq<-comm[, colnames(comm) %in% spp_morpho$spp]#PCA1 and PCA2
setdiff(colnames(comm_uniq), spp_morpho$spp)

#replacing cells with uniqueness scores
for (i in 1:nrow(spp_morpho)){
  comm_uniq[,i]<-replace(comm_uniq[,i], comm_uniq[,i] == 1, spp_morpho$uniqueness[i] )
}

## dataframe with all metrics
morpho_df = data.frame(sites = rownames(comm_uniq), morpho_mean = rowMeans(comm_uniq))
morpho_df = morpho_df %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

morpho_df<-dplyr::right_join(dplyr::select(coord_continents, x, y, country, continent) %>% 
                               mutate(x = as.character(x), y = as.character(y)),
                             select(morpho_df, x, y, morpho_mean),
                             by = c("x", "y"))
morpho_df<-dplyr::full_join(morpho_df, morpho_rich,
                            by = c("x", "y"))
morpho_df<-dplyr::left_join(morpho_df, all_rich,
                            by = c("x", "y"))

morpho_df[,6:7]<-morpho_df[,6:7] %>% replace(is.na(.), 0)#for NA in richness

morpho_df$residuals<-residuals(glm(data = morpho_df, sr~SR ))

#saving all results for clade 

falcon_morpho<-list(morpho.pca, spp_morpho, hv, spp_unique, morpho_df)


# Strigiformes ## -----------------------------------------------------------------------------------

morpho_traits<-raptor_traits[raptor_traits$Order == "Strigiformes", c(6:12)]
morpho_traits<-morpho_traits[complete.cases(morpho_traits), ]#235 spp
morpho_traits<-scale(log(morpho_traits))
names(morpho_traits)<-c("Bill length", "Bill width", "Bill depth", "Tarsus", "Wing", "Secondary length", "Tail")

# pca 
morpho.pca <- princomp(morpho_traits, cor=T)#not including body mass
summary(morpho.pca)

# visualizing dietary space
funspace_morpho<-funspace::funspace(morpho.pca, PCs = c(1, 2))
plot(funspace_morpho, 
     arrows =T, arrows.length = 1.6, arrows.label.cex = 0.7, arrows.label.pos = 1.1, arrows.label.col = "black", arrows.head = 0.05,
     pnt=T, pnt.col = "black", pnt.cex=0.2)

abline(v = 0, lty = "dashed", col = "grey70")
abline(h = 0, lty = "dashed", col = "grey70")


## species probability density in morphological space (MU PrDn)

spp_scores<-as.data.frame(morpho.pca$scores[,1:2])#only PCA1 and PCA2 
hv<-hypervolume_gaussian(spp_scores, samples.per.point = 2000, quantile.requested = 1)

spp_morpho<-as.data.frame(hypervolume_estimate_probability(hv, spp_scores));names(spp_morpho)<-"uniqueness"
spp_morpho$spp<-rownames(morpho.pca$scores)

#unique especies 
spp_unique<-spp_morpho[spp_morpho$uniqueness <= quantile(spp_morpho$uniqueness, .1),]$spp

# uniqueness per community and number of unique species

comm = comm_raptors[, colnames(comm_raptors) %in% row.names(morpho_traits)]

# convert to presence/absence data -- 
comm[comm < 1] = 0
comm[comm >= 1] = 1
comm = comm[rowSums(comm) > 0, ] # remove empty cells
comm= comm[, colSums(comm) > 0] # remove empty species
comm = comm[rowSums(comm) > 0, ] # remove empty cells
dim(comm)#231 spp

# (0) overall species richness of lineage

all_rich = data.frame(sites = rownames(comm), SR = rowSums(comm))
all_rich = all_rich %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

# (1) unique species richness

# -- PCA! and PCA2
comm_spe = comm[, colnames(comm) %in% spp_unique]
comm_spe = comm_spe[rowSums(comm_spe) > 0, ] # remove empty cells
comm_spe = comm_spe[, colSums(comm_spe) > 0] # remove empty species

morpho_rich = data.frame(sites = rownames(comm_spe), sr = rowSums(comm_spe))
morpho_rich = morpho_rich %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

morpho_rich<-morpho_rich %>% replace(is.na(.), 0)

## (2) uniqueness per community

comm_uniq<-comm[, colnames(comm) %in% spp_morpho$spp]#PCA1 and PCA2
setdiff(spp_morpho$spp, colnames(comm_uniq))

spp_morpho1<-(spp_morpho[!(spp_morpho$spp == "Glaucidium_mooreorum")& 
                                 !(spp_morpho$spp =="Otus_feae")&
                                 !(spp_morpho$spp =="Otus_moheliensis")& 
                                 !(spp_morpho$spp =="Strix_hadorami") , ])#species with traits but no distribution data


#replacing cells with uniqueness scores
for (i in 1:nrow(spp_morpho1)){
  comm_uniq[,i]<-replace(comm_uniq[,i], comm_uniq[,i] == 1, spp_morpho1$uniqueness[i] )
}

## dataframe with all metrics
morpho_df = data.frame(sites = rownames(comm_uniq), morpho_mean = rowMeans(comm_uniq))
morpho_df = morpho_df %>%
  tidyr::separate(sites, c("x", "y"), sep = "_")

morpho_df<-dplyr::right_join(dplyr::select(coord_continents, x, y, country, continent) %>% 
                               mutate(x = as.character(x), y = as.character(y)),
                             select(morpho_df, x, y, morpho_mean),
                             by = c("x", "y"))
morpho_df<-dplyr::full_join(morpho_df, morpho_rich,
                            by = c("x", "y"))
morpho_df<-dplyr::left_join(morpho_df, all_rich,
                            by = c("x", "y"))

morpho_df[,6:7]<-morpho_df[,6:7] %>% replace(is.na(.), 0)#for NA in richness

morpho_df$residuals<-residuals(glm(data = morpho_df, sr~SR ))

#saving all results for clade 

strig_morpho<-list(morpho.pca, spp_morpho, hv, spp_unique, morpho_df)



## plotting morphological uniqueness per lineage

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

unique_plot<-accip_morpho[[5]] #read same element from falcon_diet or strig_diet lists

unique_plot<-unique_plot[unique_plot$SR>0,] #delete all sites with no species
unique_plot$x<-as.numeric(unique_plot$x)
unique_plot$y<-as.numeric(unique_plot$y)

# unique richness
p_morpho = p.0 + geom_tile(data = filter(unique_plot, !is.na(continent)),
                         aes(x = x, y = y, fill = sr), inherit.aes = T) +
  labs(fill = "SR")+ ggtitle("Morphologically unique species")+theme(plot.title = element_text(hjust = 0.5, size =15), plot.margin = unit(c(-1, -1, -0.5, -1), "cm"))

# unique richness -  residuals
p_morpho_res = p.0 + geom_tile(data = filter(unique_plot, !is.na(continent)),
                             aes(x = x, y = y, fill = residuals), inherit.aes = T) +
  labs(fill = expression(SR[Res]))+ ggtitle("Morphologically unique species")+theme(plot.title = element_text(hjust = 0.5, size =15), plot.margin = unit(c(-1, -1, -0.5, -1), "cm"))


# community uniquenss

p_morphot_comm = p.1 + geom_tile(data = filter(unique_plot, !is.na(continent)),
                              aes(x = x, y = y, fill = morpho_mean*100000), inherit.aes = T) +
  labs(fill = 'PrDn ')+ ggtitle("Morphological Uniqueness") +theme(plot.title = element_text(hjust = 0.5, size =15), plot.margin = unit(c(-1, -1, -0.5, -1), "cm"))

