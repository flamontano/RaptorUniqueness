## Protection status of raptor uniqueness ##

library(dplyr)
library(tidyverse)

devtools::install_github("tidyverse/dplyr")
## all protected areas
## strict protected areas (categories I-II UICN)

pa_results<-read.csv("Data/raptors_pa_areas.csv", header=T)

known <- data.frame(
  x = c(1000, 250000),
  y = c(100, 10)
)#following Morelli et al. 2021

pa_results$target <- ifelse(pa_results$sp_distri_area_km2 <1000, 100, ifelse(pa_results$sp_distri_area_km2 >  250000, 10, approx(known$x, known$y, pa_results$sp_distri_area_km2)$y))

# species protection in PA

pa_results$protected<-ifelse(pa_results$target < pa_results$pa_overlap_distri_perc*100, 1, 0)
pa_results$strict_protected<-ifelse(pa_results$target < pa_results$pa_strict_overlap_distri_perc *100, 1, 0)

colSums(pa_results[,12:13])# 176 protected, only 19 strictly protected
colMeans(pa_results[12:13] != 0)# percentage protected


# separate protection per lineage
raptor_traits$sp<-rownames(raptor_traits)

protected_accip<-pa_results %>% left_join(raptor_traits, by = dplyr::join_by(sp == sp))%>% select(sp, protected, strict_protected, Order)%>%filter(Order == "Accipitriformes")
protected_falcon<-pa_results %>% left_join(raptor_traits, by = join_by(sp == sp))%>% select(sp, protected, strict_protected, Order)%>%filter(Order == "Falconiformes")
protected_strig<-pa_results %>% left_join(raptor_traits, by = join_by(sp == sp))%>% select(sp, protected, strict_protected, Order)%>%filter(Order == "Strigiformes")

## I. Summary statistics per lineage
# Accipitridae
colSums(protected_accip[,2:3])# 94 protected, 7 strictly protected
colMeans(protected_accip[2:3] != 0) # percentage protected
#Falconidae
colSums(protected_falcon[,2:3])# 28 protected, 5 strictly protected
colMeans(protected_falcon[2:3] != 0) # percentage protected
#Strigidae
colSums(protected_strig[,2:3])# 54 protected, 7 strictly protected
colMeans(protected_strig[2:3] != 0) # percentage protected

#II Are protected more unique than non protected?

## Accipitriformes

diet<-accip_diet[[2]]
diet$spp<-sub(" ", "_", diet$spp)
names(diet)[names(diet) == 'uniqueness'] <- 'unique'
morpho <-accip_morpho[[2]]
ed <- accip_evo[[1]]

df_accip<-protected_accip%>%left_join(diet, by = join_by(sp == spp))%>% 
  mutate(dietuniq = ifelse(unique<= quantile(unique, .1,na.rm=TRUE),1,0))%>%
  left_join(morpho,by = join_by(sp == spp))%>% mutate(morphouniq = ifelse(uniqueness<= quantile(uniqueness, .1),1,0 ))%>%
  left_join(ed,by = join_by(sp == Species))%>% mutate(evouniq = ifelse(w>= quantile(w, .9, na.rm=TRUE),1,0 ))

# number and proportion of protected and unprotected spp
df = df_accip

# number of specialist species protected
sum(df$protected*df$dietuniq, na.rm = TRUE)#protected special species
sum(df$protected*df$dietuniq, na.rm = TRUE)/sum(df$dietuniq, na.rm = TRUE)#proportion of the unique species protected
# number of unique 12 species protected
sum(df$protected*df$morphouniq, na.rm = TRUE)#protected special species
sum(df$protected*df$morphouniq, na.rm = TRUE)/sum(df$morphouniq, na.rm = TRUE)#proportion of the unique species protected
# number of phylo unique species protected
sum(df$protected*df$evouniq, na.rm = TRUE)#protected special species
sum(df$protected*df$evouniq, na.rm = TRUE)/sum(df$evouniq, na.rm = TRUE)#proportion of the unique species protected

# number of specialist species strictly protected
sum(df$strict_protected*df$dietuniq, na.rm = TRUE)#protected special species
sum(df$strict_protected*df$dietuniq, na.rm = TRUE)/sum(df$dietuniq, na.rm = TRUE)#proportion of the unique species protected
# number of unique 12 species protected
sum(df$strict_protected*df$morphouniq, na.rm = TRUE)#protected special species
sum(df$strict_protected*df$morphouniq, na.rm = TRUE)/sum(df$morphouniq, na.rm = TRUE)#proportion of the unique species protected
# number of phylo unique species protected
sum(df$strict_protected*df$evouniq, na.rm = TRUE)#protected special species
sum(df$strict_protected*df$evouniq, na.rm = TRUE)/sum(df$evouniq, na.rm = TRUE)#proportion of the unique species protected


df$protected<-as.factor(df$protected)
prot <- subset(df, protected == "1")
noprot <- subset(df, protected == "0")

#diet
var.test(prot$unique, noprot$unique)
t.test(prot$unique, noprot$unique, var.equal = T)
#functional uniqueness
var.test(prot$uniqueness, noprot$uniqueness)
t.test(prot$uniqueness, noprot$uniqueness, var.equal = T)
# phylo uniqueness
var.test(prot$w, noprot$w)
t.test(prot$w, noprot$w, var.equal = T)

## Falconiformes

diet<-falcon_diet[[2]]
diet$spp<-sub(" ", "_", diet$spp)
names(diet)[names(diet) == 'uniqueness'] <- 'unique'
morpho <-falcon_morpho[[2]]
ed <- falcon_evo[[1]]

df_falcon<-protected_falcon%>%left_join(diet, by = join_by(sp == spp))%>% 
  mutate(dietuniq = ifelse(unique<= quantile(unique, .1,na.rm=TRUE),1,0))%>%
  left_join(morpho,by = join_by(sp == spp))%>% mutate(morphouniq = ifelse(uniqueness<= quantile(uniqueness, .1),1,0 ))%>%
  left_join(ed,by = join_by(sp == Species))%>% mutate(evouniq = ifelse(w>= quantile(w, .9, na.rm=TRUE),1,0 ))

# number and proportion of protected and unprotected spp
df = df_falcon

# number of specialist species protected
sum(df$protected*df$dietuniq, na.rm = TRUE)#protected special species
sum(df$protected*df$dietuniq, na.rm = TRUE)/sum(df$dietuniq, na.rm = TRUE)#proportion of the unique species protected
# number of unique 12 species protected
sum(df$protected*df$morphouniq, na.rm = TRUE)#protected special species
sum(df$protected*df$morphouniq, na.rm = TRUE)/sum(df$morphouniq, na.rm = TRUE)#proportion of the unique species protected
# number of phylo unique species protected
sum(df$protected*df$evouniq, na.rm = TRUE)#protected special species
sum(df$protected*df$evouniq, na.rm = TRUE)/sum(df$evouniq, na.rm = TRUE)#proportion of the unique species protected

# number of specialist species strictly protected
sum(df$strict_protected*df$dietuniq, na.rm = TRUE)#protected special species
sum(df$strict_protected*df$dietuniq, na.rm = TRUE)/sum(df$dietuniq, na.rm = TRUE)#proportion of the unique species protected
# number of unique 12 species protected
sum(df$strict_protected*df$morphouniq, na.rm = TRUE)#protected special species
sum(df$strict_protected*df$morphouniq, na.rm = TRUE)/sum(df$morphouniq, na.rm = TRUE)#proportion of the unique species protected
# number of phylo unique species protected
sum(df$strict_protected*df$evouniq, na.rm = TRUE)#protected special species
sum(df$strict_protected*df$evouniq, na.rm = TRUE)/sum(df$evouniq, na.rm = TRUE)#proportion of the unique species protected


df$protected<-as.factor(df$protected)
prot <- subset(df, protected == "1")
noprot <- subset(df, protected == "0")

#diet
var.test(prot$unique, noprot$unique)
t.test(prot$unique, noprot$unique, var.equal = T)
#functional uniqueness
var.test(prot$uniqueness, noprot$uniqueness)
t.test(prot$uniqueness, noprot$uniqueness, var.equal = T)
# phylo uniqueness
var.test(prot$w, noprot$w)
t.test(prot$w, noprot$w, var.equal = T)


## Strigiformes

diet<-strig_diet[[2]]
diet$spp<-sub(" ", "_", diet$spp)
names(diet)[names(diet) == 'uniqueness'] <- 'unique'
morpho <-strig_morpho[[2]]
ed <- strig_evo[[1]]

df_strig<-protected_strig%>%left_join(diet, by = join_by(sp == spp))%>% 
  mutate(dietuniq = ifelse(unique<= quantile(unique, .1,na.rm=TRUE),1,0))%>%
  left_join(morpho,by = join_by(sp == spp))%>% mutate(morphouniq = ifelse(uniqueness<= quantile(uniqueness, .1),1,0 ))%>%
  left_join(ed,by = join_by(sp == Species))%>% mutate(evouniq = ifelse(w>= quantile(w, .9, na.rm=TRUE),1,0 ))

# number and proportion of protected and unprotected spp
df = df_strig

# number of specialist species protected
sum(df$protected*df$dietuniq, na.rm = TRUE)#protected special species
sum(df$protected*df$dietuniq, na.rm = TRUE)/sum(df$dietuniq, na.rm = TRUE)#proportion of the unique species protected
# number of unique 12 species protected
sum(df$protected*df$morphouniq, na.rm = TRUE)#protected special species
sum(df$protected*df$morphouniq, na.rm = TRUE)/sum(df$morphouniq, na.rm = TRUE)#proportion of the unique species protected
# number of phylo unique species protected
sum(df$protected*df$evouniq, na.rm = TRUE)#protected special species
sum(df$protected*df$evouniq, na.rm = TRUE)/sum(df$evouniq, na.rm = TRUE)#proportion of the unique species protected

# number of specialist species strictly protected
sum(df$strict_protected*df$dietuniq, na.rm = TRUE)#protected special species
sum(df$strict_protected*df$dietuniq, na.rm = TRUE)/sum(df$dietuniq, na.rm = TRUE)#proportion of the unique species protected
# number of unique 12 species protected
sum(df$strict_protected*df$morphouniq, na.rm = TRUE)#protected special species
sum(df$strict_protected*df$morphouniq, na.rm = TRUE)/sum(df$morphouniq, na.rm = TRUE)#proportion of the unique species protected
# number of phylo unique species protected
sum(df$strict_protected*df$evouniq, na.rm = TRUE)#protected special species
sum(df$strict_protected*df$evouniq, na.rm = TRUE)/sum(df$evouniq, na.rm = TRUE)#proportion of the unique species protected


df$protected<-as.factor(df$protected)
prot <- subset(df, protected == "1")
noprot <- subset(df, protected == "0")

#diet
var.test(prot$unique, noprot$unique)
t.test(prot$unique, noprot$unique, var.equal = T)
#functional uniqueness
var.test(prot$uniqueness, noprot$uniqueness)
t.test(prot$uniqueness, noprot$uniqueness, var.equal = T)
# phylo uniqueness
var.test(prot$w, noprot$w)
t.test(prot$w, noprot$w, var.equal = T)




