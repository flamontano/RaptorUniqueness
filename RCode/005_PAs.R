library(sf)
library(tidyverse)

if(!file.exists("global_protected_area_simplified.rds")){
  # download from http://wcmc.io/wdpa_current_release
  sf::st_layers("WDPA_WDOECM_May2023_Public/WDPA_WDOECM_May2023_Public.gdb/")
  
  x = sf::st_read("WDPA_WDOECM_May2023_Public/WDPA_WDOECM_May2023_Public.gdb/", 
                  layer = "WDPA_WDOECM_poly_May2023")
  plot(st_geometry(x))
  any(is.na(x$DESIG_ENG))
  
  # st_is_valid(x)
  names(x)
  table(x$STATUS)
  
  x = filter(x, DESIG_ENG != 'UNESCO-MAB Biosphere Reserve',
             !STATUS %in% c('Proposed', 'Not Reported'))
  # x2 = st_simplify(x, preserveTopology = TRUE)
  
  x = st_transform(x, crs = "ESRI:54002") # World_Equidistant_Cylindrical
  
  x2 = x[1:100,]
  x2
  plot(st_geometry(x2))
  
  st_s = function(x){
    xx = try(sf::st_simplify(x, preserveTopology = T, dTolerance = 1000))
    if(inherits(xx, "try-error")) {
      return(x)
    } else {
      return(xx)
    }
  }
  
  x_simple = mutate(x, Shape = st_s(Shape))
  
  # x_simple2 = st_union(st_make_valid(x_simple)) # not a good idea
  # saveRDS(x_simple2, file = "global_protected_area_simplified_unioned.rds")
  
  saveRDS(x_simple, file = "global_protected_area_simplified.rds")
} else {
  x_simple = readRDS("global_protected_area_simplified.rds")
}


# remove areas that are < 1 sq km
x_simple = filter(x_simple, REP_AREA >= 1)


sp_list = readRDS("/Users/dli/Library/CloudStorage/Box-Box/Collaboration/Flavia_Montano-Centellas/raptors_protect_area/species.rds")

sp_list = unlist(sp_list)
# Strix_hadorami does not match with any map, it is a least concern species, ignore it?
sp_list = setdiff(sp_list, "Strix_hadorami")

sp_list_map = list.files("/Users/dli/Library/CloudStorage/Box-Box/UFL/world_birds/Raptors_map", pattern = "*.shp$", full.names = T)

raptors_pa_areas = tibble(sp = sp_list, sp_distri_area_km2 = NA,
                          pa_areas_intersect_range_km2 = NA,
                          pa_strict_areas_intersect_range_km2 = NA, 
                          pa_areas_overlap_range_km2 = NA,
                          pa_strict_areas_overlap_range_km2 = NA)

for(i in 1:length(sp_list)){
  cat(i, "out of 548 ", "\n")
  sp_match = grep(paste0(sp_list[i], "_"), sp_list_map, value = T)
  # read distribution map
  sp_distri = st_read(sp_match) |> 
    st_transform(crs = st_crs(x_simple)) |> 
    st_make_valid()
  sp_distri = sp_distri[0] # remove other columns
  sp_distri_area_km2 = st_area(sp_distri) |> units::set_units(km^2) |> 
    as.vector() |> sum()
  raptors_pa_areas$sp_distri_area_km2[i] = sp_distri_area_km2
  # protected areas intersect
  pa_intersect_range = try(st_make_valid(x_simple[st_intersects(sp_distri, x_simple)[[1]], ]))
  
  if(!inherits(pa_intersect_range, "try-error")){
    pa_strict_intersect_range = filter(pa_intersect_range, IUCN_CAT %in% c("Ia", "Ib", "II"))
    
    # st_union
    pa_intersect_range = st_union(pa_intersect_range)
    pa_strict_intersect_range = st_union(pa_strict_intersect_range)
    
    raptors_pa_areas$pa_areas_intersect_range_km2[i] = st_area(pa_intersect_range) |> 
      units::set_units(km^2) |> 
      as.vector() |> sum()
    
    raptors_pa_areas$pa_strict_areas_intersect_range_km2[i] = st_area(pa_strict_intersect_range) |> 
      units::set_units(km^2) |> 
      as.vector() |> sum()
    
    # intersection?? i.e., overlapping areas? It makes more sense to actually have the whole PA, right?
    # add here just in case
    pa_overlapping = try(st_union(st_intersection(pa_intersect_range, sp_distri)))
    pa_overlapping_2 = try(st_union(st_intersection(pa_strict_intersect_range, sp_distri)))
    
    if(!inherits(pa_overlapping, "try-error")){
      raptors_pa_areas$pa_areas_overlap_range_km2[i] = st_area(pa_overlapping) |> 
        units::set_units(km^2) |> as.vector() |> sum()
      
      raptors_pa_areas$pa_strict_areas_overlap_range_km2[i] = st_area(pa_overlapping_2) |> 
        units::set_units(km^2) |> as.vector() |> sum()
    }
  }
}

saveRDS(raptors_pa_areas, "raptors_pa_areas2.rds")

raptors_pa_areas = readRDS("raptors_pa_areas.rds")

names(raptors_pa_areas)

raptors_pa_areas = mutate(raptors_pa_areas, pa_intersect_distri_perc = pa_areas_intersect_range_km2 / sp_distri_area_km2,
                          pa_strict_intersect_distri_perc = pa_strict_areas_intersect_range_km2 / sp_distri_area_km2,
                          pa_overlap_distri_perc = pa_areas_overlap_range_km2 / sp_distri_area_km2,
                          pa_strict_overlap_distri_perc = pa_strict_areas_overlap_range_km2 / sp_distri_area_km2
)

#write_csv(raptors_pa_areas, "Data/raptors_pa_areas.csv")

