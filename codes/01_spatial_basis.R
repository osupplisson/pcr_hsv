#This code produces the adjacency matrices further used for the BYM2 structured prior

# Packages ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(spdep)
library(janitor)
library(tidyverse)
library(tidylog)
library(dbscan)
# Importing shapefile from France -----------------------------------------------------------------------------------------------------------------------------------------------------------------
geometry_dep <-
  read_sf(
    "global_raw_data/data_for_cartography/1_DONNEES_LIVRAISON_2022-08-30/ADE_3-1_SHP_WGS84G_FRA/DEPARTEMENT.shp"
  ) %>%
  clean_names()



# Extracting Corsica ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
geometry_dep <- geometry_dep %>% 
  filter(!insee_dep %in% c("971", "972", "973", "974", "976")) 


# Extracting Corsica ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
corsica <- geometry_dep %>% 
  filter(insee_dep %in% c("2A", "2B"))

#Union of Corsican distritcs
union_corsica <- st_union(corsica$geometry[1],
                          corsica$geometry[2])

new_corsica <- corsica %>%
  filter(row_number() == 1) %>%
  mutate(id = 'XXX',
         nom_m = "CORSE",
         nom = 'Corse',
         insee_dep = "20",
         insee_reg = "94",
         geometry = union_corsica)


geometry_dep <- geometry_dep %>% 
  filter(!insee_dep %in% c("2A", "2B")) %>% 
  add_row(new_corsica, .after = 19)


rm("corsica")
rm("new_corsica")


# Adding region names -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
geometry_reg <-
  read_sf(
    "global_raw_data/data_for_cartography/1_DONNEES_LIVRAISON_2022-08-30/ADE_3-1_SHP_WGS84G_FRA/REGION.shp"
  ) %>%
  clean_names() %>% 
  dplyr::select(nom_region = nom,
                insee_reg = insee_reg) %>% 
  filter(!nom_region %in% c("Guadeloupe", "Martinique", "Guyane", "La Réunion", "Mayotte")) 


library(rmapshaper)
geometry_dep_bis <- geometry_dep
geometry_dep_bis$geometry  <- ms_simplify(geometry_dep_bis$geometry,
                                          keep = 0.05,
                                          keep_shapes = T)
saveRDS(geometry_dep_bis, "pcr_hsv/clean_data/input_models/geometry_fr.RDS")
rm("geometry_dep_bis")
# Extracting centroids
centroids <- st_coordinates(st_centroid(st_geometry(geometry_dep$geometry)))


# Neighboorhood matrix for CAR models---------------------------------------------------------------------------------------------------------------
# https:/r-spatial.github.io/spdep/articles/nb.html
# Creating queen contiguity
queen_contiguity_nbobject <- poly2nb(as(geometry_dep, "Spatial"),
                                     queen = TRUE
)

# SPDED to INLA graph
spdep::nb2INLA(
  "pcr_hsv/clean_data/input_models/inla_graph_districts_queen.graph",
  queen_contiguity_nbobject
)


# Delaunay triangulation neighbours
delaunay_nbobject <- tri2nb(centroids)
spdep::nb2INLA(
  "pcr_hsv/clean_data/input_models/inla_graph_districts_delaunay.graph",
  delaunay_nbobject
)

#SOI neighbours are symmetric by design
soi_nbobject <- graph2nb(soi.graph(delaunay_nbobject, centroids))
spdep::nb2INLA(
  "pcr_hsv/clean_data/input_models/inla_graph_districts_soi.graph",
  soi_nbobject
)

#Gabriel
gab_nbobject <- graph2nb(gabrielneigh(centroids), sym=TRUE)
spdep::nb2INLA(
  "pcr_hsv/clean_data/input_models/inla_graph_districts_gabriel.graph",
  gab_nbobject
)

#Relative graph
relative_nbobject <- graph2nb(relativeneigh(centroids), sym=TRUE)
spdep::nb2INLA(
  "pcr_hsv/clean_data/input_models/inla_graph_districts_relative.graph",
  relative_nbobject
)


# Neighbors ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
nb_list <- tibble(knn = seq(1,25)) %>% 
  split(as.factor(.$knn)) %>%
  imap(~knn2nb(knearneigh(centroids, k=.x), row.names=geometry_dep$id)) %>% 
  imap(~make.sym.nb(.x))
names(nb_list) <- paste0("nb", names(nb_list))


nb_list <- nb_list %>% 
  imap(function(.x, .y){
    list <- nblag(.x,3)
    names(list) <-paste0(.y,'_', seq(1,3))
    list}) %>% 
  flatten() 

nb_list %>% 
  imap(~spdep::nb2INLA(
    paste0("pcr_hsv/clean_data/input_models/inla_graph_districts_",
           .y,".graph"),
    .x
  ))


# Queen lagged ------------------------------------------------------------
queenlag2 <- spdep::nblag(queen_contiguity_nbobject = nb, maxlag = 2)
queenlag3 <- spdep::nblag(queen_contiguity_nbobject = nb, maxlag = 3)




# Saving -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
spatial_objects <- list(
  geometry_dep = geometry_dep %>% mutate(idspace = row_number()),
  geometry_reg = geometry_reg %>% mutate(idspacereg = row_number()),
  centroids = centroids,
  nb = list('queen' = queen_contiguity_nbobject,
            'delaunay' = delaunay_nbobject,
            'soi' = soi_nbobject,
            'gab_nbobject' = gab_nbobject,
            'relative_nbobject' = relative_nbobject,
            'nb_list' = nb_list))

saveRDS(spatial_objects,
        "pcr_hsv/clean_data/input_models/spatial_objects.RDS")


