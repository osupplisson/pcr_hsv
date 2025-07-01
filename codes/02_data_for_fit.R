print("Cleaning ent")
rm(list = ls())
# Packages ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
required_packages <- c(
  "MatrixModels",
  "Matrix",
  "readr",
  "spdep",
  "sf",
  "janitor",
  "INLA",
  "loo",
  "ggdist",
  "tidyverse",
  "tidylog"
)


lapply(required_packages, library, character.only = TRUE)
set.seed(123)

# Spatial objects ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
spatial_objects <-
  readRDS("pcr_hsv/clean_data/input_models/spatial_objects.RDS")


# INLA graph --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
list_graph <- base::list.files(
  "pcr_hsv/clean_data/input_models/",
  all.files = T,
  full.names = T,
  pattern = ".graph"
) %>%
  tibble(folder = .) %>%
  split(as.factor(.$folder)) %>%
  imap( ~ inla.read.graph(.x$folder))


names(list_graph) <-
  str_remove(str_remove(
    str_remove(names(list_graph), "pcr_hsv/clean_data/input_models/"),
    "inla_graph_districts_"
  ), ".graph")


#
# list_graph$regions <- spatial_objects$sparse_matrix$matrix_custom_region
# list_graph$regionscities <- spatial_objects$sparse_matrix$matrix_custom_region_cities

# Data --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load the dataset
raw_data <-
  readRDS(file = "pcr_hsv/clean_data/input_models/raw_data.RDS")

raw_data <- raw_data %>%
  dplyr::select(year,
                year_week,
                insee_dep,
                sexe,
                age,
                hspcr_2,
                nature,
                n_dossier) 

spatial_objects$geometry_dep$geometry <- NULL
spatial_objects$geometry_reg$geometry <- NULL

table(raw_data$hspcr_2)


# Outcome -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
raw_data <- raw_data %>%
  mutate(
    observed_qual = case_when(
      hspcr_2 %in% c("hspn") ~ "negative",
      hspcr_2 %in% c("hspp") ~ "HSV-1/2",
      str_detect(hspcr_2, "1") ~ "HSV-1",
      str_detect(hspcr_2, "2") ~ "HSV-2"
    )
  ) %>%
  dplyr::select(year,
                year_week,
                insee_dep,
                sexe,
                age,
                hspcr_2,
                observed_qual,
                nature,
                n_dossier) %>%
  left_join(spatial_objects$geometry_dep %>% dplyr::select(insee_dep, insee_dep),
            by = "insee_dep")



# N dossier ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
id <- raw_data %>%
  distinct(n_dossier) %>%
  mutate(id = row_number(.))

raw_data <- raw_data %>%
  left_join(id, by = "n_dossier") %>%
  select(-n_dossier) %>%
  mutate(nature = factor(nature),
         females = ifelse(sexe == "2", 1, 0))

#Creating outcomes
raw_data <- raw_data %>%
  mutate(
    hsv1 = ifelse(observed_qual %in% c("HSV-1", "HSV-1/2"), 1, 0),
    hsv2 = ifelse(observed_qual %in% c("HSV-2", "HSV-1/2"), 1, 0),
    N = 1
  ) %>%
  mutate(groupsex = ifelse(females == 1, 2, 1),
         groupyear =   case_when(year == "2022" ~ 1,
                                 year == "2023" ~ 2,
                                 year == "2024" ~ 3,
                                 year == "2025" ~ 4),
         sexe = ifelse(females == 1, "Females", "Males")) %>%
  mutate(idage = age + 1) %>%
  dplyr::select(-hspcr_2,
                -id) %>%
  left_join(spatial_objects$geometry_dep, by = "insee_dep") %>%
  left_join(spatial_objects$geometry_reg,
            by = c("insee_reg"))




aggregated <- raw_data %>%
  dplyr::select(-observed_qual,
                -nature) 

aggregated <- aggregated%>% 
  group_by_at(vars(
    colnames(aggregated) %>% tibble("x" = .) %>% filter(!x %in% c("id", "hsv1", "hsv2", "N")) %>% .$x
  )) %>%
  summarise(N = sum(N),
            hsv1 = sum(hsv1),
            hsv2 = sum(hsv2)) %>%
  ungroup()


dfdate <- tibble(date = seq(as.Date("2021-01-01"), as.Date("2025-01-31"), "1 day"))  %>% 
  mutate(year = isoyear(date),
         week = isoweek(date),
         year_week = paste0(year,"_",week)) %>%
  arrange(year, week) %>% 
  distinct(year_week, .keep_all = T) %>% 
  mutate(idtime = row_number())

aggregated <- aggregated %>% 
  left_join(dfdate, by = c("year_week", "year"))

min(aggregated$week)
max(aggregated$date)
dfdate <- dfdate %>% 
  filter(idtime >= min(aggregated$idtime)) %>% 
  filter(idtime <= max(aggregated$idtime)) %>% 
  mutate(idtime = idtime - min(idtime) + 1)
aggregated <- aggregated %>% 
  mutate(idtime = idtime - min(idtime)+1)

# Save --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
save.image("pcr_hsv/clean_data/input_models/data_for_fit.rda")

