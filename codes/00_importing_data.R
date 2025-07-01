#Cleaning R objects ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rm(list = ls())

# Packages ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(readxl)
library(xlsx)
library(janitor)
library(tidyverse)
library(tidylog)

# Importing data ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
file <- read_excel("pcr_hsv\\raw_data\\extraction_janvier_2025.xlsx") %>%
  clean_names()


file <- file %>% 
  mutate(dateprelev = as.Date(dateprelev)) %>% 
  mutate(year = isoyear(dateprelev),
         week = isoweek(dateprelev)) %>% 
  mutate(year_week = paste0(year,"_",week)) 

table(file$year)

file <- file %>% 
  filter(year %in% c("2022","2023", "2024"))

min(file$dateprelev)
max(file$dateprelev)

#DISTRICTS
file <- file %>% 
  mutate(insee_dep = str_sub(code_postal_patient, start = 1, end = 2))


# Removing NA -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Keeping only selected nature
table(file$nature,file$year)
table(file$nature)

nature_to_keep <- c("ANAL",
                    "GENIT",
                    "PU",
                    "PV",
                    "VULV")

#Start of the flowchart
file <- file %>% 
  filter(nature %in% nature_to_keep) 


file <- file %>% 
  filter(!is.na(as.numeric(age))) %>% 
  filter(as.numeric(age) >= 0) %>% 
  filter(!is.na(sexe)) %>%
  filter(!is.na(insee_dep)) %>% 
  filter(!insee_dep %in% c(".", "AM", "ML", "TD", "99", "98", "97", "00", ".", "00", "KW")) %>% 
  filter(hspcr_2 != "££1") %>%
  filter(hspcr_2 != "inhib") %>% 
  filter(!(nature == "PU" & sexe == "2"))%>% 
  filter(!(nature %in% c("PV", "VULV") & sexe == "1"))

nrow(file)
table(file$nature,file$year)
table(file$insee_dep,file$nature)
table(file$sexe,file$nature)
max(file$dateprelev)
# save --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS(file, "pcr_hsv\\clean_data\\input_models\\raw_data.RDS")

