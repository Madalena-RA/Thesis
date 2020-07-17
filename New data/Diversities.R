#
# Diversities - Bird communities and landscapes

#The output: dataframe with :
# 1 - State
# 2 - Route Number
# 3 - Route unique ID
# 4 - Redundancy (q=0, q=1)
# 6 - Representativeness (q=0, q=1)
# 8 - Alpha diversity (q=0, q=1)
# 10 - Landscape redundancy (q=0, q=1)
# 12 - Landscape representativeness (q=0, q=1)
# 14 - Landscape alpha (q=0, q=1)
# 15 - BCR

# ----------------------------------------------libraries---------------------------------------------------
library(tidyverse)
library(dplyr)
library(rdiversity)

setwd("~/Github/Thesis/New data")

# -----------------------------------------------  Data  --------------------------------------------------- 
Community_Data_2015_2017_max_tbl <- readRDS("Community_Data_2015_2017_max_tbl.rds")
#new data - each route is ID with a unique code StateNum_RouteNum

Community_Data_2015_2017_max_matrix <- readRDS("Community_Data_2015_2017_max_matrix.rds")
#This matrix contains rows as species and columns as routes (unique ID)


#Create output dataframe
output_data <- data.frame("Route"=colnames(Community_Data_2015_2017_max_matrix)) %>%
  separate(Route, c("State", "RouteNum"), "_") %>%  #creates State Num and Route Num columns 
  mutate("Unique_State_Route_ID" = paste0(State, "_", RouteNum)) #creates Unique_State_Route_ID


# ------------------------------------------ Diversity Calculations ---------------------------------------

#creates metacommunity object (required from rdiversity)
meta.data <- metacommunity(Community_Data_2015_2017_max_matrix)

# 1 - Redundancy (q=0 and q=1) ----
redq0 <- raw_sub_rho(meta.data, 0)
redq1 <- raw_sub_rho(meta.data, 1)

#cleaning rdiversity outputs to match the output dataframe
redq0 <- redq0 %>% select(partition_name, diversity) %>% 
  dplyr::rename("Unique_State_Route_ID"=partition_name, "Redundancy q0"=diversity)

redq1 <- redq1 %>%  select(partition_name, diversity) %>% 
  dplyr::rename("Unique_State_Route_ID" = partition_name, "Redundancy q1"=diversity)

#merging information with output dataframe
output_data <- output_data %>% merge(redq0, by="Unique_State_Route_ID") %>% 
  merge(redq1, by="Unique_State_Route_ID")


# 2 - Representativeness (q=0 and q=1) ---- 
repq0 <- norm_sub_rho(meta.data, 0)
repq1 <- norm_sub_rho(meta.data, 1)

#cleaning rdiversity outputs to match the output dataframe
repq0 <- repq0 %>% select(partition_name, diversity) %>% 
  dplyr::rename("Unique_State_Route_ID"=partition_name, "Representativeness q0"=diversity)

repq1 <- repq1 %>%  select(partition_name, diversity) %>% 
  dplyr::rename("Unique_State_Route_ID" = partition_name, "Representativeness q1"=diversity)

#merging information with output dataframe
output_data <- output_data %>% merge(repq0, by="Unique_State_Route_ID") %>% 
  merge(repq1, by="Unique_State_Route_ID")

# 3 - Alpha Diversity (q=0, q=1) ----
alphaq0 <- norm_sub_alpha(meta.data, 0)
alphaq1 <- norm_sub_alpha(meta.data, 1)

#cleaning rdiversity outputs to match the output dataframe
alphaq0 <- alphaq0 %>% select(partition_name, diversity) %>% 
  dplyr::rename("Unique_State_Route_ID"=partition_name, "Alpha q0"=diversity)

alphaq1 <- alphaq1 %>%  select(partition_name, diversity) %>% 
  dplyr::rename("Unique_State_Route_ID" = partition_name, "Alpha q1"=diversity)


#merging information with output dataframe
output_data <- output_data %>% merge(alphaq0, by="Unique_State_Route_ID") %>% 
  merge(alphaq1, by="Unique_State_Route_ID")

# ------------------------------------------- Landscape data ---------------------------------------------

#11 - Open Water (WATER)
#12 - Ice/snow (ICE)
#21 - Developed, Open Space (URBAN_OPEN)
#22 - Developed, Low Intensity (URBAN_LOW)
#23 - Developed, Medium Intensity (URBAN_MEDIUM)
#24 - Developed, High Intensity (URBAN_HIGH)
#31 - Barren Land (BARREN)
#41 - Deciduous Forest (DECIDUOUS_FOREST)
#42 - Evergreen Forest (EVERGREEN_FOREST)
#43 - Mixed Forest (MIXED_FOREST)
#52 - Shrub/scrub (SHRUBLAND)
#71 - Grassland/Herbaceous (GRASSLAND)
#81 - Pasture (PASTURE)
#82 - Cultivated Crops (CROPLAND)
#90 - Woody Wetland (WOODY_WETLAND)
#95 - Emergent Herbaceous Wetland (HERB_WETLAND)

#Data: Unique_Route_ID + landscape codes

#rename the columns

#landscape_data <- #dataframe %>% 
  mutate("Unique_State_Route_ID" = paste0(State, "_", RouteNum)) #confirm name of state and route num columns
  #Renaming the columns to the corret landscape
  rename("WATER"=#, "ICE"=#, "URBAN_OPEN"=#, "URBAN_LOW"=#,
         "URBAN_MEDIUM"=#, "URBAN_HIGH"=#, "BARREN"=#, "DECIDUOUS_FOREST"=9,
         "EVERGREEN_FOREST"=#,"MIXED_FOREST"=#, "SHRUBLAND"=#,
         "GRASSLAND"=#,"PASTURE"=#, "CROPLAND"=#,"WODDY_WETLAND"=#,
         "HERB_WETLAND"=#, %>% 
  #
  mutate(total = rowSums(select(.,-RouteName))) %>% 
  mutate("Urban"=DevOpen+DevLow+DevMed+DevHigh) %>% 
  mutate("Forest"=DeciduousForest+EvergreenForest+MixedForest) %>% 
  merge(bcr, by="RouteName")








