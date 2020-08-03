#
# Diversities - Bird communities and Landscapes

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

# ---------------------------------------------- libraries -------------------------------------------------
library(tidyverse)
library(dplyr)
library(rdiversity)
library(sp)
library(sf)

setwd("~/Github/Thesis/New data")

# -----------------------------------------------  Data  --------------------------------------------------- 
Community_Data_2015_2017_max_tbl <- readRDS("Community_Data_2015_2017_max_tbl.rds")
#new data - each route is ID with a unique code StateNum_RouteNum

Community_Data_2015_2017_max_matrix <- readRDS("Community_Data_2015_2017_max_matrix.rds")
#This matrix contains rows as species and columns as routes (unique ID)
dim(Community_Data_2015_2017_max_matrix)


#REVIEW THIS
routes.birds <- as.factor(colnames(Community_Data_2015_2017_max_matrix))
routes.land <- as.factor(landscape_data$U_S_R_I)


dif <- setdiff(routes.land, routes.birds)

# ------------------------------------------ Diversity Calculations ---------------------------------------

#Create output dataframe for bird communities diversities
bird_diversity <- data.frame("Route"=colnames(Community_Data_2015_2017_max_matrix)) %>%
  separate(Route, c("State", "RouteNum"), "_") %>%  #creates State Num and Route Num columns 
  mutate("Unique_State_Route_ID" = paste0(State, "_", RouteNum)) #creates Unique_State_Route_ID



#creates metacommunity object (required from rdiversity)
meta.data <- metacommunity(Community_Data_2015_2017_max_matrix)

# 1 - Redundancy (q=0 and q=1) ----
redq0 <- raw_sub_rho(meta.data, 0)
redq1 <- raw_sub_rho(meta.data, 1)

#cleaning rdiversity outputs to match the output dataframe
redq0 <- redq0 %>% select(partition_name, diversity) %>% 
  dplyr::rename("Unique_State_Route_ID"=partition_name, "Redundancy_q0"=diversity)

redq1 <- redq1 %>%  select(partition_name, diversity) %>% 
  dplyr::rename("Unique_State_Route_ID" = partition_name, "Redundancy_q1"=diversity)

#merging information with output dataframe
bird_diversity <- bird_diversity %>% merge(redq0, by="Unique_State_Route_ID") %>% 
  merge(redq1, by="Unique_State_Route_ID")


# 2 - Representativeness (q=0 and q=1) ---- 
repq0 <- norm_sub_rho(meta.data, 0)
repq1 <- norm_sub_rho(meta.data, 1)

#cleaning rdiversity outputs to match the output dataframe
repq0 <- repq0 %>% select(partition_name, diversity) %>% 
  dplyr::rename("Unique_State_Route_ID"=partition_name, "Representativeness_q0"=diversity)

repq1 <- repq1 %>%  select(partition_name, diversity) %>% 
  dplyr::rename("Unique_State_Route_ID" = partition_name, "Representativeness_q1"=diversity)

#merging information with output dataframe
bird_diversity <- bird_diversity %>% merge(repq0, by="Unique_State_Route_ID") %>% 
  merge(repq1, by="Unique_State_Route_ID")

# 3 - Alpha Diversity (q=0, q=1) ----
alphaq0 <- norm_sub_alpha(meta.data, 0)
alphaq1 <- norm_sub_alpha(meta.data, 1)

#cleaning rdiversity outputs to match the output dataframe
alphaq0 <- alphaq0 %>% select(partition_name, diversity) %>% 
  dplyr::rename("Unique_State_Route_ID"=partition_name, "Alpha_q0"=diversity)

alphaq1 <- alphaq1 %>%  select(partition_name, diversity) %>% 
  dplyr::rename("Unique_State_Route_ID" = partition_name, "Alpha_q1"=diversity)


#merging information with output dataframe
bird_diversity <- bird_diversity %>% merge(alphaq0, by="Unique_State_Route_ID") %>% 
  merge(alphaq1, by="Unique_State_Route_ID")

bird_diversity <- bird_diversity %>% select(-State, -RouteNum)

#check for missing values
sum(is.na(bird_diversity))


# ------------------------------------------- Landscape data -------------------------------------------------
#
# Landscape codes: 
#11 - Open Water (WATER)
#12 - Ice/snow (ICE)
#21 - Developed, Open Space (URBAN_OPEN)
#22 - Developed, Low Intensity (URBAN_LOW)
#23 - Developed, Medium Intensity (URBAN_MEDIUM)
#24 - Developed, High Intensity (URBAN_HIGH) #survey was conducted in secondary roads
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

#read data file
landscape_data <- read.csv("landscape_data_buffers_2016.csv")
dim(landscape_data)

#rename the columns

landscape.16 <- landscape_data %>% 
  rename("WATER"=3, "ICE"=4, "URBAN_OPEN"=5, "URBAN_LOW"=6,
         "URBAN_MEDIUM"=7, "URBAN_HIGH"=8, "BARREN"=9, "DECIDUOUS_FOREST"=10,
         "EVERGREEN_FOREST"=11,"MIXED_FOREST"=12, "SHRUBLAND"=13,
         "GRASSLAND"=14,"PASTURE"=15, "CROPLAND"=16,"WO0DY_WETLAND"=17,
         "HERB_WETLAND"=18) %>% 
  mutate(total = rowSums(landscape_data[, c(3:18)])) 
  

#create metacommunity obeject
#metacommunity object input needs to have routes as columns and landscape as rows

rownames(landscape.16)<- landscape.16$U_S_R_I 
landscape.16 <- landscape.16 %>%  select(-U_S_R_I)

meta.landscape <- metacommunity(t(landscape.16))




#---------------------------------------- Landscape Diversities --------------------------------------------------

#Create output dataframe for bird communities diversities
land_diversity <- data.frame("Unique_State_Route_ID"=rownames(landscape.16)) 



# 1 - Redundancy (q=0 and q=1) ----
land_redq0 <- raw_sub_rho(meta.landscape, 0)
land_redq1 <- raw_sub_rho(meta.landscape, 1)

#cleaning rdiversity outputs to match the output dataframe
land_redq0 <- land_redq0 %>% select(partition_name, diversity) %>% 
  dplyr::rename("Unique_State_Route_ID"=partition_name, "Landscape_Redundancy_q0"=diversity)

land_redq1 <- land_redq1 %>%  select(partition_name, diversity) %>% 
  dplyr::rename("Unique_State_Route_ID"=partition_name, "Landscape_Redundancy_q1"=diversity)

#merging information with output dataframe
land_diversity <- land_diversity %>% merge(land_redq0, by="Unique_State_Route_ID") %>% 
  merge(land_redq1, by="Unique_State_Route_ID")


# 2 - Representativeness (q=0 and q=1) ---- 
land_repq0 <- norm_sub_rho(meta.landscape, 0)
land_repq1 <- norm_sub_rho(meta.landscape, 1)

#cleaning rdiversity outputs to match the output dataframe
land_repq0 <- land_repq0 %>% select(partition_name, diversity) %>% 
  dplyr::rename("Unique_State_Route_ID"=partition_name, "Landscape_Representativeness_q0"=diversity)

land_repq1 <- land_repq1 %>%  select(partition_name, diversity) %>% 
  dplyr::rename("Unique_State_Route_ID" = partition_name, "Landscape_Representativeness_q1"=diversity)

#merging information with output dataframe
land_diversity <- land_diversity %>% merge(land_repq0, by="Unique_State_Route_ID") %>% 
  merge(land_repq1, by="Unique_State_Route_ID")

# 3 - Alpha Diversity (q=0, q=1) ----
land_alphaq0 <- norm_sub_alpha(meta.landscape, 0)
land_alphaq1 <- norm_sub_alpha(meta.landscape, 1)

#cleaning rdiversity outputs to match the output dataframe
land_alphaq0 <- land_alphaq0 %>% select(partition_name, diversity) %>% 
  dplyr::rename("Unique_State_Route_ID"=partition_name, "Landscape_Alpha_q0"=diversity)

land_alphaq1 <- land_alphaq1 %>%  select(partition_name, diversity) %>% 
  dplyr::rename("Unique_State_Route_ID" = partition_name, "Landscape_Alpha_q1"=diversity)


#merging information with output dataframe
land_diversity <- land_diversity %>% merge(land_alphaq0, by="Unique_State_Route_ID") %>% 
  merge(land_alphaq1, by="Unique_State_Route_ID")

sum(is.na(output_data))



# ------------------------------------------- Clean output data -------------------------------------------

bird.routes <- bird_diversity$Unique_State_Route_ID
land.routes <- land_diversity$Unique_State_Route_ID
different.land_bird <- as.factor(setdiff(land.routes, bird.routes)) #land has 725 different routes from bird data
different.bird_land <- as.factor(setdiff(bird.routes, land.routes)) #bird has 524 different routes from land data

test<- left_join(x=bird_diversity, land_diversity, by=c("Unique_State_Route_ID"="Unique_State_Route_ID"))
# takes all the routes from x (bird) and matches with land routes. 
# I believe this is done wrong because not all bird routes exist in the land cover and somehow the ouput
#has the same rows as the bird data 


test.1<- inner_join(x=bird_diversity, land_diversity, by=c("Unique_State_Route_ID"="Unique_State_Route_ID"))
# inner_join only matches rows present in both dataframes
# 2054 routes 

write.csv(test.1, "~/Github/Thesis/New data/Diversities.csv")


sum(test.1$Alpha_q0==62.2)




# P L O T S ---- 

#Redundancy
ggplot(test.1, aes(Landscape_Redundancy_q0, Redundancy_q0))+
  geom_point(alpha=0.2) + 
  geom_smooth(method="glm")+
  geom_rug() +
  ggtitle("Communities Redundancy against Landscape Redundancy", "q=0")


ggplot(test.1, aes(Landscape_Redundancy_q1, Redundancy_q1))+
  geom_point(alpha=0.2) + 
  geom_smooth(method="glm")+
  geom_rug()+
  ggtitle("Communities Redundancy against Landscape Redundancy", "q=1")


#Representativeness
ggplot(test.1, aes(Landscape_Representativeness_q0, Representativeness_q0))+
  geom_point(alpha=0.2) + 
  geom_smooth(method="glm")+
  geom_rug() +
  ggtitle("Communities Representativeness against Landscape Representativeness", "q=0")



ggplot(test.1, aes(Landscape_Representativeness_q1, Representativeness_q1))+
  geom_point(alpha=0.2) + 
  geom_smooth(method="glm")+
  geom_rug()+
  ggtitle("Communities Representativeness against Landscape Representativeness", "q=1")


#Alpha
ggplot(test.1, aes(Landscape_Alpha_q0, Alpha_q0))+
  geom_point(alpha=0.2) + 
  geom_smooth(method="glm")+
  geom_rug() +
  ggtitle("Communities Species Richeness against Landscape Richeness", "Alpha q=0")



ggplot(test.1, aes(Landscape_Alpha_q1, Alpha_q1))+
  geom_point(alpha=0.2) + 
  geom_smooth(method="glm")+
  geom_rug()+
  ggtitle("Communities Simpson Diversity against Landscape Simpson Diversity", "Alpha q=1")



####

ggplot(test.1, aes(Landscape_Alpha_q0, Landscape_Redundancy_q1))+
  geom_point(alpha=0.2) + 
  geom_smooth(method="glm")+
  geom_rug() +
  ggtitle("Landscape Richness against Landscape Redundancy","q=1")











