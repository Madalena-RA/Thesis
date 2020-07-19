# 
# Modeling
#
# 
library(dplyr) #Data Manipulation
library(lme4) #for the model

# 
# ---------------------------------------- Data ----------------------------------------------------

Diversities <- read.csv("~/Github/Thesis/New data/Diversities.csv")
# Diversities if a dataframe with information on the 
#Redudancy (q=0,q=1), Representativeness (q=0,q=1) and Alpha (q=0,q=1) for 
#bird communities data and landscape data for 2504 routes with unique ID
# Diversities dataframe contains only the routes that had information available
# for both landcover and bird community data 



landscape_data_buffers_2016 <- read.csv("~/Github/Thesis/New data/landscape_data_buffers_2016.csv")
# landscape_data_buffers_2016 is a dataframe created on QGIS with the values of each type
#of landcover in the 23km buffers centred in the centroid of each route present in the 
#Routes_Compiled shapefile. 3229 total routes with a unique ID that matches the ID from Diversities

#Landscape data has more routes than Diversity 

#We need to select the 2504 diversity routes from the landscape data
#
Landscape <-  landscape_data_buffers_2016 %>% 
  rename("WATER"=3, "ICE"=4, "URBAN_OPEN"=5, "URBAN_LOW"=6,
         "URBAN_MEDIUM"=7, "URBAN_HIGH"=8, "BARREN"=9, "DECIDUOUS_FOREST"=10,
         "EVERGREEN_FOREST"=11,"MIXED_FOREST"=12, "SHRUBLAND"=13,
         "GRASSLAND"=14,"PASTURE"=15, "CROPLAND"=16,"WOODY_WETLAND"=17,
         "HERB_WETLAND"=18) %>% 
  mutate(total = rowSums(landscape_data_buffers_2016[, c(3:18)])) 

#Change the values of each landcover to proportions

rownames(Landscape)<- Landscape$U_S_R_I 
Landscape <- Landscape %>%  select(-U_S_R_I)

Landscape.proportions <- Landscape[,1:17] / Landscape[,18] 
Landscape.proportions$U_S_R_I <- rownames(Landscape.proportions)

#Merge the landscape information with the indices dataframe
data <- inner_join(Diversities, Landscape.proportions, by=c("Unique_State_Route_ID"="U_S_R_I"))
data <- data %>% rename("ID"=X)
#
# Histograms
hist(data$Redundancy_q0,50)
hist(data$Redundancy_q1,50)
hist(data$Representativeness_q0,50)
hist(data$Representativeness_q1,50)
hist(data$Alpha_q0,50)
hist(data$Alpha_q1,50)

hist(data$Landscape_Redundancy_q0, 50)
hist(data$Landscape_Redundancy_q1, 50)
hist(data$Landscape_Representativeness_q0, 50)
hist(data$Landscape_Representativeness_q1, 50)
hist(data$Landscape_Alpha_q0, 50)
hist(data$Landscape_Alpha_q1, 50)

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

# ----------------------------------------------- MODELS --------------------------------------------------
# REDUNDANCY MODELS ????????
Red.model1 <- glm(Redundancy_q0 ~ WATER + ICE + 
                    URBAN_OPEN + URBAN_LOW + URBAN_MEDIUM + URBAN_HIGH+
                    BARREN + DECIDUOUS_FOREST + EVERGREEN_FOREST + MIXED_FOREST +
                    SHRUBLAND + GRASSLAND + PASTURE + CROPLAND +
                    WOODY_WETLAND + HERB_WETLAND +
                    Landscape_Redundancy_q0+I(Landscape_Redundancy_q0^2)+
                    Landscape_Redundancy_q1 + I(Landscape_Redundancy_q1^2) +
                    Landscape_Alpha_q0 + I(Landscape_Alpha_q0^2) +
                    Landscape_Alpha_q1 + I(Landscape_Alpha_q1^2) +
                    Alpha_q0 + Alpha_q1 +
                    (1|ID), data=data)
summary(Red.model1)

Red.model2 <- glm(Redundancy_q1 ~ WATER + ICE + 
                    URBAN_OPEN + URBAN_LOW + URBAN_MEDIUM + URBAN_HIGH+
                    BARREN + DECIDUOUS_FOREST + EVERGREEN_FOREST + MIXED_FOREST +
                    SHRUBLAND + GRASSLAND + PASTURE + CROPLAND +
                    WOODY_WETLAND + HERB_WETLAND, data=data)
summary(Red.model2)
































