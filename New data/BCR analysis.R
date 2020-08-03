# 
#
# Calculates the Averages of the indices while grouping the routes by BCR
# 
# Creates two new shapefiles: 
# - Averages_by_routes
#   each route has the averages of their BCR
# - Averages_by_BCR 
#   

# ---------------------------------------Set up------------------------------------------------------------

#
# 
library(dplyr)
library(rgdal) #spatial objects manipulation


Diversities <- read.csv("~/Github/Thesis/New data/Diversities.csv")
routes <- read.csv("~/GLASGOW/Dissertation/data/routes.csv")

#Adding the Unique ID for each route

#Making sure that state are two-digit numbers and Routes are three-digit numbers
routes$StateNum <- str_pad(routes$StateNum, width = 2,pad = "0")
routes$Route <- str_pad(routes$Route, width = 3,pad = "0")

routes <- routes %>% 
  mutate("U_S_R_I"=paste0(StateNum, "_", Route))

#only BCR and U_S_R_I
BCR <- routes %>% select(U_S_R_I, BCR)

Diversities_BCR <- left_join(Diversities, BCR, by=c("Unique_State_Route_ID"="U_S_R_I"))

write.csv(Diversities_BCR, "~/Github/Thesis/New data/Diversities_BCR.csv")



#---------------------------------------- Averages by BCR -------------------------------------------------------------------

#Calculates the Average of every column by BCR 
#The final dataframe has a row for each of the 30 BCR and their respective average value 
#for each variable
Averages <- Diversities_BCR %>% select(-X, -Unique_State_Route_ID) %>% 
  group_by(BCR) %>% summarise_all(list(mean)) 

#Prefixing the columns with "Mean." to avoid mistakes with the indices previously calculated
colnames(Averages) <- paste0("Mean.",colnames(Averages))

#Fixing the BCR column (there is no mean for BCR)
Averages <- Averages %>% rename("BCR"=Mean.BCR)


#Adding the new information for every route 
#This will only be use to create the Averages_by_routes shapefile
Averages.routes <- left_join(Diversities_BCR, Averages, by=c("BCR"="BCR"), )




#----------------------------------- Spatial Data ----------------------------------------------

#reads the shapefile
Routes <- readOGR("~/Github/Thesis/New data/Routes_Compiled.shp", GDAL1_integer64_policy = T)

#Extracts the slot of the shapefile that contains the information as a dataframe
# which is easier to manipulate
routes.data <- Routes@data 

#Changes the information present in the dataframe.  Adds information ("Averages")
routes.data<- left_join(routes.data, Averages.routes, by=c("U_S_R_I"="Unique_State_Route_ID"))

#redefines the dataframe as the slot of the shapefile
Routes@data<- routes.data

#Saves the new shapefile (Averages_by_Routes)
writeOGR(obj= Routes,dsn ="~/Github/Thesis/New data",layer = "Averages_by_routes", driver="ESRI Shapefile", check_exists = FALSE )


#reads BCR shapefile
BCR.shapefile <- readOGR("~/Github/Thesis/New data/USA_only_BCR.shp", GDAL1_integer64_policy = T)

#Addinf the averages to the shapefile - by BCR
BCR.data <- BCR.shapefile@data
BCR.data <- BCR.data %>% select(-AVRed, -AVRep)

BCR.data <- left_join(BCR.data, Averages, by=c("BCR"="BCR"))

BCR.shapefile@data <-BCR.data

#saves the new shapefile (Averages_by_BCR)
writeOGR(obj= BCR.shapefile,dsn ="~/Github/Thesis/New data",layer = "Averages_by_BCR", driver="ESRI Shapefile", check_exists = FALSE )

