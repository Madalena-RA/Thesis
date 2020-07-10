#Spatial data - routes

#Creates a routes and BCR shapefile with a unique identifier for each route
#The unique Identifier is in the format RouteName_State_BCR 

#useful packages
library(tigris)
library(leaflet)
library(sp)
library(rgdal)



#Routes 3731 observations
routes.shp <- rgdal::readOGR("C:/Users/madal/Documents/GLASGOW/Dissertation/data/USBBS.routes.All.shp", GDAL1_integer64_policy =T)
routes.shp.data <- routes.shp@data

#BCR
BCR.shp <- readOGR("C:/Users/madal/Documents/GLASGOW/Dissertation/data/BCR_data/BCR_Terrestrial_master.shp", GDAL1_integer64_policy = T)
BCR.shp.data <- BCR.shp@data

#States
states.shp <- tigris::states(cb = FALSE, resolution = "500k", year = NULL)
states.shp.data <- states.shp@data


#Centroid
#gives the coordinates of the centroids for each route
#important to match each route to a BCR and a state
centroids <- gCentroid(routes.shp, byid=T) #Formal SpatialPoints with coordinates of the centroids 

#all have the same projection
routes.shp@proj4string
BCR.shp@proj4string
centroids@proj4string
states.shp@proj4string

#' *note: the over function() conserves the order of the rows therfore we can use cbind with no problem*

#Matching each route to a BCR
points_in_BCR <- sp::over(centroids, BCR.shp) #over() creates a dataframe with the same order of rows as the 1st argumet

#Matching each route to a State
points_in_state <- sp::over(centroids, states.shp)

routes.shp.data <- cbind(routes.shp.data, points_in_BCR$BCR)
routes.shp.data <- cbind(routes.shp.data, points_in_state$NAME)

#rename the columns
routes.shp.data <- routes.shp.data %>% 
                    dplyr::rename("BCR"=11, "state"=12) %>% 
                      mutate("state"=toupper(state))

#Each State has a number code
States_Num_Name <- read.csv("~/Github/Thesis/States_Num_Name.csv", sep=";")

States_Num_Name <- States_Num_Name %>% 
  filter(CountryNum=="840") %>% 
  mutate("state"=str_trim(State.Prov.TerrName)) %>% 
  select(-State.Prov.TerrName) %>% 
  mutate("state"=toupper(state))

#Match the code to the state
routes.shp.data <- merge(routes.shp.data, States_Num_Name, by="state", all.x=T, all.y=F)

#creates the unique identifier for each route
routes.shp.data$RteName_St_BCR<- paste0(routes.shp.data$RTENAME,
                                   "_", routes.shp.data$StateNum, "_", routes.shp.data$BCR)


#Check if there are duplicates
duplicates <- count(routes.shp.data, "RteName_St_BCR")
duplicated <- duplicates %>% subset(freq >1) 
duplicated <- duplicated$RteName_St_BCR #name of the routes present more than once in the routes data


multipled_routes <- routes.shp.data %>% dplyr::mutate(duplicated= if_else(RteName_St_BCR%in% duplicated, "yes","no"))

#only routes with unique names
unique<- subset(multipled_routes, duplicated=="no")
doubled <- subset(multipled_routes, duplicated=="yes")

#creates empty matrix to match number of routes in the shapefile
blanc <- matrix(nrow=769, ncol=16, NA) 
blanc <- data.frame(blanc) 
colnames(blanc)=colnames(unique)
 
#now the subseted dataframe has the same dimension as the original shapefile (no error messages)
unique <- rbind(unique, blanc)


#Now, the shape file only contains information for the routes that have unique names
routes.shp@data <- unique


writeOGR(obj= routes.shp,dsn ="C:/Users/madal/Documents/GLASGOW/Dissertation/data/Clean",layer = "unique_routes", driver="ESRI Shapefile", check_exists = FALSE )









