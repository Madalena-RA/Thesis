#Sampling effort 



#libraries
library(sp)
library(rgdal)
library(rgeos)
library(reshape2)
library(stringr)
library(ggplot2)

setwd("~/Github/Thesis/New data")


bcr.shp <- readOGR("BCR_Terrestrial_master.shp")
bcr.shp.data <-bcr.shp@data 

bcr.shp.data$Area <- gArea(spgeom = bcr.shp, byid = TRUE)

bcr.usa <- bcr.shp.data %>%
  subset(COUNTRY=="USA") %>% 
  subset(PROVINCE_S != "ALASKA") %>% 
  subset(PROVINCE_S != "HAWAIIAN ISLANDS")

bcr.usa.size <- bcr.usa %>% select(BCR, Area) 

bcr.usa.size <- aggregate(bcr.usa.size$Area, by=list(BCR=bcr.usa.size$BCR), FUN=sum)
bcr.usa.size <- bcr.usa.size %>% dplyr::rename("Area"=x)


#ADD BCR 
usbbs.routes <- read.csv("~/GLASGOW/Dissertation/data/routes.csv")

#standard state number and Route number
usbbs.routes$StateNum <- str_pad(usbbs.routes$StateNum, width = 2,pad = "0")
usbbs.routes$Route <- str_pad(usbbs.routes$Route, width = 3,pad = "0")


#creates the same unique code for each route
usbbs.routes$Unique_State_Route_ID <- paste0(usbbs.routes$StateNum,sep="_", usbbs.routes$Route)

#selects only BCR and route Code
bcr.code <- usbbs.routes %>%   select(BCR, Unique_State_Route_ID)

number.bcr <- count(x = bcr.code, BCR = BCR)





#BCR, Area, number of routes
bcr.usa.size <- left_join(bcr.usa.size, number.bcr, by=c("BCR"="BCR"), all.y=TRUE, all.x=TRUE)


bcr.usa.size <- bcr.size %>% mutate("proportion"=n/Area)
mean(bcr.size$proportion)

bcr.size <- bcr.size[-c(36,37), ] 

ggplot(bcr.size, aes(BCR, proportion))+
  geom_point()
  
  


ggplot(bcr.size, aes(BCR, Area, color=as.factor(BCR))) +
  geom_point(aes(size=n), alpha=0.5) +
  xlim(c(0,40))+
  ggtitle("BCR Area and number of routes")

prop <- na.omit(bcr.size$proportion)
