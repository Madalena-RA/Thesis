
library(tidyverse)
library(dplyr)
library(betapart)
library(rdiversity)

#working directory
setwd("~/GLASGOW/Dissertation/data")





hist_01_23km_LC <- read.csv("~/GLASGOW/Dissertation/data/hist_01_23km_LC.csv")

hist_01_23km_LC <- hist_01_23km_LC %>% 
  rename("RouteName"=1, "OpenWater"=3, "Snow"=4, "DevOpen"=5,
         "DevLow"=6, "DevMed"=7, "DevHigh"=8, "Barren"=9,
         "DeciduousForest"=10,"EvergreenForest"=11, "MixedForest"=12,
         "Shrub"=13,"Grassland"=14, "Pasture"=15,"Crops"=16,
         "WoodyWetlands"=17,"HerbaceousWetland"=18) %>% 
  mutate(total = rowSums(select(.,-RouteName))) %>% 
  mutate("Urban"=DevOpen+DevLow+DevMed+DevHigh) %>% 
  mutate("Forest"=DeciduousForest+EvergreenForest+MixedForest) %>% 
  merge(bcr, by="RouteName")

habitat.2001 <-hist_01_23km_LC %>% select(-total, -Forest, -Urban, -X23HIST_0)

hist_01_23km_LC[,c(2:19)] <- (hist_01_23km_LC[,c(2:19)]/hist_01_23km_LC[,c(2:19)]$total)*100

Hist_01_23km_LC <- index.df %>% select(-Red_bcr.2016,-Rep_bcr.2016,-Red_usa.2016,-Rep_usa.2016) %>% 
  merge(hist_01_23km_LC, by="RouteName")



# Redundancy of habitat 2001 ----
habitat.df.01 <- data.frame()

for (r in bcr.usa) {
  data <- habitat.2001 %>% filter(BCR==r) %>% column_to_rownames("RouteName") %>% select(-BCR)
  meta.data <- metacommunity(t(data))
  redundancy <- raw_sub_rho(meta.data,1)
  redundancy.clean <- redundancy %>% rename("RouteName"=partition_name) %>% select(RouteName, diversity) %>% 
    merge(bcr,by="RouteName")
  habitat.df.01 <- rbind(habitat.df.01, redundancy.clean)
}

habitat.df.01 <- habitat.df.01 %>%  rename("Red_habitat"=diversity) %>% select(-BCR) %>% 
  merge(Hist_01_23km_LC, by="RouteName")

ggplot(habitat.df.01, aes(Red_habitat,Red_bcr.2001)) +
  geom_point() +
  geom_smooth(method = "loess") +
  geom_rug() +
  ggtitle("Redundancy of bird community as a function of the redundancy of the landscape", "within regions, 2001")


# Representativeness of habitat 2001 ----
habitat.rep.df.01 <- data.frame()

for (r in bcr.usa) {
  data <- habitat.2001 %>% filter(BCR==r) %>% column_to_rownames("RouteName") %>% select(-BCR)
  meta.data <- metacommunity(t(data))
  Representativeness <- norm_sub_rho(meta.data,1)
  representativeness.clean <- Representativeness %>% rename("RouteName"=partition_name) %>% select(RouteName, diversity) %>% 
    merge(bcr,by="RouteName")
  habitat.rep.df.01 <- rbind(habitat.rep.df.01, representativeness.clean)
}

habitat.rep.df.01 <- habitat.rep.df.01 %>%  rename("Rep_habitat"=diversity) %>% select(-BCR) %>% 
  merge(Hist_01_23km_LC, by="RouteName")

ggplot(habitat.rep.df.01, aes(Rep_habitat,Rep_bcr.2001)) +
  geom_point() +
  geom_smooth(method = "loess") +
  geom_rug() +
  ggtitle("Representativeness of bird community as a function of the Representativeness of the landscape", "within regions, 2001")





# PLOTS ---- 
ggplot(Hist_01_23km_LC, aes(x=Forest, Red_bcr.2001)) +
  geom_point(alpha=0.2, col="darkgreen") +
  ylab("Redundancy")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Redundancy q=1 as a function of proportion of Forest, 2001","Redundancy within regions")

ggplot(Hist_01_23km_LC, aes(x=DeciduousForest, Red_bcr.2001)) +
  geom_point(alpha=0.2, col="darkgreen") +
  ylab("Redundancy")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Redundancy q=1 as a function of proportion of Deciduous Forest, 2001","Redundancy within regions")

ggplot(Hist_01_23km_LC, aes(x=EvergreenForest, Red_bcr.2001)) +
  geom_point(alpha=0.2, col="darkgreen") +
  ylab("Redundancy")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Redundancy q=1 as a function of proportion of Evergreen Forest, 2001","Redundancy within regions")

ggplot(Hist_01_23km_LC, aes(x=MixedForest, Red_bcr.2001)) +
  geom_point(alpha=0.2, col="darkgreen") +
  ylab("Redundancy")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Redundancy q=1 as a function of proportion of Mixed Forest, 2001","Redundancy within regions")

#urban
ggplot(Hist_01_23km_LC, aes(x=Urban, Red_bcr.2001)) +
  geom_point(alpha=0.2, col="black") +
  ylab("Redundancy")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Redundancy q=1 as a function of proportion of Urban , 2001","Redundancy within regions")


ggplot(Hist_01_23km_LC, aes(x=DevOpen, Red_bcr.2001)) +
  geom_point(alpha=0.2, col="black") +
  ylab("Redundancy")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Redundancy q=1 as a function of proportion of Urban Open Space, 2001","Redundancy within regions")

ggplot(Hist_01_23km_LC, aes(x=DevLow, Red_bcr.2001)) +
  geom_point(alpha=0.2, col="black") +
  ylab("Redundancy")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Redundancy q=1 as a function of proportion of Urban Low, 2001","Redundancy within regions")

ggplot(Hist_01_23km_LC, aes(x=DevMed, Red_bcr.2001)) +
  geom_point(alpha=0.2, col="black") +
  ylab("Redundancy")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Redundancy q=1 as a function of proportion of Urban Medium, 2001","Redundancy within regions")

ggplot(Hist_01_23km_LC, aes(x=DevHigh, Red_bcr.2001)) +
  geom_point(alpha=0.2, col="black") +
  ylab("Redundancy")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Redundancy q=1 as a function of proportion of Urban High, 2001","Redundancy within regions")

#agriculture
ggplot(Hist_01_23km_LC, aes(x=Pasture, Red_bcr.2001)) +
  geom_point(alpha=0.2, col="coral") +
  ylab("Redundancy")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Redundancy q=1 as a function of proportion of Pasture, 2001","Redundancy within regions")

ggplot(Hist_01_23km_LC, aes(x=Crops, Red_bcr.2001)) +
  geom_point(alpha=0.2, col="coral") +
  ylab("Redundancy")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Redundancy q=1 as a function of proportion of Cropland, 2001","Redundancy within regions")


#grassland and shrubland

ggplot(Hist_01_23km_LC, aes(x=Shrub, Red_bcr.2001)) +
  geom_point(alpha=0.2, col="brown") +
  ylab("Redundancy")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Redundancy q=1 as a function of proportion of Shrubland, 2001","Redundancy within regions")

ggplot(Hist_01_23km_LC, aes(x=Grassland, Red_bcr.2001)) +
  geom_point(alpha=0.2, col="brown") +
  ylab("Redundancy")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Redundancy q=1 as a function of proportion of Grassland, 2001","Redundancy within regions")


#Wetlands
ggplot(Hist_01_23km_LC, aes(x=WoodyWetlands, Red_bcr.2001)) +
  geom_point(alpha=0.2, col="mediumpurple4") +
  ylab("Redundancy")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Redundancy q=1 as a function of proportion of Woody Wetland, 2001","Redundancy within regions")

ggplot(Hist_01_23km_LC, aes(x=HerbaceousWetland, Red_bcr.2001)) +
  geom_point(alpha=0.2, col="mediumpurple4") +
  ylab("Redundancy")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Redundancy q=1 as a function of proportion of Herbaceous Wetland, 2001","Redundancy within regions")



### REP ----

ggplot(Hist_01_23km_LC, aes(x=Forest, Rep_bcr.2001)) +
  geom_point(alpha=0.2, col="darkgreen") +
  ylab("Representativeness")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Representativeness q=1 as a function of proportion of Forest, 2001","Representativeness within regions")

ggplot(Hist_01_23km_LC, aes(x=DeciduousForest, Rep_bcr.2001)) +
  geom_point(alpha=0.2, col="darkgreen") +
  ylab("Representativeness")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Representativeness q=1 as a function of proportion of Deciduous Forest, 2001","Representativeness within regions")

ggplot(Hist_01_23km_LC, aes(x=EvergreenForest, Rep_bcr.2001)) +
  geom_point(alpha=0.2, col="darkgreen") +
  ylab("Representativeness")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Representativeness q=1 as a function of proportion of Evergreen Forest, 2001","Representativeness within regions")

ggplot(Hist_01_23km_LC, aes(x=MixedForest, Rep_bcr.2001)) +
  geom_point(alpha=0.2, col="darkgreen") +
  ylab("Representativeness")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Representativeness q=1 as a function of proportion of Mixed Forest, 2001","Representativeness within regions")

#urban
ggplot(Hist_01_23km_LC, aes(x=Urban, Rep_bcr.2001)) +
  geom_point(alpha=0.2, col="black") +
  ylab("Representativeness")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Representativeness q=1 as a function of proportion of Urban , 2001","Representativeness within regions")


ggplot(Hist_01_23km_LC, aes(x=DevOpen, Rep_bcr.2001)) +
  geom_point(alpha=0.2, col="black") +
  ylab("Representativeness")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Representativeness q=1 as a function of proportion of Urban Open Space, 2001","Representativeness within regions")

ggplot(Hist_01_23km_LC, aes(x=DevLow, Rep_bcr.2001)) +
  geom_point(alpha=0.2, col="black") +
  ylab("Representativeness")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Representativeness q=1 as a function of proportion of Urban Low, 2001","Representativeness within regions")

ggplot(Hist_01_23km_LC, aes(x=DevMed, Rep_bcr.2001)) +
  geom_point(alpha=0.2, col="black") +
  ylab("Representativeness")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Representativeness q=1 as a function of proportion of Urban Medium, 2001","Representativeness within regions")

ggplot(Hist_01_23km_LC, aes(x=DevHigh, Rep_bcr.2001)) +
  geom_point(alpha=0.2, col="black") +
  ylab("Representativeness")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Representativeness q=1 as a function of proportion of Urban High, 2001","Representativeness within regions")

#agriculture
ggplot(Hist_01_23km_LC, aes(x=Pasture, Rep_bcr.2001)) +
  geom_point(alpha=0.2, col="coral") +
  ylab("Representativeness")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Representativeness q=1 as a function of proportion of Pasture, 2001","Representativeness within regions")

ggplot(Hist_01_23km_LC, aes(x=Crops, Rep_bcr.2001)) +
  geom_point(alpha=0.2, col="coral") +
  ylab("Representativeness")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Representativeness q=1 as a function of proportion of Cropland, 2001","Representativeness within regions")


#grassland and shrubland

ggplot(Hist_01_23km_LC, aes(x=Shrub, Rep_bcr.2001)) +
  geom_point(alpha=0.2, col="brown") +
  ylab("Representativeness")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Representativeness q=1 as a function of proportion of Shrubland, 2001","Representativeness within regions")

ggplot(Hist_01_23km_LC, aes(x=Grassland, Rep_bcr.2001)) +
  geom_point(alpha=0.2, col="brown") +
  ylab("Representativeness")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Representativeness q=1 as a function of proportion of Grassland, 2001","Representativeness within regions")


#Wetlands
ggplot(Hist_01_23km_LC, aes(x=WoodyWetlands, Rep_bcr.2001)) +
  geom_point(alpha=0.2, col="mediumpurple4") +
  ylab("Representativeness")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Representativeness q=1 as a function of proportion of Woody Wetland, 2001","Representativeness within regions")

ggplot(Hist_01_23km_LC, aes(x=HerbaceousWetland, Rep_bcr.2001)) +
  geom_point(alpha=0.2, col="mediumpurple4") +
  ylab("Representativeness")+
  geom_rug()+
  geom_smooth(method="loess", col="black") + 
  ggtitle("Representativeness q=1 as a function of proportion of Herbaceous Wetland, 2001","Representativeness within regions")







c <- c()
for (b in bcr.usa) {
  x<- sum(Hist_01_23km_LC$BCR==b)
  c<-c(c,x)
}
table<-cbind(bcr.usa,c)
print(order(table, bcr.usa))


