#'
#'*Analysis Beta diversities*
#'
#' 
#'

library(tidyverse)
library(dplyr)
library(betapart)
library(rdiversity)

#working directory
setwd("~/GLASGOW/Dissertation/data")

#abundance data

load("vegan.spmatrix_segments_t1.rda") #data from 2001
load("vegan.spmatrix_segments_t2.rda") #data from 2016
load("data.segments.4km.rda")
load("C:/Users/madal/Documents/GLASGOW/Dissertation/data/lc_01_16_23km.rda")
routes <- read.csv("C:/Users/madal/Downloads/routes/routes.csv")
states <- read.csv("~/GLASGOW/Dissertation/data/States_Num_Name.csv", sep=";")
bcr <- routes %>% select(RouteName, BCR)

transects <- mutate(data.segments.4km, RouteName= substr(data.segments.4km$partition,1,nchar(as.vector(data.segments.4km$partition))-2)) %>% 
  select(partition, RouteName)
transects.names <- unique(transects$RouteName)

# Clean data ----

columns.names.t1 <- colnames(spmatrix.t1)
columns.names.t2 <- colnames(spmatrix.t2)

missing.sp.2 <- setdiff(columns.names.t1, columns.names.t2)                            #species that matrix 1 has and 2 doesnt
missing.sp.1 <- setdiff(columns.names.t2, columns.names.t1)                            #species that matrix 2 has and 1 doesnt


#add new columns (with the missing species) set to zero

missing.species.t1 <- matrix(0,nrow = 4094,ncol = 24)                                  #creates empty matrix for the missing species
colnames(missing.species.t1) <- missing.sp.1                                           #adds missing species names
spmatrix.t1.full <- cbind(spmatrix.t1, missing.species.t1)                      #binds the two matrices
spmatrix.t1.full <- spmatrix.t1.full[, order(as.integer(colnames(spmatrix.t1.full)))]  #orders the columns

missing.species.t2 <- matrix(0, nrow = 4094, ncol=14)
colnames(missing.species.t2) <- missing.sp.2
spmatrix.t2.full <- cbind(spmatrix.t2, missing.species.t2)
spmatrix.t2.full <- spmatrix.t2.full[, order(as.integer(colnames(spmatrix.t2.full)))] 

spmatrix.t1_binary <- spmatrix.t1.full
spmatrix.t1_binary[spmatrix.t1_binary > 0] = 1      #binary data from 2001

spmatrix.t2_binary <- spmatrix.t2.full
spmatrix.t2_binary[spmatrix.t2_binary > 0] = 1      #binary data from 2016



#1.1  - NEW HABITAT DATA BY TRANSECTS (wrong) ----

#landscape proportions : sum
#delta temp: mean
#deleted all alpha diversity indices 

data.transects.4km <- data.segments.4km %>% 
  merge(transects, by="partition") %>% 
  select(-partition) 

temperature <- data.transects.4km %>% select(delta.temp, RouteName)
temperature <- aggregate(delta.temp~RouteName, data=temperature, FUN = mean)


data.transects.4km <- aggregate(.~RouteName, data=data.transects.4km, sum)
data.transects.4km <- data.transects.4km %>% select(-delta.temp)


data.transects.4km <- merge(data.transects.4km[,c(1:18)],temperature, by= "RouteName")


## 1.2 - 23 km buffers ---- 

lc.01.16.23km <- lc.01.16.23km %>%  rename("RouteName"=partition)

#transects 2001 landscape
lc.01.23km <- lc.01.16.23km %>% select(RouteName, water.t1,icesnow.t1,urban.t1,deciduousforest.t1,evergreenforest.t1,
                                       mixedforest.t1,grassland.t1,cropland.t1,shrubland.t1,wetland.t1, delta.urban)
#transects 2016 landscape
lc.16.23km <- lc.01.16.23km %>% select(RouteName, water.t2,icesnow.t2,urban.t2,deciduousforest.t2,evergreenforest.t2,
                                       mixedforest.t2,grassland.t2,cropland.t2,shrubland.t2,wetland.t2, delta.urban)



#'##################################################################################
# Distinctiveness approach (all indices) ----
#
#'##################################################################################

meta.2001 <- as.data.frame(spmatrix.t1.full) %>% 
  mutate("partition"=rownames(spmatrix.t1.full)) %>%  #creates new column with partitions
  merge(transects, by="partition") %>%                #adds transects names
  select(-partition)                                  #deletes partition column

meta.2001 <- aggregate(.~RouteName, data = meta.2001, sum) #adds observations per RouteName
rownames(meta.2001) <- meta.2001$RouteName                 #defines routenames as rownames
meta.2001 <- select(meta.2001, -RouteName)                 #deletes routenames as colums
meta.2001 <- metacommunity(t(meta.2001))                   #transforms object into a metacommunity object


meta.2016 <- as.data.frame(spmatrix.t2.full) %>% 
  mutate("partition"=rownames(spmatrix.t2.full)) %>% 
  merge(transects, by="partition") %>% 
  select(-partition) 

meta.2016 <- aggregate(.~RouteName, data = meta.2016, sum)
rownames(meta.2016) <- meta.2016$RouteName
meta.2016 <- select(meta.2016, -RouteName)
meta.2016 <- metacommunity(t(meta.2016))


# MEASURES 2001 ----

measures_2001 <- data.frame("RouteName"=transects.names,
                            "Red_q0"=0,"Red_q1"=0,"Red_q2"=0,"Red_qinf"=0,
                            "Dist_q0"=0,"Dist_q1"=0,"Dist_q2"=0,"Dist_qinf"=0,
                            "Rep_q0"=0,"Rep_q1"=0,"Rep_q2"=0,"Rep_qinf"=0,
                            "Eff_q0"=0,"Eff_q1"=0,"Eff_q2"=0,"Eff_qinf"=0)



# redundancy 2001 ----

red_q0.t1 <- raw_sub_rho(meta.2001, 0)
red_q1.t1 <- raw_sub_rho(meta.2001, 1)
red_q2.t1 <- raw_sub_rho(meta.2001, 2)
red_qinf.t1 <- raw_sub_rho(meta.2001, Inf)

measures_2001$Red_q0 <- red_q0.t1$diversity
measures_2001$Red_q1 <- red_q1.t1$diversity
measures_2001$Red_q2 <- red_q2.t1$diversity
measures_2001$Red_qinf <- red_qinf.t1$diversity

# distinctiveness 2001 ----

dist_q0.t1 <- raw_sub_beta(meta.2001, 0)
dist_q1.t1 <- raw_sub_beta(meta.2001, 1)
dist_q2.t1 <- raw_sub_beta(meta.2001, 2)
dist_qinf.t1 <- raw_sub_beta(meta.2001, Inf)

measures_2001$Dist_q0 <- dist_q0.t1$diversity
measures_2001$Dist_q1 <- dist_q1.t1$diversity
measures_2001$Dist_q2 <- dist_q2.t1$diversity
measures_2001$Dist_qinf <- dist_qinf.t1$diversity

# representativeness 2001 ----

rep_q0.t1 <- norm_sub_rho(meta.2001, 0)
rep_q1.t1 <- norm_sub_rho(meta.2001, 1)
rep_q2.t1 <- norm_sub_rho(meta.2001, 2)
rep_qinf.t1 <- norm_sub_rho(meta.2001, Inf)

measures_2001$Rep_q0 <- rep_q0.t1$diversity
measures_2001$Rep_q1 <- rep_q1.t1$diversity
measures_2001$Rep_q2 <- rep_q2.t1$diversity
measures_2001$Rep_qinf <- rep_qinf.t1$diversity

# effective number of distinct subcommunities 2001 ----

eff_q0.t1 <- norm_sub_beta(meta.2001, 0)
eff_q1.t1 <- norm_sub_beta(meta.2001, 1)
eff_q2.t1 <- norm_sub_beta(meta.2001, 2)
eff_qinf.t1 <- norm_sub_beta(meta.2001, Inf)

measures_2001$Eff_q0 <- eff_q0.t1$diversity
measures_2001$Eff_q1 <- eff_q1.t1$diversity
measures_2001$Eff_q2 <- eff_q2.t1$diversity
measures_2001$Eff_qinf <- eff_qinf.t1$diversity


# MEASURES 2016 ----

measures_2016 <- data.frame("RouteName"=transects.names,
                            "Red_q0"=0,"Red_q1"=0,"Red_q2"=0,"Red_qinf"=0,
                            "Dist_q0"=0,"Dist_q1"=0,"Dist_q2"=0,"Dist_qinf"=0,
                            "Rep_q0"=0,"Rep_q1"=0,"Rep_q2"=0,"Rep_qinf"=0,
                            "Eff_q0"=0,"Eff_q1"=0,"Eff_q2"=0,"Eff_qinf"=0)



# redundancy 2016 ----

red_q0.t2 <- raw_sub_rho(meta.2016, 0)
red_q1.t2 <- raw_sub_rho(meta.2016, 1)
red_q2.t2 <- raw_sub_rho(meta.2016, 2)
red_qinf.t2 <- raw_sub_rho(meta.2016, Inf)

measures_2016$Red_q0 <- red_q0.t2$diversity
measures_2016$Red_q1 <- red_q1.t2$diversity
measures_2016$Red_q2 <- red_q2.t2$diversity
measures_2016$Red_qinf <- red_qinf.t2$diversity

# distinctiveness 2016 ----

dist_q0.t2 <- raw_sub_beta(meta.2016, 0)
dist_q1.t2 <- raw_sub_beta(meta.2016, 1)
dist_q2.t2 <- raw_sub_beta(meta.2016, 2)
dist_qinf.t2 <- raw_sub_beta(meta.2016, Inf)

measures_2016$Dist_q0 <- dist_q0.t2$diversity
measures_2016$Dist_q1 <- dist_q1.t2$diversity
measures_2016$Dist_q2 <- dist_q2.t2$diversity
measures_2016$Dist_qinf <- dist_qinf.t2$diversity

# representativeness 2016 ----

rep_q0.t2 <- norm_sub_rho(meta.2016, 0)
rep_q1.t2 <- norm_sub_rho(meta.2016, 1)
rep_q2.t2 <- norm_sub_rho(meta.2016, 2)
rep_qinf.t2 <- norm_sub_rho(meta.2016, Inf)

measures_2016$Rep_q0 <- rep_q0.t2$diversity
measures_2016$Rep_q1 <- rep_q1.t2$diversity
measures_2016$Rep_q2 <- rep_q2.t2$diversity
measures_2016$Rep_qinf <- rep_qinf.t2$diversity

#effective number of distinct subcommunities  2016 ----

eff_q0.t2 <- norm_sub_beta(meta.2016, 0)
eff_q1.t2 <- norm_sub_beta(meta.2016, 1)
eff_q2.t2 <- norm_sub_beta(meta.2016, 2)
eff_qinf.t2 <- norm_sub_beta(meta.2016, Inf)

measures_2016$Eff_q0 <- eff_q0.t2$diversity
measures_2016$Eff_q1 <- eff_q1.t2$diversity
measures_2016$Eff_q2 <- eff_q2.t2$diversity
measures_2016$Eff_qinf <- eff_qinf.t2$diversity



transects_2001 <- measures_2001 %>%
  merge(lc.01.23km, by="RouteName") %>%                                         #adds landcover    
  merge(bcr, by="RouteName") %>%                                                #adds BCRs
  mutate("forest.t1"=deciduousforest.t1 + evergreenforest.t1 + mixedforest.t1)  #creates new forest column

transects_2016 <- measures_2016 %>%  
  merge(lc.16.23km, by="RouteName") %>% 
  merge(bcr, by="RouteName")  %>% 
  mutate("forest.t2"=deciduousforest.t2 + evergreenforest.t2 + mixedforest.t2)



#errors!!!
#measures_2001 <- merge(measures_2001, data.transects.4km, by="RouteName")
#measures_2001 <- mutate(measures_2001, delta.total=(abs(delta.urban) + abs(delta.forest) +abs(delta.cropland)+abs(delta.grassland) +abs(delta.wetland))/2)
#measures_2001 <- merge(measures_2001, bcr, by="RouteName")


#measures_2016 <- merge(measures_2016, data.transects.4km, by="RouteName")
#measures_2016 <- mutate(measures_2016, delta.total=(abs(delta.urban) + abs(delta.forest) +abs(delta.cropland)+abs(delta.grassland) +abs(delta.wetland))/2)
#measures_2016 <- merge(measures_2016, bcr, by="RouteName")



#'######################################################################################
# 2 - P L O T S ----
#
#'######################################################################################

ggplot(measures_2001, aes(as.factor(BCR),Eff_q0)) +
  geom_boxplot() +
  ggtitle("Effective number distinct subcommunities by BCR in 2001", "norm_sub_beta; q=0") +
  ylab("Effective number")+ 
  ylim(0,22)+
  geom_hline(yintercept =  mean(measures_2001$Eff_q0), color="red")

ggplot(measures_2016, aes(as.factor(BCR),Eff_q0)) +
  geom_boxplot() +
  ggtitle("Effective number distinct subcommunities by BCR in 2016", "norm_sub_beta; q=0") +
  ylab("Effective number")+ 
  ylim(0,22)+
  geom_hline(yintercept =  mean(measures_2016$Eff_q0), color="red")


#REDUNDANCY ----
#Urban 2001 ----
ggplot(transects_2001, aes(urban.t1, Red_q0)) +
  geom_point(alpha=0.2, col="black") +
  geom_smooth(method="loess", colour="black") +
  geom_rug() + 
  scale_x_log10()+
  ggtitle("Redundancy q=0 as a function of urban in 2001", "transects") +
  ylab("Redundancy")  +
  xlab("log(urban.t1)") +
  ylim(0,1200)

ggplot(transects_2001, aes(urban.t1, Red_q1)) +
  geom_point(alpha=0.2, col="black") +
  geom_smooth(method="loess", colour="black") +
  geom_rug() + 
  scale_x_log10()+
  ggtitle("Redundancy q=1 as a function of urban in 2001", "transects") +
  ylab("Redundancy")  +
  xlab("log(urban.t1)") +
  ylim(0,1000)

ggplot(transects_2001, aes(urban.t1, Red_q2)) +
  geom_point(alpha=0.2, col="black") +
  geom_smooth(method="loess", colour="black") +
  geom_rug() + 
  scale_x_log10()+
  ggtitle("Redundancy q=2 as a function of urban in 2001", "transects") +
  ylab("Redundancy")  +
  xlab("log(urban.t1)") +
  ylim(0,800)

ggplot(transects_2001, aes(urban.t1, Red_qinf)) +
  geom_point(alpha=0.2, col="black") +
  geom_smooth(method="loess", colour="black") +
  geom_rug() + 
  scale_x_log10()+
  ggtitle("Redundancy q=Infas a function of urban in 2001", "transects") +
  ylab("Redundancy")  +
  xlab("log(urban.t1)") +
  ylim(0,300)


#Urban 2016 ----

ggplot(transects_2016, aes(urban.t2, Red_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgray") +
  geom_rug() + 
  ggtitle("Redundancy q=0 as a function of urban in 2016", "transects") +
  ylab("Redundancy") +
  ylim(0,1200)

ggplot(measures_2016, aes(urban.t2, Red_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgray") +
  geom_rug() + 
  ggtitle("Redundancy q=1 as a function of urban in 2016") +
  ylab("Redundancy") +
  ylim(0,1200)

ggplot(measures_2016, aes(urban.t2, Red_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgray") +
  geom_rug() + 
  ggtitle("Redundancy q=2 as a function of urban in 2016") +
  ylab("Redundancy") +
  ylim(0,1200)

ggplot(measures_2016, aes(urban.t2, Red_qinf)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgray") +
  geom_rug() + 
  ggtitle("Redundancy q=inf as a function of urban in 2016") +
  ylab("Redundancy") +
  ylim(0,1200)


#Forest 2001 ----

ggplot(transects_2001, aes(forest.t1, Red_q0)) +
  geom_point(alpha=0.2, col="darkgreen") +
  geom_smooth(method="loess", color="black") +
  geom_rug() + 
  ggtitle("Redundancy q=0 as a function of forest in 2001", "transects") +
  ylab("Redundancy")  +
  ylim(0,1200)

ggplot(transects_2001, aes(forest.t1, Red_q1)) +
  geom_point(alpha=0.2, col="darkgreen") +
  geom_smooth(method="loess", color="black") +
  geom_rug() + 
  ggtitle("Redundancy q=1 as a function of forest in 2001", "transects") +
  ylab("Redundancy")  +
  ylim(0,1000)

ggplot(transects_2001, aes(forest.t1, Red_q2)) +
  geom_point(alpha=0.2, col="darkgreen") +
  geom_smooth(method="loess", color="black") +
  geom_rug() + 
  ggtitle("Redundancy q=2 as a function of forest in 2001", "transects") +
  ylab("Redundancy")  +
  ylim(0,800)

ggplot(transects_2001, aes(forest.t1, Red_qinf)) +
  geom_point(alpha=0.2, col="darkgreen") +
  geom_smooth(method="loess", color="black") +
  geom_rug() + 
  ggtitle("Redundancy q=inf as a function of forest in 2001", "transects") +
  ylab("Redundancy")  +
  ylim(0,300)

#Forest 2016 ----
ggplot(transects_2016, aes(forest.t2, Red_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgreen") +
  geom_rug() + 
  ggtitle("Redundancy q=0 as a function of forest in 2016","transects") +
  ylab("Redundancy")  +
  ylim(0,1200)


ggplot(measures_2016, aes(forest.t2, Red_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgreen") +
  geom_rug() + 
  ggtitle("Redundancy q=1 as a function of forest in 2016") +
  ylab("Redundancy") +
  ylim(0,1200)

ggplot(measures_2016, aes(forest.t2, Red_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgreen") +
  geom_rug() + 
  ggtitle("Redundancy q=2 as a function of forest in 2016") +
  ylab("Redundancy") +
  ylim(0,1200)

ggplot(measures_2016, aes(forest.t2, Red_qinf)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgreen") +
  geom_rug() + 
  ggtitle("Redundancy q=inf as a function of forest in 2016") +
  ylab("Redundancy") +
  ylim(0,1200)


#Grassland 2001 ----
ggplot(transects_2001, aes(grassland.t1, Red_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Redundancy q=0 as a function of grassland in 2001", "transects") +
  ylab("Redundancy")  +
  ylim(0,1200)


ggplot(measures_2001, aes(grassland.t1, Red_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Redundancy q=1 as a function of grassland in 2001") +
  ylab("Redundancy") +
  ylim(0,1200)

ggplot(measures_2001, aes(grassland.t1, Red_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Redundancy q=2 as a function of grassland in 2001") +
  ylab("Redundancy") +
  ylim(0,1200)

ggplot(measures_2001, aes(grassland.t1, Red_qinf)) +
  geom_point() +
  geom_smooth(method="loess",color="lightgreen") +
  geom_rug() + 
  ggtitle("Redundancy q=inf as a function of grassland in 2001") +
  ylab("Redundancy") +
  ylim(0,1200)

# Grassland 2016 ----
ggplot(transects_2016, aes(grassland.t2, Red_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Redundancy q=0 as a function of grassland in 2016", "transects") +
  ylab("Redundancy")  +
  ylim(0,1200)


ggplot(measures_2016, aes(grassland.t2, Red_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Redundancy q=1 as a function of grassland in 2016") +
  ylab("Redundancy") +
  ylim(0,1200)

ggplot(measures_2016, aes(grassland.t2, Red_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Redundancy q=2 as a function of grassland in 2016") +
  ylab("Redundancy") +
  ylim(0,1200)

ggplot(measures_2016, aes(grassland.t2, Red_qinf)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Redundancy q=inf as a function of grassland in 2016") +
  ylab("Redundancy") +
  ylim(0,1200)

#cropland 2001 ----
ggplot(transects_2001, aes(cropland.t1, Red_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Redundancy q=0 as a function of cropland in 2001", "transects") +
  ylab("Redundancy")  +
  ylim(0,1200)


ggplot(measures_2001, aes(cropland.t1, Red_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Redundancy q=1 as a function of cropland in 2001") +
  ylab("Redundancy") +
  ylim(0,1200)

ggplot(measures_2001, aes(cropland.t1, Red_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Redundancy q=2 as a function of cropland in 2001") +
  ylab("Redundancy") +
  ylim(0,1200)

ggplot(measures_2001, aes(cropland.t1, Red_qinf)) +
  geom_point() +
  geom_smooth(method="loess",color="yellow") +
  geom_rug() + 
  ggtitle("Redundancy q=inf as a function of cropland in 2001") +
  ylab("Redundancy") +
  ylim(0,1200)

# cropland 2016 ----
ggplot(transects_2016, aes(cropland.t2, Red_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Redundancy q=0 as a function of cropland in 2016", "transects") +
  ylab("Redundancy")  +
  ylim(0,1200)


ggplot(measures_2016, aes(cropland.t2, Red_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Redundancy q=1 as a function of cropland in 2016") +
  ylab("Redundancy") +
  ylim(0,1200)

ggplot(measures_2016, aes(cropland.t2, Red_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Redundancy q=2 as a function of cropland in 2016") +
  ylab("Redundancy") +
  ylim(0,1200)

ggplot(measures_2016, aes(cropland.t2, Red_qinf)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Redundancy q=inf as a function of cropland in 2016") +
  ylab("Redundancy") +
  ylim(0,1200)

#Wetland 2001 ----
ggplot(transects_2001, aes(wetland.t1, Red_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Redundancy q=0 as a function of wetland in 2001","transects") +
  ylab("Redundancy")  +
  ylim(0,1200)


ggplot(measures_2001, aes(wetland.t1, Red_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Redundancy q=1 as a function of wetland in 2001") +
  ylab("Redundancy") +
  ylim(0,1200)

ggplot(measures_2001, aes(wetland.t1, Red_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Redundancy q=2 as a function of wetland in 2001") +
  ylab("Redundancy") +
  ylim(0,1200)

ggplot(measures_2001, aes(wetland.t1, Red_qinf)) +
  geom_point() +
  geom_smooth(method="loess",color="purple") +
  geom_rug() + 
  ggtitle("Redundancy q=inf as a function of wetland in 2001") +
  ylab("Redundancy") +
  ylim(0,1200)

# wetland 2016 ----
ggplot(transects_2016, aes(wetland.t2, Red_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Redundancy q=0 as a function of wetland in 2016","transects") +
  ylab("Redundancy")  +
  ylim(0,1200)


ggplot(measures_2016, aes(wetland.t2, Red_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Redundancy q=1 as a function of wetland in 2016") +
  ylab("Redundancy") +
  ylim(0,1200)

ggplot(measures_2016, aes(wetland.t2, Red_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Redundancy q=2 as a function of wetland in 2016") +
  ylab("Redundancy") +
  ylim(0,1200)

ggplot(measures_2016, aes(wetland.t2, Red_qinf)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Redundancy q=inf as a function of wetland in 2016") +
  ylab("Redundancy") +
  ylim(0,1200)

#DISTINCTIVENESS ----
#Urban 2001 ----
ggplot(transects_2001, aes(urban.t1, Dist_q0)) +
  geom_point(alpha=0.2, col="black") +
  geom_smooth(method="loess", colour="black") +
  geom_rug() + 
  scale_x_log10()+
  ggtitle("Distinctiveness q=0 as a function of urban in 2001", "transects") +
  ylab("Distinctiveness")  +
  xlab("log(urban.t1)") +
  ylim(0,0.02) 

ggplot(transects_2001, aes(urban.t1, Dist_q1)) +
  geom_point(alpha=0.2, col="black") +
  geom_smooth(method="loess", colour="black") +
  geom_rug() + 
  scale_x_log10()+
  ggtitle("Distinctiveness q=1 as a function of urban in 2001", "transects") +
  ylab("Distinctiveness")  +
  xlab("log(urban.t1)") +
  ylim(0,0.1) 

ggplot(transects_2001, aes(urban.t1, Dist_q2)) +
  geom_point(alpha=0.2, col="black") +
  geom_smooth(method="loess", colour="black") +
  geom_rug() + 
  scale_x_log10()+
  ggtitle("Distinctiveness q=2 as a function of urban in 2001", "transects") +
  ylab("Distinctiveness")  +
  xlab("log(urban.t1)") +
  ylim(0,0.3) 

ggplot(transects_2001, aes(urban.t1, Dist_qinf)) +
  geom_point(alpha=0.2, col="black") +
  geom_smooth(method="loess", colour="black") +
  geom_rug() + 
  scale_x_log10()+
  ggtitle("Distinctiveness q=inf as a function of urban in 2001", "transects") +
  ylab("Distinctiveness")  +
  xlab("log(urban.t1)") +
  ylim(0,1) 


#Urban 2016 ----

ggplot(transects_2016, aes(urban.t2, Dist_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgray") +
  geom_rug() + 
  ggtitle("Distinctiveness q=0 as a function of urban in 2016", "transects") +
  ylab("Distinctiveness") 
  #ylim(0,25)

ggplot(measures_2016, aes(urban.t2, Dist_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgray") +
  geom_rug() + 
  ggtitle("Distinctiveness q=1 as a function of urban in 2016") +
  ylab("Distinctiveness") +
  ylim(0,25)

ggplot(measures_2016, aes(urban.t2, Dist_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgray") +
  geom_rug() + 
  ggtitle("Distinctiveness q=2 as a function of urban in 2016") +
  ylab("Distinctiveness") +
  ylim(0,25)

ggplot(measures_2016, aes(urban.t2, Dist_qinf)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgray") +
  geom_rug() + 
  ggtitle("Distinctiveness q=inf as a function of urban in 2016") +
  ylab("Distinctiveness") +
  ylim(0,25)


#Forest 2001 ----
ggplot(transects_2001, aes(forest.t1, Dist_q0)) +
  geom_point(alpha=0.2, col="darkgreen") +
  geom_smooth(method="loess", color="black") +
  geom_rug() + 
  ggtitle("Distinctiveness q=0 as a function of forest in 2001", "transects") +
  ylab("Distinctiveness")  +
  ylim(0,0.02)

ggplot(transects_2001, aes(forest.t1, Dist_q1)) +
  geom_point(alpha=0.2, col="darkgreen") +
  geom_smooth(method="loess", color="black") +
  geom_rug() + 
  ggtitle("Distinctiveness q=1 as a function of forest in 2001", "transects") +
  ylab("Distinctiveness")  +
  ylim(0,0.1)

ggplot(transects_2001, aes(forest.t1, Dist_q2)) +
  geom_point(alpha=0.2, col="darkgreen") +
  geom_smooth(method="loess", color="black") +
  geom_rug() + 
  ggtitle("Distinctiveness q=2 as a function of forest in 2001", "transects") +
  ylab("Distinctiveness")  +
  ylim(0,0.3)

ggplot(transects_2001, aes(forest.t1, Dist_qinf)) +
  geom_point(alpha=0.2, col="darkgreen") +
  geom_smooth(method="loess", color="black") +
  geom_rug() + 
  ggtitle("Distinctiveness q=inf as a function of forest in 2001", "transects") +
  ylab("Distinctiveness")  +
  ylim(0,1)



#Forest 2016 ----
ggplot(measures_2016, aes(forest.t2, Dist_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgreen") +
  geom_rug() + 
  ggtitle("Distinctiveness q=0 as a function of forest in 2016") +
  ylab("Distinctiveness")  +
  ylim(0,25)


ggplot(measures_2016, aes(forest.t2, Dist_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgreen") +
  geom_rug() + 
  ggtitle("Distinctiveness q=1 as a function of forest in 2016") +
  ylab("Distinctiveness") +
  ylim(0,25)

ggplot(measures_2016, aes(forest.t2, Dist_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgreen") +
  geom_rug() + 
  ggtitle("Distinctiveness q=2 as a function of forest in 2016") +
  ylab("Distinctiveness") +
  ylim(0,25)

ggplot(measures_2016, aes(forest.t2, Dist_qinf)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgreen") +
  geom_rug() + 
  ggtitle("Distinctiveness q=inf as a function of forest in 2016") +
  ylab("Distinctiveness") +
  ylim(0,25)


#Grassland 2001 ----
ggplot(measures_2001, aes(grassland.t1, Dist_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Distinctiveness q=0 as a function of grassland in 2001") +
  ylab("Distinctiveness")  +
  ylim(0,25)


ggplot(measures_2001, aes(grassland.t1, Dist_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Distinctiveness q=1 as a function of grassland in 2001") +
  ylab("Distinctiveness") +
  ylim(0,25)

ggplot(measures_2001, aes(grassland.t1, Dist_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Distinctiveness q=2 as a function of grassland in 2001") +
  ylab("Distinctiveness") +
  ylim(0,25)

ggplot(measures_2001, aes(grassland.t1, Dist_qinf)) +
  geom_point() +
  geom_smooth(method="loess",color="lightgreen") +
  geom_rug() + 
  ggtitle("Distinctiveness q=inf as a function of grassland in 2001") +
  ylab("Distinctiveness") +
  ylim(0,25)

# Grassland 2016 ----
ggplot(measures_2016, aes(grassland.t2, Dist_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Distinctiveness q=0 as a function of grassland in 2016") +
  ylab("Distinctiveness")  +
  ylim(0,25)


ggplot(measures_2016, aes(grassland.t2, Dist_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Distinctiveness q=1 as a function of grassland in 2016") +
  ylab("Distinctiveness") +
  ylim(0,25)

ggplot(measures_2016, aes(grassland.t2, Dist_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Distinctiveness q=2 as a function of grassland in 2016") +
  ylab("Distinctiveness") +
  ylim(0,25)

ggplot(measures_2016, aes(grassland.t2, Dist_qinf)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Distinctiveness q=inf as a function of grassland in 2016") +
  ylab("Distinctiveness") +
  ylim(0,25)

#cropland 2001 ----
ggplot(measures_2001, aes(cropland.t1, Dist_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Distinctiveness q=0 as a function of cropland in 2001") +
  ylab("Distinctiveness")  +
  ylim(0,25)


ggplot(measures_2001, aes(cropland.t1, Dist_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Distinctiveness q=1 as a function of cropland in 2001") +
  ylab("Distinctiveness") +
  ylim(0,25)

ggplot(measures_2001, aes(cropland.t1, Dist_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Distinctiveness q=2 as a function of cropland in 2001") +
  ylab("Distinctiveness") +
  ylim(0,25)

ggplot(measures_2001, aes(cropland.t1, Dist_qinf)) +
  geom_point() +
  geom_smooth(method="loess",color="yellow") +
  geom_rug() + 
  ggtitle("Distinctiveness q=inf as a function of cropland in 2001") +
  ylab("Distinctiveness") +
  ylim(0,25)

# cropland 2016 ----
ggplot(measures_2016, aes(cropland.t2, Dist_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Distinctiveness q=0 as a function of cropland in 2016") +
  ylab("Distinctiveness")  +
  ylim(0,25)


ggplot(measures_2016, aes(cropland.t2, Dist_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Distinctiveness q=1 as a function of cropland in 2016") +
  ylab("Distinctiveness") +
  ylim(0,25)

ggplot(measures_2016, aes(cropland.t2, Dist_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Distinctiveness q=2 as a function of cropland in 2016") +
  ylab("Distinctiveness") +
  ylim(0,25)

ggplot(measures_2016, aes(cropland.t2, Dist_qinf)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Distinctiveness q=inf as a function of cropland in 2016") +
  ylab("Distinctiveness") +
  ylim(0,25)

#Wetland 2001 ----
ggplot(measures_2001, aes(wetland.t1, Dist_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Distinctiveness q=0 as a function of wetland in 2001") +
  ylab("Distinctiveness")  +
  ylim(0,25)


ggplot(measures_2001, aes(wetland.t1, Dist_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Distinctiveness q=1 as a function of wetland in 2001") +
  ylab("Distinctiveness") +
  ylim(0,25)

ggplot(measures_2001, aes(wetland.t1, Dist_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Distinctiveness q=2 as a function of wetland in 2001") +
  ylab("Distinctiveness") +
  ylim(0,25)

ggplot(measures_2001, aes(wetland.t1, Dist_qinf)) +
  geom_point() +
  geom_smooth(method="loess",color="purple") +
  geom_rug() + 
  ggtitle("Distinctiveness q=inf as a function of wetland in 2001") +
  ylab("Distinctiveness") +
  ylim(0,25)

# wetland 2016 ----
ggplot(measures_2016, aes(wetland.t2, Dist_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Distinctiveness q=0 as a function of wetland in 2016") +
  ylab("Distinctiveness")  +
  ylim(0,25)


ggplot(measures_2016, aes(wetland.t2, Dist_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Distinctiveness q=1 as a function of wetland in 2016") +
  ylab("Distinctiveness") +
  ylim(0,25)

ggplot(measures_2016, aes(wetland.t2, Dist_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Distinctiveness q=2 as a function of wetland in 2016") +
  ylab("Distinctiveness") +
  ylim(0,25)

ggplot(measures_2016, aes(wetland.t2, Dist_qinf)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Distinctiveness q=inf as a function of wetland in 2016") +
  ylab("Distinctiveness") +
  ylim(0,25)

#REPRESENTATIVENESS ----
#Urban 2001 ----
ggplot(transects_2001, aes(urban.t1, Rep_q0)) +
  geom_point(alpha=0.2, col="black") +
  geom_smooth(method="loess", colour="black") +
  geom_rug() + 
  scale_x_log10()+
  ggtitle("Representativeness q=0 as a function of urban in 2001", "transects") +
  ylab("Representativeness")  +
  xlab("log(urban.t1)") +
  ylim(0,1) 

ggplot(transects_2001, aes(urban.t1, Rep_q1)) +
  geom_point(alpha=0.2, col="black") +
  geom_smooth(method="loess", colour="black") +
  geom_rug() + 
  scale_x_log10()+
  ggtitle("Representativeness q=1 as a function of urban in 2001", "transects") +
  ylab("Representativeness")  +
  xlab("log(urban.t1)") +
  ylim(0,0.7) 

ggplot(transects_2001, aes(urban.t1, Rep_q2)) +
  geom_point(alpha=0.2, col="black") +
  geom_smooth(method="loess", colour="black") +
  geom_rug() + 
  scale_x_log10()+
  ggtitle("Representativeness q=2 as a function of urban in 2001", "transects") +
  ylab("Representativeness")  +
  xlab("log(urban.t1)") +
  ylim(0,0.6) 

ggplot(transects_2001, aes(urban.t1, Rep_qinf)) +
  geom_point(alpha=0.2, col="black") +
  geom_smooth(method="loess", colour="black") +
  geom_rug() + 
  scale_x_log10()+
  ggtitle("Representativeness q=inf as a function of urban in 2001", "transects") +
  ylab("Representativeness")  +
  xlab("log(urban.t1)") +
  ylim(0,0.3) 



#Urban 2016 ----

ggplot(transects_2016, aes(urban.t2, Rep_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgray") +
  geom_rug() + 
  ggtitle("Representativeness q=0 as a function of urban in 2016", "transects") +
  ylab("Representativeness") +
  ylim(0,0.8)

ggplot(measures_2016, aes(urban.t2, Rep_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgray") +
  geom_rug() + 
  ggtitle("Representativeness q=1 as a function of urban in 2016") +
  ylab("Representativeness") +
  ylim(0,0.8)

ggplot(measures_2016, aes(urban.t2, Rep_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgray") +
  geom_rug() + 
  ggtitle("Representativeness q=2 as a function of urban in 2016") +
  ylab("Representativeness") +
  ylim(0,0.8)

ggplot(measures_2016, aes(urban.t2, Rep_qinf)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgray") +
  geom_rug() + 
  ggtitle("Representativeness q=inf as a function of urban in 2016") +
  ylab("Representativeness") +
  ylim(0,0.8)


#Forest 2001 ----
ggplot(transects_2001, aes(forest.t1, Rep_q0)) +
  geom_point(alpha=0.2, col="darkgreen") +
  geom_smooth(method="loess", color="black") +
  geom_rug() + 
  ggtitle("Representativeness q=0 as a function of forest in 2001", "transects") +
  ylab("Representativeness")  +
  ylim(0,1)

ggplot(transects_2001, aes(forest.t1, Rep_q1)) +
  geom_point(alpha=0.2, col="darkgreen") +
  geom_smooth(method="loess", color="black") +
  geom_rug() + 
  ggtitle("Representativeness q=1 as a function of forest in 2001", "transects") +
  ylab("Representativeness")  +
  ylim(0,0.7)

ggplot(transects_2001, aes(forest.t1, Rep_q2)) +
  geom_point(alpha=0.2, col="darkgreen") +
  geom_smooth(method="loess", color="black") +
  geom_rug() + 
  ggtitle("Representativeness q=2 as a function of forest in 2001", "transects") +
  ylab("Representativeness")  +
  ylim(0,0.6)

ggplot(transects_2001, aes(forest.t1, Rep_qinf)) +
  geom_point(alpha=0.2, col="darkgreen") +
  geom_smooth(method="loess", color="black") +
  geom_rug() + 
  ggtitle("Representativeness q=inf as a function of forest in 2001", "transects") +
  ylab("Representativeness")  +
  ylim(0,0.3)



#Forest 2016 ----
ggplot(measures_2016, aes(forest.t2, Rep_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgreen") +
  geom_rug() + 
  ggtitle("Representativeness q=0 as a function of forest in 2016") +
  ylab("Representativeness")  +
  ylim(0,0.8)


ggplot(measures_2016, aes(forest.t2, Rep_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgreen") +
  geom_rug() + 
  ggtitle("Representativeness q=1 as a function of forest in 2016") +
  ylab("Representativeness") +
  ylim(0,0.8)

ggplot(measures_2016, aes(forest.t2, Rep_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgreen") +
  geom_rug() + 
  ggtitle("Representativeness q=2 as a function of forest in 2016") +
  ylab("Representativeness") +
  ylim(0,0.8)

ggplot(measures_2016, aes(forest.t2, Rep_qinf)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgreen") +
  geom_rug() + 
  ggtitle("Representativeness q=inf as a function of forest in 2016") +
  ylab("Representativeness") +
  ylim(0,0.8)


#Grassland 2001 ----
ggplot(measures_2001, aes(grassland.t1, Rep_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Representativeness q=0 as a function of grassland in 2001") +
  ylab("Representativeness")  +
  ylim(0,0.8)


ggplot(measures_2001, aes(grassland.t1, Rep_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Representativeness q=1 as a function of grassland in 2001") +
  ylab("Representativeness") +
  ylim(0,0.8)

ggplot(measures_2001, aes(grassland.t1, Rep_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Representativeness q=2 as a function of grassland in 2001") +
  ylab("Representativeness") +
  ylim(0,0.8)

ggplot(measures_2001, aes(grassland.t1, Rep_qinf)) +
  geom_point() +
  geom_smooth(method="loess",color="lightgreen") +
  geom_rug() + 
  ggtitle("Representativeness q=inf as a function of grassland in 2001") +
  ylab("Representativeness") +
  ylim(0,0.8)

# Grassland 2016 ----
ggplot(measures_2016, aes(grassland.t2, Rep_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Representativeness q=0 as a function of grassland in 2016") +
  ylab("Representativeness")  +
  ylim(0,0.8)


ggplot(measures_2016, aes(grassland.t2, Rep_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Representativeness q=1 as a function of grassland in 2016") +
  ylab("Representativeness") +
  ylim(0,0.8)

ggplot(measures_2016, aes(grassland.t2, Rep_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Representativeness q=2 as a function of grassland in 2016") +
  ylab("Representativeness") +
  ylim(0,0.8)

ggplot(measures_2016, aes(grassland.t2, Rep_qinf)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Representativeness q=inf as a function of grassland in 2016") +
  ylab("Representativeness") +
  ylim(0,0.8)

#cropland 2001 ----
ggplot(measures_2001, aes(cropland.t1, Rep_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Representativeness q=0 as a function of cropland in 2001") +
  ylab("Representativeness")  +
  ylim(0,0.8)


ggplot(measures_2001, aes(cropland.t1, Rep_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Representativeness q=1 as a function of cropland in 2001") +
  ylab("Representativeness") +
  ylim(0,0.8)

ggplot(measures_2001, aes(cropland.t1, Rep_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Representativeness q=2 as a function of cropland in 2001") +
  ylab("Representativeness") +
  ylim(0,0.8)

ggplot(measures_2001, aes(cropland.t1, Rep_qinf)) +
  geom_point() +
  geom_smooth(method="loess",color="yellow") +
  geom_rug() + 
  ggtitle("Representativeness q=inf as a function of cropland in 2001") +
  ylab("Representativeness") +
  ylim(0,0.8)

# cropland 2016 ----
ggplot(measures_2016, aes(cropland.t2, Rep_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Representativeness q=0 as a function of cropland in 2016") +
  ylab("Representativeness")  +
  ylim(0,0.8)


ggplot(measures_2016, aes(cropland.t2, Rep_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Representativeness q=1 as a function of cropland in 2016") +
  ylab("Representativeness") +
  ylim(0,0.8)

ggplot(measures_2016, aes(cropland.t2, Rep_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Representativeness q=2 as a function of cropland in 2016") +
  ylab("Representativeness") +
  ylim(0,0.8)

ggplot(measures_2016, aes(cropland.t2, Rep_qinf)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Representativeness q=inf as a function of cropland in 2016") +
  ylab("Representativeness") +
  ylim(0,0.8)

#Wetland 2001 ----
ggplot(measures_2001, aes(wetland.t1, Rep_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Representativeness q=0 as a function of wetland in 2001") +
  ylab("Representativeness")  +
  ylim(0,0.8)


ggplot(measures_2001, aes(wetland.t1, Rep_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Representativeness q=1 as a function of wetland in 2001") +
  ylab("Representativeness") +
  ylim(0,0.8)

ggplot(measures_2001, aes(wetland.t1, Rep_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Representativeness q=2 as a function of wetland in 2001") +
  ylab("Representativeness") +
  ylim(0,0.8)

ggplot(measures_2001, aes(wetland.t1, Rep_qinf)) +
  geom_point() +
  geom_smooth(method="loess",color="purple") +
  geom_rug() + 
  ggtitle("Representativeness q=inf as a function of wetland in 2001") +
  ylab("Representativeness") +
  ylim(0,0.8)

# wetland 2016 ----
ggplot(measures_2016, aes(wetland.t2, Rep_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Representativeness q=0 as a function of wetland in 2016") +
  ylab("Representativeness")  +
  ylim(0,0.8)


ggplot(measures_2016, aes(wetland.t2, Rep_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Representativeness q=1 as a function of wetland in 2016") +
  ylab("Representativeness") +
  ylim(0,0.8)

ggplot(measures_2016, aes(wetland.t2, Rep_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Representativeness q=2 as a function of wetland in 2016") +
  ylab("Representativeness") +
  ylim(0,0.8)

ggplot(measures_2016, aes(wetland.t2, Rep_qinf)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Representativeness q=inf as a function of wetland in 2016") +
  ylab("Representativeness") +
  ylim(0,0.8)

#EFFECTIVE NUMBER OF DISTINCT SUBCOMMUNITIES ----
#Urban 2001 ----
ggplot(transects_2001, aes(urban.t1, Eff_q0)) +
  geom_point(alpha=0.2, col="black") +
  geom_smooth(method="loess", colour="black") +
  geom_rug() + 
  scale_x_log10()+
  ggtitle("Eff. nr. distinct communitites q=0 as a function of urban in 2001", "transects") +
  ylab("Eff. nr. distinct communitites")  +
  xlab("log(urban.t1)") +
  ylim(0,15) 

ggplot(transects_2001, aes(urban.t1, Eff_q1)) +
  geom_point(alpha=0.2, col="black") +
  geom_smooth(method="loess", colour="black") +
  geom_rug() + 
  scale_x_log10()+
  ggtitle("Eff. nr. distinct communitites q=1 as a function of urban in 2001", "transects") +
  ylab("Eff. nr. distinct communitites")  +
  xlab("log(urban.t1)") +
  ylim(0,55) 

ggplot(transects_2001, aes(urban.t1, Eff_q2)) +
  geom_point(alpha=0.2, col="black") +
  geom_smooth(method="loess", colour="black") +
  geom_rug() + 
  scale_x_log10()+
  ggtitle("Eff. nr. distinct communitites q=2 as a function of urban in 2001", "transects") +
  ylab("Eff. nr. distinct communitites")  +
  xlab("log(urban.t1)") +
  ylim(0,150) 

ggplot(transects_2001, aes(urban.t1, Eff_qinf)) +
  geom_point(alpha=0.2, col="black") +
  geom_smooth(method="loess", colour="black") +
  geom_rug() + 
  scale_x_log10()+
  ggtitle("Eff. nr. distinct communitites q=inf as a function of urban in 2001", "transects") +
  ylab("Eff. nr. distinct communitites")  +
  xlab("log(urban.t1)") +
  ylim(0,2000) 

#Urban 2016 ----

ggplot(transects_2016, aes(urban.t2, Eff_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgray") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=0 as a function of urban in 2016", "transects") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)

ggplot(measures_2016, aes(urban.t2, Eff_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgray") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=1 as a function of urban in 2016") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)

ggplot(measures_2016, aes(urban.t2, Eff_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgray") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=2 as a function of urban in 2016") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)

ggplot(measures_2016, aes(urban.t2, Eff_qinf)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgray") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=inf as a function of urban in 2016") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)


#Forest 2001 ----
ggplot(transects_2001, aes(forest.t1, Eff_q0)) +
  geom_point(alpha=0.2, col="darkgreen") +
  geom_smooth(method="loess", color="black") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=0 as a function of forest in 2001", "transects") +
  ylab("Eff. nr. distinct communitites")  +
  ylim(0,15)

ggplot(transects_2001, aes(forest.t1, Eff_q1)) +
  geom_point(alpha=0.2, col="darkgreen") +
  geom_smooth(method="loess", color="black") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=1 as a function of forest in 2001", "transects") +
  ylab("Eff. nr. distinct communitites")  +
  ylim(0,55)

ggplot(transects_2001, aes(forest.t1, Eff_q2)) +
  geom_point(alpha=0.2, col="darkgreen") +
  geom_smooth(method="loess", color="black") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=2 as a function of forest in 2001", "transects") +
  ylab("Eff. nr. distinct communitites")  +
  ylim(0,150)

ggplot(transects_2001, aes(forest.t1, Eff_qinf)) +
  geom_point(alpha=0.2, col="darkgreen") +
  geom_smooth(method="loess", color="black") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q= inf as a function of forest in 2001", "transects") +
  ylab("Eff. nr. distinct communitites")  +
  ylim(0,2000)



#Forest 2016 ----
ggplot(measures_2016, aes(forest.t2, Eff_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgreen") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=0 as a function of forest in 2016") +
  ylab("Eff. nr. distinct communitites")  +
  ylim(0,25)


ggplot(measures_2016, aes(forest.t2, Eff_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgreen") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=1 as a function of forest in 2016") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)

ggplot(measures_2016, aes(forest.t2, Eff_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgreen") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=2 as a function of forest in 2016") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)

ggplot(measures_2016, aes(forest.t2, Eff_qinf)) +
  geom_point() +
  geom_smooth(method="loess", color="darkgreen") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=inf as a function of forest in 2016") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)


#Grassland 2001 ----
ggplot(measures_2001, aes(grassland.t1, Eff_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=0 as a function of grassland in 2001") +
  ylab("Eff. nr. distinct communitites")  +
  ylim(0,25)


ggplot(measures_2001, aes(grassland.t1, Eff_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=1 as a function of grassland in 2001") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)

ggplot(measures_2001, aes(grassland.t1, Eff_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=2 as a function of grassland in 2001") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)

ggplot(measures_2001, aes(grassland.t1, Eff_qinf)) +
  geom_point() +
  geom_smooth(method="loess",color="lightgreen") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=inf as a function of grassland in 2001") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)

# Grassland 2016 ----
ggplot(measures_2016, aes(grassland.t2, Eff_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=0 as a function of grassland in 2016") +
  ylab("Eff. nr. distinct communitites")  +
  ylim(0,25)


ggplot(measures_2016, aes(grassland.t2, Eff_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=1 as a function of grassland in 2016") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)

ggplot(measures_2016, aes(grassland.t2, Eff_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=2 as a function of grassland in 2016") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)

ggplot(measures_2016, aes(grassland.t2, Eff_qinf)) +
  geom_point() +
  geom_smooth(method="loess", color="lightgreen") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=inf as a function of grassland in 2016") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)

#cropland 2001 ----
ggplot(measures_2001, aes(cropland.t1, Eff_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=0 as a function of cropland in 2001") +
  ylab("Eff. nr. distinct communitites")  +
  ylim(0,25)


ggplot(measures_2001, aes(cropland.t1, Eff_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=1 as a function of cropland in 2001") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)

ggplot(measures_2001, aes(cropland.t1, Eff_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=2 as a function of cropland in 2001") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)

ggplot(measures_2001, aes(cropland.t1, Eff_qinf)) +
  geom_point() +
  geom_smooth(method="loess",color="yellow") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=inf as a function of cropland in 2001") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)

# cropland 2016 ----
ggplot(measures_2016, aes(cropland.t2, Eff_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=0 as a function of cropland in 2016") +
  ylab("Eff. nr. distinct communitites")  +
  ylim(0,25)


ggplot(measures_2016, aes(cropland.t2, Eff_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=1 as a function of cropland in 2016") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)

ggplot(measures_2016, aes(cropland.t2, Eff_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=2 as a function of cropland in 2016") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)

ggplot(measures_2016, aes(cropland.t2, Eff_qinf)) +
  geom_point() +
  geom_smooth(method="loess", color="yellow") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=inf as a function of cropland in 2016") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)

#Wetland 2001 ----
ggplot(measures_2001, aes(wetland.t1, Eff_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=0 as a function of wetland in 2001") +
  ylab("Eff. nr. distinct communitites")  +
  ylim(0,25)


ggplot(measures_2001, aes(wetland.t1, Eff_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=1 as a function of wetland in 2001") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)

ggplot(measures_2001, aes(wetland.t1, Eff_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=2 as a function of wetland in 2001") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)

ggplot(measures_2001, aes(wetland.t1, Eff_qinf)) +
  geom_point() +
  geom_smooth(method="loess",color="purple") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=inf as a function of wetland in 2001") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)

# wetland 2016 ----
ggplot(measures_2016, aes(wetland.t2, Eff_q0)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=0 as a function of wetland in 2016") +
  ylab("Eff. nr. distinct communitites")  +
  ylim(0,25)


ggplot(measures_2016, aes(wetland.t2, Eff_q1)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=1 as a function of wetland in 2016") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)

ggplot(measures_2016, aes(wetland.t2, Eff_q2)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=2 as a function of wetland in 2016") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)

ggplot(measures_2016, aes(wetland.t2, Eff_qinf)) +
  geom_point() +
  geom_smooth(method="loess", color="purple") +
  geom_rug() + 
  ggtitle("Eff. nr. distinct communitites q=inf as a function of wetland in 2016") +
  ylab("Eff. nr. distinct communitites") +
  ylim(0,25)

#' ########################################################################################
# 3 - Distinctiveness by species by SEGMENTS -----
#' 
#' ########################################################################################


meta.sp.2001 <- metacommunity(t(spmatrix.t1.full))
meta.sp.2016 <- metacommunity(t(spmatrix.t2.full))

# 2001 ----

segments_2001 <- data.frame("partition"=rownames(spmatrix.t1.full),
                        "Red_q0"=0,"Red_q1"=0,"Red_q2"=0,"Red_qinf"=0,
                        "Dist_q0"=0,"Dist_q1"=0,"Dist_q2"=0,"Dist_qinf"=0,
                        "Rep_q0"=0,"Rep_q1"=0,"Rep_q2"=0,"Rep_qinf"=0,
                        "Eff_q0"=0,"Eff_q1"=0,"Eff_q2"=0,"Eff_qinf"=0)

# A) redundancy 2001 ----

red_q0.t1 <- raw_sub_rho(meta.sp.2001, 0)
red_q1.t1 <- raw_sub_rho(meta.sp.2001, 1)
red_q2.t1 <- raw_sub_rho(meta.sp.2001, 2)
red_qinf.t1 <- raw_sub_rho(meta.sp.2001, Inf)

segments_2001$Red_q0 <- red_q0.t1$diversity
segments_2001$Red_q1 <- red_q1.t1$diversity
segments_2001$Red_q2 <- red_q2.t1$diversity
segments_2001$Red_qinf <- red_qinf.t1$diversity

# B) distinctiveness 2001 ----

dist_q0.t1 <- raw_sub_beta(meta.sp.2001, 0)
dist_q1.t1 <- raw_sub_beta(meta.sp.2001, 1)
dist_q2.t1 <- raw_sub_beta(meta.sp.2001, 2)
dist_qinf.t1 <- raw_sub_beta(meta.sp.2001, Inf)

segments_2001$Dist_q0 <- dist_q0.t1$diversity
segments_2001$Dist_q1 <- dist_q1.t1$diversity
segments_2001$Dist_q2 <- dist_q2.t1$diversity
segments_2001$Dist_qinf <- dist_qinf.t1$diversity

# C) representativeness 2001 ----

rep_q0.t1 <- norm_sub_rho(meta.sp.2001, 0)
rep_q1.t1 <- norm_sub_rho(meta.sp.2001, 1)
rep_q2.t1 <- norm_sub_rho(meta.sp.2001, 2)
rep_qinf.t1 <- norm_sub_rho(meta.sp.2001, Inf)

segments_2001$Rep_q0 <- rep_q0.t1$diversity
segments_2001$Rep_q1 <- rep_q1.t1$diversity
segments_2001$Rep_q2 <- rep_q2.t1$diversity
segments_2001$Rep_qinf <- rep_qinf.t1$diversity

# D) effective number of distinct subcommunities 2001 ----

eff_q0.t1 <- norm_sub_beta(meta.sp.2001, 0)
eff_q1.t1 <- norm_sub_beta(meta.sp.2001, 1)
eff_q2.t1 <- norm_sub_beta(meta.sp.2001, 2)
eff_qinf.t1 <- norm_sub_beta(meta.sp.2001, Inf)

segments_2001$Eff_q0 <- eff_q0.t1$diversity
segments_2001$Eff_q1 <- eff_q1.t1$diversity
segments_2001$Eff_q2 <- eff_q2.t1$diversity
segments_2001$Eff_qinf <- eff_qinf.t1$diversity

segments_2001 <- rename_at(.tbl = segments_2001, vars(Red_q0:Eff_qinf), paste0,".sp")

# 2016 ----

segments_2016 <- data.frame("partition"=rownames(spmatrix.t2.full),
                            "Red_q0"=0,"Red_q1"=0,"Red_q2"=0,"Red_qinf"=0,
                            "Dist_q0"=0,"Dist_q1"=0,"Dist_q2"=0,"Dist_qinf"=0,
                            "Rep_q0"=0,"Rep_q1"=0,"Rep_q2"=0,"Rep_qinf"=0,
                            "Eff_q0"=0,"Eff_q1"=0,"Eff_q2"=0,"Eff_qinf"=0)

# A) redundancy 2016 ----

red_q0.t2 <- raw_sub_rho(meta.sp.2016, 0)
red_q1.t2 <- raw_sub_rho(meta.sp.2016, 1)
red_q2.t2 <- raw_sub_rho(meta.sp.2016, 2)
red_qinf.t2 <- raw_sub_rho(meta.sp.2016, Inf)

segments_2016$Red_q0 <- red_q0.t2$diversity
segments_2016$Red_q1 <- red_q1.t2$diversity
segments_2016$Red_q2 <- red_q2.t2$diversity
segments_2016$Red_qinf <- red_qinf.t2$diversity

# B) distinctiveness 2016 ----

dist_q0.t2 <- raw_sub_beta(meta.sp.2016, 0)
dist_q1.t2 <- raw_sub_beta(meta.sp.2016, 1)
dist_q2.t2 <- raw_sub_beta(meta.sp.2016, 2)
dist_qinf.t2 <- raw_sub_beta(meta.sp.2016, Inf)

segments_2016$Dist_q0 <- dist_q0.t2$diversity
segments_2016$Dist_q1 <- dist_q1.t2$diversity
segments_2016$Dist_q2 <- dist_q2.t2$diversity
segments_2016$Dist_qinf <- dist_qinf.t2$diversity

# C) representativeness 2001 ----

rep_q0.t2 <- norm_sub_rho(meta.sp.2016, 0)
rep_q1.t2 <- norm_sub_rho(meta.sp.2016, 1)
rep_q2.t2 <- norm_sub_rho(meta.sp.2016, 2)
rep_qinf.t2 <- norm_sub_rho(meta.sp.2016, Inf)

segments_2016$Rep_q0 <- rep_q0.t2$diversity
segments_2016$Rep_q1 <- rep_q1.t2$diversity
segments_2016$Rep_q2 <- rep_q2.t2$diversity
segments_2016$Rep_qinf <- rep_qinf.t2$diversity

# D) effective number of distinct subcommunities 2001 ----

eff_q0.t2 <- norm_sub_beta(meta.sp.2016, 0)
eff_q1.t2 <- norm_sub_beta(meta.sp.2016, 1)
eff_q2.t2 <- norm_sub_beta(meta.sp.2016, 2)
eff_qinf.t2 <- norm_sub_beta(meta.sp.2016, Inf)

segments_2016$Eff_q0 <- eff_q0.t2$diversity
segments_2016$Eff_q1 <- eff_q1.t2$diversity
segments_2016$Eff_q2 <- eff_q2.t2$diversity
segments_2016$Eff_qinf <- eff_qinf.t2$diversity

segments_2016 <- rename_at(.tbl = segments_2016, vars(Red_q0:Eff_qinf), paste0,".sp")









#' ########################################################################################
# 4 - Distinctiveness by landscape by SEGMENTS -----
#' 
#' ########################################################################################

landscapes.2001 <- data.segments.4km %>%  select(partition,urban.t1, forest.t1, grassland.t1,cropland.t1, wetland.t1) 
rownames(landscapes.2001) <- landscapes.2001$partition
landscapes.2001 <- landscapes.2001 %>%  select(-partition)

landscapes.2016 <- data.segments.4km %>%  select(partition,urban.t2, forest.t2, grassland.t2,cropland.t2, wetland.t2) 
rownames(landscapes.2016) <- landscapes.2016$partition
landscapes.2016 <- landscapes.2016 %>%  select(-partition)

meta.land.2001 <- metacommunity(t(landscapes.2001))
meta.land.2016 <- metacommunity(t(landscapes.2016))

#  2001 ----

land_2001 <- data.frame("partition"=rownames(landscapes.2001),
                        "Red_q0"=0,"Red_q1"=0,"Red_q2"=0,"Red_qinf"=0,
                        "Dist_q0"=0,"Dist_q1"=0,"Dist_q2"=0,"Dist_qinf"=0,
                        "Rep_q0"=0,"Rep_q1"=0,"Rep_q2"=0,"Rep_qinf"=0,
                        "Eff_q0"=0,"Eff_q1"=0,"Eff_q2"=0,"Eff_qinf"=0)

# A) redundancy 2001 ----

red_q0.t1 <- raw_sub_rho(meta.land.2001, 0)
red_q1.t1 <- raw_sub_rho(meta.land.2001, 1)
red_q2.t1 <- raw_sub_rho(meta.land.2001, 2)
red_qinf.t1 <- raw_sub_rho(meta.land.2001, Inf)

land_2001$Red_q0 <- red_q0.t1$diversity
land_2001$Red_q1 <- red_q1.t1$diversity
land_2001$Red_q2 <- red_q2.t1$diversity
land_2001$Red_qinf <- red_qinf.t1$diversity

# B) distinctiveness 2001 ----

dist_q0.t1 <- raw_sub_beta(meta.land.2001, 0)
dist_q1.t1 <- raw_sub_beta(meta.land.2001, 1)
dist_q2.t1 <- raw_sub_beta(meta.land.2001, 2)
dist_qinf.t1 <- raw_sub_beta(meta.land.2001, Inf)

land_2001$Dist_q0 <- dist_q0.t1$diversity
land_2001$Dist_q1 <- dist_q1.t1$diversity
land_2001$Dist_q2 <- dist_q2.t1$diversity
land_2001$Dist_qinf <- dist_qinf.t1$diversity

# C) representativeness 2001 ----

rep_q0.t1 <- norm_sub_rho(meta.land.2001, 0)
rep_q1.t1 <- norm_sub_rho(meta.land.2001, 1)
rep_q2.t1 <- norm_sub_rho(meta.land.2001, 2)
rep_qinf.t1 <- norm_sub_rho(meta.land.2001, Inf)

land_2001$Rep_q0 <- rep_q0.t1$diversity
land_2001$Rep_q1 <- rep_q1.t1$diversity
land_2001$Rep_q2 <- rep_q2.t1$diversity
land_2001$Rep_qinf <- rep_qinf.t1$diversity

# D) effective number of distinct subcommunities 2001 ----

eff_q0.t1 <- norm_sub_beta(meta.land.2001, 0)
eff_q1.t1 <- norm_sub_beta(meta.land.2001, 1)
eff_q2.t1 <- norm_sub_beta(meta.land.2001, 2)
eff_qinf.t1 <- norm_sub_beta(meta.land.2001, Inf)

land_2001$Eff_q0 <- eff_q0.t1$diversity
land_2001$Eff_q1 <- eff_q1.t1$diversity
land_2001$Eff_q2 <- eff_q2.t1$diversity
land_2001$Eff_qinf <- eff_qinf.t1$diversity

land_2001 <- rename_at(.tbl = land_2001, vars(Red_q0:Eff_qinf), paste0,".L")
land_sp_2001 <- merge(land_2001, segments_2001, by="partition")



# 2016 ----

land_2016 <- data.frame("partition"=rownames(spmatrix.t2.full),
                        "Red_q0"=0,"Red_q1"=0,"Red_q2"=0,"Red_qinf"=0,
                        "Dist_q0"=0,"Dist_q1"=0,"Dist_q2"=0,"Dist_qinf"=0,
                        "Rep_q0"=0,"Rep_q1"=0,"Rep_q2"=0,"Rep_qinf"=0,
                        "Eff_q0"=0,"Eff_q1"=0,"Eff_q2"=0,"Eff_qinf"=0)

# A) redundancy 2016 ----

red_q0.t2 <- raw_sub_rho(meta.land.2016, 0)
red_q1.t2 <- raw_sub_rho(meta.land.2016, 1)
red_q2.t2 <- raw_sub_rho(meta.land.2016, 2)
red_qinf.t2 <- raw_sub_rho(meta.land.2016, Inf)

land_2016$Red_q0 <- red_q0.t2$diversity
land_2016$Red_q1 <- red_q1.t2$diversity
land_2016$Red_q2 <- red_q2.t2$diversity
land_2016$Red_qinf <- red_qinf.t2$diversity

# B) distinctiveness 2016 ----

dist_q0.t2 <- raw_sub_beta(meta.land.2016, 0)
dist_q1.t2 <- raw_sub_beta(meta.land.2016, 1)
dist_q2.t2 <- raw_sub_beta(meta.land.2016, 2)
dist_qinf.t2 <- raw_sub_beta(meta.land.2016, Inf)

land_2016$Dist_q0 <- dist_q0.t2$diversity
land_2016$Dist_q1 <- dist_q1.t2$diversity
land_2016$Dist_q2 <- dist_q2.t2$diversity
land_2016$Dist_qinf <- dist_qinf.t2$diversity

# C) representativeness 2016 ----

rep_q0.t2 <- norm_sub_rho(meta.land.2016, 0)
rep_q1.t2 <- norm_sub_rho(meta.land.2016, 1)
rep_q2.t2 <- norm_sub_rho(meta.land.2016, 2)
rep_qinf.t2 <- norm_sub_rho(meta.land.2016, Inf)

land_2016$Rep_q0 <- rep_q0.t2$diversity
land_2016$Rep_q1 <- rep_q1.t2$diversity
land_2016$Rep_q2 <- rep_q2.t2$diversity
land_2016$Rep_qinf <- rep_qinf.t2$diversity

# D) effective number of distinct subcommunities 2016 ----

eff_q0.t2 <- norm_sub_beta(meta.land.2016, 0)
eff_q1.t2 <- norm_sub_beta(meta.land.2016, 1)
eff_q2.t2 <- norm_sub_beta(meta.land.2016, 2)
eff_qinf.t2 <- norm_sub_beta(meta.land.2016, Inf)

land_2016$Eff_q0 <- eff_q0.t2$diversity
land_2016$Eff_q1 <- eff_q1.t2$diversity
land_2016$Eff_q2 <- eff_q2.t2$diversity
land_2016$Eff_qinf <- eff_qinf.t2$diversity

land_2016 <- rename_at(.tbl = land_2016, vars(Red_q0:Eff_qinf), paste0,".L")
land_sp_2016 <- merge(land_2016, segments_2016, by="partition")


# P L O T S ----

#Redundancy
ggplot(land_2001, aes(Red_q0.L, Red_q0)) +
  geom_point()+
  geom_smooth(method = "glm") +
  geom_rug() + 
  scale_x_log10()+
  ylab("Redundancy of sp community") +
  xlab("Redundancy of landscape") +
  ylim(c(-100,1200)) +
  ggtitle("Redundancy of community as a function of distinctiveness of habitat 2001")


ggplot(land_2016, aes(Red_q0.L, Red_q0)) +
  geom_point()+
  geom_smooth(method = "glm") +
  geom_rug() + 
  ylab("Redundancy of sp community") +
  xlab("Redundancy of landscape") +
  ylim(c(-100,1200)) +
  ggtitle("Redundancy of community as a function of distinctiveness of habitat 2016")




#Distinctiveness
ggplot(land_2001, aes(Dist_q0.L, Dist_q0)) +
  geom_point()+
  geom_smooth(method = "glm") +
  geom_rug() + 
  ylab("Distinctiveness of sp community") +
  xlab("Distinctiveness of landscape") +
  ylim(c(0,0.015)) +
  ggtitle("Distinctiveness of community as a function of distinctiveness of habitat 2001")
  

ggplot(land_2016, aes(Dist_q0.L, Dist_q0)) +
  geom_point()+
  geom_smooth(method = "glm") +
  geom_rug() + 
  ylab("Distinctiveness of sp community") +
  xlab("Distinctiveness of landscape") +
  ylim(c(0,0.015)) +
  ggtitle("Distinctiveness of community as a function of distinctiveness of habitat 2016")


#Representativeness
ggplot(land_2001, aes(Rep_q0.L, Rep_q0)) +
  geom_point()+
  geom_smooth(method = "glm") +
  geom_rug() + 
  ylab("Representativeness of sp community") +
  xlab("Representativeness of landscape") +
  ylim(c(0,1)) +
  ggtitle("Representativeness of community as a function of distinctiveness of habitat 2001")


ggplot(land_2016, aes(Rep_q0.L, Rep_q0)) +
  geom_point()+
  geom_smooth(method = "glm") +
  geom_rug() + 
  ylab("Representativeness of sp community") +
  xlab("Representativeness of landscape") +
  ylim(c(0,1)) +
  ggtitle("Representativeness of community as a function of distinctiveness of habitat 2016")

#Effective number

ggplot(land_2001, aes(Eff_q0.L, Eff_q0.sp)) +
  geom_point()+
  geom_smooth(method = "glm") +
  geom_rug() + 
  ylab("Eff. numbr distinct  sp community") +
  xlab("Eff numbr distinct  landscape") +
  ggtitle("Eff. numbr distinct subcommunities as a function of the Eff. numbr distinct habitat 2001 ")

ggplot(land_sp_2016, aes(Eff_q0.L, Eff_q0.sp)) +
  geom_point()+
  geom_smooth(method = "glm") +
  geom_rug() + 
  ylab("Eff. numbr distinct  sp community") +
  xlab("Eff numbr distinct  landscape") +
  ggtitle("Eff. numbr distinct subcommunities as a function of the Eff. numbr distinct habitat 2016 ", "By segments")




# P L O T S - By segments

#redundancy

ggplot(land_sp_2001, aes(Red_q0.L, Red_q0.sp)) +
  geom_point()+
  geom_smooth(method = "glm") +
  geom_rug() + 
  ylab("Redundancy of sp community") +
  xlab("Redundancy of landscape") +
  ggtitle("Redundancy of community as a function of redundancy of habitat 2001 ", "By segments")

ggplot(land_sp_2016, aes(Red_q0.L, Red_q0.sp)) +
  geom_point()+
  geom_smooth(method = "glm") +
  geom_rug() + 
  ylab("Redundancy of sp community") +
  xlab("Redundancy of landscape") +
  ggtitle("Redundancy of community as a function of redundancy of habitat 2016 ", "By segments")


#distinctiveness

ggplot(land_sp_2001, aes(Dist_q0.L, Dist_q0.sp)) +
  geom_point()+
  geom_smooth(method = "glm") +
  geom_rug() + 
  ylab("distinctiveness of sp community") +
  xlab("distinctiveness of landscape") +
  ggtitle("distinctiveness of community as a function of distinctiveness of habitat 2001 ", "By segments")

ggplot(land_sp_2016, aes(Dist_q0.L, Dist_q0.sp)) +
  geom_point()+
  geom_smooth(method = "glm") +
  geom_rug() + 
  ylab("distinctiveness of sp community") +
  xlab("distinctiveness of landscape") +
  ggtitle("distinctiveness of community as a function of distinctiveness of habitat 2016 ", "By segments")


#Representativeness

ggplot(land_sp_2001, aes(Rep_q0.L, Rep_q0.sp)) +
  geom_point()+
  geom_smooth(method = "glm") +
  geom_rug() + 
  ylab("Representativeness of sp community") +
  xlab("Representativeness of landscape") +
  ggtitle("Representativeness of community as a function of Representativeness of habitat 2001 ", "By segments")

ggplot(land_sp_2016, aes(Rep_q0.L, Rep_q0.sp)) +
  geom_point()+
  geom_smooth(method = "glm") +
  geom_rug() + 
  ylab("Representativeness of sp community") +
  xlab("Representativeness of landscape") +
  ggtitle("Representativeness of community as a function of Representativeness of habitat 2016 ", "By segments")


#Effective number

ggplot(land_sp_2001, aes(Eff_q0.L, Eff_q0.sp)) +
  geom_point()+
  geom_smooth(method = "glm") +
  geom_rug() + 
  ylab("Eff. numbr distinct  sp community") +
  xlab("Eff numbr distinct  landscape") +
  ggtitle("Eff. numbr distinct subcommunities as a function of the Eff. numbr distinct habitat 2001 ", "By segments")

ggplot(land_sp_2016, aes(Eff_q0.L, Eff_q0.sp)) +
  geom_point()+
  geom_smooth(method = "glm") +
  geom_rug() + 
  ylab("Eff. numbr distinct  sp community") +
  xlab("Eff numbr distinct  landscape") +
  ggtitle("Eff. numbr distinct subcommunities as a function of the Eff. numbr distinct habitat 2016 ", "By segments")


#' ##############################################################################################
# 5 - Temperature ----
#
#'###############################################################################################

#' biodiversity of species community by segments

deltas <- data.frame("partition"= rownames(spmatrix.t1.full),
                     "delta.Red_q0"=0, "delta.Dist_q0"=0,
                     "delta.Rep_q0"=0, "delta.Eff_q0"=0)

red_2001 <- raw_sub_rho(meta.sp.2001, 0)
red_2016 <- raw_sub_rho(meta.sp.2016, 0)

dist_2001 <- raw_sub_beta(meta.sp.2001, 0)
dist_2016 <- raw_sub_beta(meta.sp.2016, 0)

rep_2001 <- norm_sub_rho(meta.sp.2001, 0)
rep_2016 <- norm_sub_rho(meta.sp.2016, 0)

eff_2001 <- norm_sub_beta(meta.sp.2001, 0)
eff_2016 <- norm_sub_beta(meta.sp.2016, 0)

deltas$delta.Red_q0 <- red_2016$diversity - red_2001$diversity
deltas$delta.Dist_q0 <- dist_2016$diversity - dist_2001$diversity
deltas$delta.Rep_q0 <- rep_2016$diversity - rep_2001$diversity
deltas$delta.Eff_q0 <- eff_2016$diversity - eff_2001$diversity


temperature <- data.segments.4km %>%  select(partition, delta.temp)

deltas.total <- data.segments.4km %>%  mutate(delta.total=(abs(delta.urban) + abs(delta.forest) +abs(delta.cropland)+abs(delta.grassland) +abs(delta.wetland))/2) %>% 
  select(partition, delta.total)

deltas <- deltas %>%  merge( temperature, by= "partition") %>%  merge(deltas.total, by="partition")

ggplot(deltas, aes(delta.temp, delta.Red_q0)) +
  geom_point() +
  geom_smooth(method="loess")+
  geom_rug()+
  ylab("Redundancy q=0") 

ggplot(deltas, aes(delta.temp, delta.Dist_q0)) +
  geom_point() +
  geom_smooth(method="loess")+
  geom_rug()+
  ylab("Distinctiveness q=0")

ggplot(deltas, aes(delta.temp, delta.Rep_q0)) +
  geom_point() +
  geom_smooth(method="loess")+
  geom_rug()+
  ylab("Representativeness q=0")

ggplot(deltas, aes(delta.temp, delta.Eff_q0)) +
  geom_point() +
  geom_smooth(method="loess")+
  geom_rug()+
  ylab("Effective nr distinct communities q=0")


#Delta total

ggplot(deltas, aes(delta.total, delta.Red_q0)) +
  geom_point() +
  geom_smooth(method="loess")+
  geom_rug()+
  ylab("Redundancy q=0") 

ggplot(deltas, aes(delta.total, delta.Dist_q0)) +
  geom_point() +
  geom_smooth(method="loess")+
  geom_rug()+
  ylab("Distinctiveness q=0")

ggplot(deltas, aes(delta.total, delta.Rep_q0)) +
  geom_point() +
  geom_smooth(method="loess")+
  geom_rug()+
  ylab("Representativeness q=0")

ggplot(deltas, aes(delta.total, delta.Eff_q0)) +
  geom_point() +
  geom_smooth(method="loess")+
  geom_rug()+
  ylab("Effective nr distinct communities q=0")





# ~~~ REDUNDANCY AND DISTINCTIVENESS ~~~  ----


index <- data.frame("RouteName"=transects.names,
                    "Redun_q1.2001"=0,
                    "Redun_q1.2016"=0,
                    "delta.redundancy"=0,
                    "Repre_q1.2001"=0,
                    "Repre_q1.2016"=0,
                    "delta.representativeness"=0)

#redundancy
red_q1.2001 <- raw_sub_rho(meta.2001, 1)
red_q1.2016 <- raw_sub_rho(meta.2016, 1)

index$Redun_q1.2001 <- red_q1.2001$diversity
index$Redun_q1.2016 <- red_q1.2016$diversity
index$delta.redundancy <- index$Redun_q1.2016 - index$Redun_q1.2001

#representativeness
rep_q1.2001 <- norm_sub_rho(meta.2001,1)
rep_q1.2016 <- norm_sub_rho(meta.2016,1)

index$Repre_q1.2001 <- rep_q1.2001$diversity
index$Repre_q1.2016 <- rep_q1.2016$diversity
index$delta.representativeness <- index$Repre_q1.2016 - index$Repre_q1.2001








