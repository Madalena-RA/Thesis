#packages
library(tidyverse)
library(vegan)
library(betapart)
library(ggplot2)
library(rdiversity)

#working directory
setwd("~/GLASGOW/Dissertation/data")

source("USBBS_Diversity_Scripts/fun_berger_parker.R")

#'##########################################################################################
#'
#' *Import sp matrices and setup diversity datafile*
#'
#
#abundance data
load("vegan.spmatrix_segments_t1.rda") #data from 2001
load("vegan.spmatrix_segments_t2.rda") #data from 2016
load("data.segments.4km.rda")
routes <- read.csv("C:/Users/madal/Downloads/routes/routes.csv")
states <- read.csv("~/GLASGOW/Dissertation/data/States_Num_Name.csv", sep=";")


#they differ in the number of species (columns)
dim(spmatrix.t1)
dim(spmatrix.t2)


#'#########################################################################################
#' 
#' *Distinctiveness approach*
#'
#'

#' 
#' *Raw subcommunity beta diversity for timepoint 1 and 2*
#
# Raw beta ----


meta.t1<-metacommunity(t(spmatrix.t1)) #transform into metacommunity object
meta.t2<-metacommunity(t(spmatrix.t2)) #transform into metacommunity object

raw.sub.beta <- data.frame("partition" = rownames(spmatrix.t1),
                           "b.q0.t1"=0, "b.q1.t1"=0, "b.q2.t1"=0, "b.qinf.t1"=0,
                           "b.q0.t2"=0, "b.q1.t2"=0, "b.q2.t2"=0, "b.qinf.t2"=0, 
                           "delta.b.q0"=0,"delta.b.q1"=0,"delta.b.q2"=0,"delta.b.qinf"=0)

q0.t1<-raw_sub_beta(meta.t1, 0)      
q1.t1<-raw_sub_beta(meta.t1, 1)
q2.t1<-raw_sub_beta(meta.t1, 2)
qinf.t1<-raw_sub_beta(meta.t1, Inf)

q0.t2<-raw_sub_beta(meta.t2, 0)      
q1.t2<-raw_sub_beta(meta.t2, 1)
q2.t2<-raw_sub_beta(meta.t2, 2)
qinf.t2<-raw_sub_beta(meta.t2, Inf)


raw.sub.beta$b.q0.t1 <- q0.t1$diversity
raw.sub.beta$b.q1.t1 <- q1.t1$diversity
raw.sub.beta$b.q2.t1 <- q2.t1$diversity
raw.sub.beta$b.qinf.t1 <- qinf.t1$diversity

raw.sub.beta$b.q0.t2 <- q0.t2$diversity
raw.sub.beta$b.q1.t2 <- q1.t2$diversity
raw.sub.beta$b.q2.t2 <- q2.t2$diversity
raw.sub.beta$b.qinf.t2 <- qinf.t2$diversity

raw.sub.beta$delta.b.q0 <- q0.t2$diversity - q0.t1$diversity
raw.sub.beta$delta.b.q1 <- q1.t2$diversity - q1.t1$diversity
raw.sub.beta$delta.b.q2 <- q2.t2$diversity - q2.t1$diversity
raw.sub.beta$delta.b.qinf <- qinf.t2$diversity - qinf.t1$diversity


raw.sub.beta <- merge(raw.sub.beta, data.segments.4km, by="partition")

#' **Raw subcommunity beta diversity plots**
#' **2001 raw beta plots**

ggplot(raw.sub.beta, aes(urban.t1,b.q0.t1)) + 
  geom_point() +
  geom_smooth(method="loess")+
  geom_rug() +
  scale_x_sqrt()

ggplot(raw.sub.beta, aes(urban.t1,b.q1.t1)) + 
  geom_point() +
  geom_smooth(method="loess") +
  geom_rug() +
  scale_x_sqrt()

ggplot(raw.sub.beta, aes(urban.t1,b.q2.t1)) + 
  geom_point() +
  geom_smooth(method="loess") +
  geom_rug() +
  scale_x_sqrt()

ggplot(raw.sub.beta, aes(urban.t1,b.qinf.t1)) + 
  geom_point() +
  geom_smooth(method="loess") +
  geom_rug() +
  scale_x_sqrt()

#' **2016 raw beta plots**

ggplot(raw.sub.beta, aes(urban.t2,b.q0.t2)) + 
  geom_point() +
  geom_smooth(method="loess")+
  geom_rug() +
  scale_x_sqrt()

ggplot(raw.sub.beta, aes(urban.t2,b.q1.t2)) + 
  geom_point() +
  geom_smooth(method="loess") +
  geom_rug() +
  scale_x_sqrt()

ggplot(raw.sub.beta, aes(urban.t2,b.q2.t2)) + 
  geom_point() +
  geom_smooth(method="loess") +
  geom_rug() +
  scale_x_sqrt()

ggplot(raw.sub.beta, aes(urban.t2,b.qinf.t2)) + 
  geom_point() +
  geom_smooth(method="loess") +
  geom_rug() +
  scale_x_sqrt()

#' **Delta raw beta plots**
ggplot(raw.sub.beta, aes(forest.t2,b.q0.t2)) + 
  geom_point() +
  geom_smooth(method="loess")+
  geom_rug() 

ggplot(raw.sub.beta, aes(delta.urban,delta.b.q1)) + 
  geom_point() +
  geom_smooth(method="loess") +
  geom_rug() +
  scale_x_sqrt()

ggplot(raw.sub.beta, aes(delta.urban,delta.b.q2)) + 
  geom_point() +
  geom_smooth(method="loess") +
  geom_rug() +
  scale_x_sqrt()

ggplot(raw.sub.beta, aes(delta.urban,delta.b.qinf)) + 
  geom_point() +
  geom_smooth(method="loess") +
  geom_rug() +
  scale_x_sqrt()



#' 
#' *Normalised subcommunity beta diversity for timepoint 1 and 2*
#
# Normalised beta ----


norm.sub.beta <- data.frame("partition" = rownames(spmatrix.t1),
                           "Nb.q0.t1"=0, "Nb.q1.t1"=0, "Nb.q2.t1"=0, "Nb.qinf.t1"=0,
                           "Nb.q0.t2"=0, "Nb.q1.t2"=0, "Nb.q2.t2"=0, "Nb.qinf.t2"=0, 
                           "delta.Nb.q0"=0,"delta.Nb.q1"=0,"delta.Nb.q2"=0,"delta.Nb.qinf"=0)

Nb.q0.t1<-norm_sub_beta(meta.t1, 0)      
Nb.q1.t1<-norm_sub_beta(meta.t1, 1)
Nb.q2.t1<-norm_sub_beta(meta.t1, 2)
Nb.qinf.t1<-norm_sub_beta(meta.t1, Inf)

Nb.q0.t2<-norm_sub_beta(meta.t2, 0)      
Nb.q1.t2<-norm_sub_beta(meta.t2, 1)
Nb.q2.t2<-norm_sub_beta(meta.t2, 2)
Nb.qinf.t2<-norm_sub_beta(meta.t2, Inf)


norm.sub.beta$Nb.q0.t1 <- Nb.q0.t1$diversity
norm.sub.beta$Nb.q1.t1 <- Nb.q1.t1$diversity
norm.sub.beta$Nb.q2.t1 <- Nb.q2.t1$diversity
norm.sub.beta$Nb.qinf.t1 <- Nb.qinf.t1$diversity

norm.sub.beta$Nb.q0.t2 <- Nb.q0.t2$diversity
norm.sub.beta$Nb.q1.t2 <- Nb.q1.t2$diversity
norm.sub.beta$Nb.q2.t2 <- Nb.q2.t2$diversity
norm.sub.beta$Nb.qinf.t2 <- Nb.qinf.t2$diversity

norm.sub.beta$delta.Nb.q0 <- Nb.q0.t2$diversity - Nb.q0.t1$diversity
norm.sub.beta$delta.Nb.q1 <- Nb.q1.t2$diversity - Nb.q1.t1$diversity
norm.sub.beta$delta.Nb.q2 <- Nb.q2.t2$diversity - Nb.q2.t1$diversity
norm.sub.beta$delta.Nb.qinf <- Nb.qinf.t2$diversity - Nb.qinf.t1$diversity


norm.sub.beta <- merge(norm.sub.beta, data.segments.4km, by="partition")
norm.sub.beta <-  merge(norm.sub.beta,routename, by="partition")
norm.sub.beta <- merge(norm.sub.beta, bcr, by= "RouteName")

#' **Normalised subcommunity beta diversity plots**
#' **2001 normalised beta plots**

ggplot(norm.sub.beta, aes(urban.t1,Nb.q0.t1)) + 
  geom_point() +
  geom_smooth(method="loess")+
  geom_rug() 

ggplot(norm.sub.beta, aes(urban.t1,Nb.q1.t1)) + 
  geom_point() +
  geom_smooth(method="loess") +
  geom_rug() 

ggplot(norm.sub.beta, aes(urban.t1,Nb.q2.t1)) + 
  geom_point() +
  geom_smooth(method="loess") +
  geom_rug() 

ggplot(norm.sub.beta, aes(urban.t1,Nb.qinf.t1)) + 
  geom_point() +
  geom_smooth(method="loess") +
  geom_rug() 

#' **2016 Normalised beta plots**

ggplot(norm.sub.beta, aes(urban.t2,Nb.q0.t2)) + 
  geom_point() +
  geom_smooth(method="loess")+
  geom_rug() 

ggplot(norm.sub.beta, aes(urban.t2,Nb.q1.t2)) + 
  geom_point() +
  geom_smooth(method="loess") +
  geom_rug() 

ggplot(norm.sub.beta, aes(urban.t2,Nb.q2.t2)) + 
  geom_point() +
  geom_smooth(method="loess") +
  geom_rug() 

ggplot(norm.sub.beta, aes(urban.t2,Nb.qinf.t2)) + 
  geom_point() +
  geom_smooth(method="loess") +
  geom_rug()

#' **Delta raw beta plots**
ggplot(norm.sub.beta, aes(delta.urban,delta.Nb.q0)) + 
  geom_point() +
  geom_smooth(method="loess")+
  geom_rug() 

ggplot(norm.sub.beta, aes(delta.urban,delta.Nb.q1)) + 
  geom_point() +
  geom_smooth(method="loess") +
  geom_rug() 

ggplot(norm.sub.beta, aes(delta.urban,delta.Nb.q2)) + 
  geom_point() +
  geom_smooth(method="loess") +
  geom_rug() 

ggplot(norm.sub.beta, aes(delta.urban,delta.Nb.qinf)) + 
  geom_point() +
  geom_smooth(method="loess") +
  geom_rug() 


#' **Variation plots**
#' 
ggplot(norm.sub.beta, aes(as.factor(BCR),Nb.q0.t1)) + 
  geom_boxplot() +
  scale_y_continuous(limits = c(0,60))+
  ggtitle("Variation of normalised beta diversity q=0 in 2001 by segments")

ggplot(norm.sub.beta, aes(as.factor(BCR),Nb.q0.t2)) + 
  geom_boxplot() +
  scale_y_continuous(limits = c(0,60))+
  ggtitle("Variation of normalised beta diversity q=0 in 2016 by segments")

summary(norm.sub.beta$Nb.q0.t1)


#'#########################################################################################
#' 
#' **Dissimilarity approach** 
#'
#'

  # Sorensen dissimilarity index

#dissimilarity$Jac.turnover <- jaccard$beta.jtu       # Turnover component of Jaccard 
#dissimilarity$Simpson.dissim <- sorensen$beta.sim    # Turnover component of Sorensen (Simpson dissimilarity index)
#dissimilarity$Jac.nestedness <- jaccard$beta.jne     # Nestedness component of Jaccard
#dissimilarity$Sor.nestedness <- sorensen$beta.sne    # Nestedness component of Sorensen

#'
#'
#' *Jaccard = b+c/a+b+c*dim(spmatrix.t1)
dim(spmatrix.t2)

#matrix 1 has less 10 species than matrix 2

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

#same size
dim(spmatrix.t1.full)
dim(spmatrix.t2.full)

n1 <- colnames(spmatrix.t1.full)
n2 <- colnames(spmatrix.t2.full)

length(setdiff(n1,n2)) #do not differ in species 


#'##########################################################################################
#' 
#' **Jaccard and Sorensen Indices**
#'
# Jaccard and Sorensen by segments ----
#' 
#turning both matrices into binary data
#presence/absence data
spmatrix.t1_binary <- spmatrix.t1.full
spmatrix.t1_binary[spmatrix.t1_binary > 0] = 1

spmatrix.t2_binary <- spmatrix.t2.full
spmatrix.t2_binary[spmatrix.t2_binary > 0] = 1


#dissimilarity <- data.frame("partition" = rownames(spmatrix.t1), "Jaccard dissim"=0, "Sorensen dissim"=0,  "Jac turnover"=0, "Simpson dissim"=0, "Jac nestedness"=0, "Sor nestedness"=0)

dissimilarity <- data.frame("partition" = rownames(spmatrix.t1), "Jaccard dissim"=0, "Sorensen dissim"=0, "Jaccard turnover"=0)
jaccard <- beta.temp(x = spmatrix.t1_binary, y = spmatrix.t2_binary, index.family = "jaccard")
sorensen <- beta.temp(x = spmatrix.t1_binary, y = spmatrix.t2_binary, index.family = "sorensen")

dissimilarity$Jaccard.dissim <- jaccard$beta.jac     # Jaccard dissimilarity index
dissimilarity$Sorensen.dissim <- sorensen$beta.sor 
dissimilarity$Jaccard.turnover <- jaccard$beta.jtu



#Check how the Jaccard and Sorensen are calculated ----

#' Jaccard = b+c/a+b+c
#' Sorensen = b+c/2a+b+c



#' ########################################################################################
#' 
#' 
#' 
#' Merge the dissimilarity dataframe with data.segments
#' 
#' 
new.dissimilarity <- merge(dissimilarity, data.segments.4km, by="partition")

new.dissimilarity <-  mutate(new.dissimilarity, RouteName= substr(new.dissimilarity$partition,1,nchar(as.vector(new.dissimilarity$partition))-2))

length(unique(new.dissimilarity$RouteName))

#left_join

routes <- select(routes,StateNum, RouteName)

new.dissimilarity<- merge(new.dissimilarity, routes, by="RouteName")

usa <- states %>% filter(CountryNum=="840") %>% rename(state=State.Prov.TerrName)

new.dissimilarity<- merge(new.dissimilarity, usa, by="StateNum")

length(unique(new.dissimilarity$state)) #47?? missing Rhode Island smallest state
new.dissimilarity<- mutate(new.dissimilarity, delta.total= (abs(delta.urban) + abs(delta.forest) +abs(delta.cropland)+abs(delta.grassland) +abs(delta.wetland))/2)


all <- merge(new.dissimilarity, norm.sub.beta, by ="partition")


#' Plots - Jaccard index
ggplot(all, aes(delta.Nb.q0,Jaccard.turnover)) +
  geom_point() +
  geom_smooth(method="loess") +
  geom_rug() 

ggplot(new.dissimilarity, aes(delta.total,delta.q0)) +
  geom_point() +
  scale_x_sqrt()+
  geom_smooth(method="loess") +
  geom_rug()

ggplot(new.dissimilarity, aes(delta.forest,Jaccard.dissim )) +
  geom_point() +
  scale_x_sqrt()+
  geom_smooth(method="loess") 



ggplot(new.dissimilarity, aes(urban.t2,Jaccard.dissim )) +
  geom_point() +
  geom_smooth(method="loess")+
  scale_x_sqrt()

ggplot(new.dissimilarity, aes(delta.forest,Jaccard.dissim)) +
  geom_point() + 
  geom_smooth(method="loess")

ggplot(new.dissimilarity, aes(delta.grassland,Jaccard.dissim)) +
  geom_point()+ 
  geom_smooth(method="loess")

ggplot(new.dissimilarity, aes(delta.cropland,Jaccard.dissim)) +
  geom_point() + 
  geom_smooth(method="loess")

ggplot(new.dissimilarity, aes(delta.wetland,Jaccard.dissim)) +
  geom_point() + 
  geom_smooth(method="loess") 

ggplot(new.dissimilarity, aes(delta.temp,Jaccard.dissim)) +
  geom_point() +
  geom_smooth(method="loess")


#' Plots - Sorenses index

ggplot(new.dissimilarity, aes(delta.urban,Sorensen.dissim )) +
  geom_point() +
  geom_smooth(method="loess")+
  scale_x_sqrt()

ggplot(new.dissimilarity, aes(urban.t2,Sorensen.dissim )) +
  geom_point() +
  geom_smooth(method="loess")+
  scale_x_sqrt()

ggplot(new.dissimilarity, aes(delta.forest,Sorensen.dissim)) +
  geom_point() + 
  geom_smooth(method="loess")

ggplot(new.dissimilarity, aes(delta.grassland,Sorensen.dissim)) +
  geom_point()+ 
  geom_smooth(method="loess")

ggplot(new.dissimilarity, aes(delta.cropland,Sorensen.dissim)) +
  geom_point() + 
  geom_smooth(method="loess")

ggplot(new.dissimilarity, aes(delta.wetland,Sorensen.dissim)) +
  geom_point() + 
  geom_smooth(method="loess") 

ggplot(new.dissimilarity, aes(delta.temp,Sorensen.dissim, color=as.factor(StateNum))) +
  geom_point() +
  geom_smooth(method="loess")


#'##########################################################################################
#' 
#' **Bray-Curtis Index**
#' 
# Bray Curtis ----
#' 

rownames(spmatrix.t2.full) <- paste(rownames(spmatrix.t2.full), "t2", sep = ".")
spmatrix.t1t2 <- cbind(spmatrix.t1.full, spmatrix.t2.full)
dim(spmatrix.t1t2)

bray <- vegdist(spmatrix.t1t2, "bray")
distmat <- as.matrix(bray) 








































#'########################################################################################
#'
#' ** Analyse the three "outliers" **
#' 
# "Outliers" ----
#' 

high_dissimilarity <- subset(new.dissimilarity[new.dissimilarity$Jaccard.dissim>0.8,])
#Featherville 5
#Fox Park 5
#Warner Vally 1

write.csv(high_dissimilarity,"~\\GLASGOW\\Dissertation\\data\\high_dissimilarity.analysed.csv", row.names = TRUE) #save file


high_dissimilarity_t1 <- t(data.frame("FEATHERVILLE 5"=spmatrix.t1.full["FEATHERVILLE 5",] ,
                                      "FOX PARK 5" = spmatrix.t1.full["FOX PARK 5",],
                                      "WARNER VALLY 1" = spmatrix.t1.full["WARNER VALLY 1",]))

write.csv(high_dissimilarity_t1,"~\\GLASGOW\\Dissertation\\data\\high_dissimilarity.t1.csv", row.names = TRUE) #save file

high_dissimilarity_t2 <- t(data.frame("FEATHERVILLE 5"=spmatrix.t2.full["FEATHERVILLE 5",] ,
                                      "FOX PARK 5" = spmatrix.t2.full["FOX PARK 5",],
                                      "WARNER VALLY 1" = spmatrix.t2.full["WARNER VALLY 1",]))

write.csv(high_dissimilarity_t2,"~\\GLASGOW\\Dissertation\\data\\high_dissimilarity.t2.csv", row.names = TRUE) #save file


high_dissimilarity_delta <- t(data.frame("FEATHERVILLE 5"=spmatrix.t2.full["FEATHERVILLE 5",] -spmatrix.t1.full["FEATHERVILLE 5",] ,
                                         "FOX PARK 5" = spmatrix.t2.full["FOX PARK 5",] - spmatrix.t1.full["FOX PARK 5",],
                                         "WARNER VALLY 1" = spmatrix.t2.full["WARNER VALLY 1",] - spmatrix.t1.full["WARNER VALLY 1",]))

#lost species ---- 
species <- colnames(high_dissimilarity_delta)

featherville.lost.species<- as.data.frame(t(high_dissimilarity_delta)) %>%  mutate(species = species) %>% filter(FEATHERVILLE.5 <0) %>% select(FEATHERVILLE.5, species)
foxpark.lost.species<- as.data.frame(t(high_dissimilarity_delta)) %>%  mutate(species = species) %>% filter(FOX.PARK.5 <0) %>% select(FOX.PARK.5, species)
warner.lost.species<- as.data.frame(t(high_dissimilarity_delta)) %>%  mutate(species = species) %>% filter(WARNER.VALLY.1 <0) %>% select(WARNER.VALLY.1, species)

 
# Distinctiveness by habitat cover ----

#habitat matrices

habitat.t1 <- data.segments.4km %>% select( partition, water.t1, urban.t1, forest.t1, grassland.t1, cropland.t1,wetland.t1) 

habitat.t2 <- data.segments.4km %>% select(partition, water.t2, urban.t2, forest.t2, grassland.t2, cropland.t2,wetland.t2) %>% 
  rename(water=water.t2, urban=urban.t2, forest=forest.t2, grassland=grassland.t2, cropland=cropland.t2, wetland=wetland.t2)


meta.habitat.t1 <- metacommunity(t(habitat.t1[,c(2:7)]))
meta.habitat.t2 <- metacommunity(t(habitat.t2[,c(2:7)]))



habitat.norm.beta <- data.frame("partition" = habitat.t1$partition,
                            "habitat.q0.t1"=0, "habitat.q0.t2"=0, "habitat.q0.delta"=0, "Jacctur.habitat"=0)

habitat.q0.t1<-norm_sub_beta(meta.habitat.t1, 0)      
habitat.q0.t2<-norm_sub_beta(meta.habitat.t2, 0)    

habitat.norm.beta$habitat.q0.t1 <- habitat.q0.t1$diversity
habitat.norm.beta$habitat.q0.t2 <- habitat.q0.t2$diversity
habitat.norm.beta$habitat.q0.delta <- habitat.q0.t2$diversity -  habitat.q0.t1$diversity

jaccard.hab <- beta.temp(x = meta.habitat.t1, y = meta.habitat.t2, index.family = "jaccard")
habitat.norm.beta$Jacctur.habitat <- jaccard.hab$beta.jtu

habitat.bird.diss <- merge(habitat.norm.beta, norm.sub.beta, by="partition")


ggplot(habitat.bird.diss, aes(habitat.q0.delta, delta.Nb.q0)) +
  geom_point()



# BCR ----
#' 
routename <- new.dissimilarity %>% select(partition, RouteName)  
bcr <- routes %>% select(RouteName, BCR)

biomes <- biomes %>%  
  merge(dissimilarity, data.segments.4km,"partition") 
  
biomes <- merge(biomes,routename,by="partition")
biomes <- merge(biomes, bcr, by= "RouteName")
sort(unique(biomes$BCR))

#BCRs Jaccard turnover
ggplot(biomes, aes(as.factor(BCR), Jaccard.turnover, color=as.factor(BCR))) +
  geom_boxplot()+
  ggtitle("Variation of the Jaccard turnover index by BCR ")

ggplot(biomes, aes(as.factor(BCR), Jaccard.turnover)) +
  geom_boxplot()+
  ggtitle("Variation of the Jaccard turnover index by BCR")

#Delta Urban
ggplot(biomes, aes(delta.urban, Jaccard.turnover)) +
  geom_point()+
  facet_wrap(~BCR) +
  geom_smooth(method="glm") +
  ggtitle("Jaccard turnover as a function of delta urban by BCR - glm")

ggplot(biomes, aes(delta.urban, Jaccard.turnover)) +
  geom_point()+
  facet_wrap(~BCR) +
  geom_smooth(method="loess") +
  ggtitle("Jaccard turnover as a function of delta urban by BCR - loess")

ggplot(biomes, aes(delta.urban, Jaccard.turnover)) +
  geom_point()+
  geom_smooth(method="loess") +
  ggtitle("Jaccard turnover as a function of delta urban - loess")

#Delta Forest
ggplot(biomes, aes(delta.forest, Jaccard.turnover)) +
  geom_point()+
  facet_wrap(~BCR) +
  geom_smooth(method="glm") +
  ggtitle("Jaccard turnover as a function of delta forest by BCR - glm")

ggplot(biomes, aes(delta.forest, Jaccard.turnover)) +
  geom_point()+
  facet_wrap(~BCR) +
  geom_smooth(method="loess") +
  ggtitle("Jaccard turnover as a function of delta forest by BCR - loess")

ggplot(biomes, aes(delta.forest, Jaccard.turnover)) +
  geom_point()+
  geom_smooth(method="loess") +
  ggtitle("Jaccard turnover as a function of delta forest - loess")


#Delta temperature
ggplot(biomes, aes(delta.temp, Jaccard.turnover)) +
  geom_point()+
  facet_wrap(~BCR) +
  geom_smooth(method="glm") +
  ggtitle("Jaccard turnover as a function of delta temperature by BCR - glm")

ggplot(biomes, aes(delta.temp, Jaccard.turnover)) +
  geom_point()+
  facet_wrap(~BCR) +
  geom_smooth(method="loess") +
  ggtitle("Jaccard turnover as a function of delta temperature by BCR - loess")

ggplot(biomes, aes(delta.temp, Jaccard.turnover)) +
  geom_point()+
  geom_smooth(method="loess") +
  ggtitle("Jaccard turnover as a function of delta temperature - loess")


# Analysis by transects ----








