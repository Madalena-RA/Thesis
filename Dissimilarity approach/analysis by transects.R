#'
#'*Analysis by transects*
#'
#' 
#'

library(tidyverse)
library(dplyr)
library(betapart)

#working directory
setwd("~/GLASGOW/Dissertation/data")

#abundance data
load("vegan.spmatrix_segments_t1.rda") #data from 2001
load("vegan.spmatrix_segments_t2.rda") #data from 2016
load("data.segments.4km.rda")
routes <- read.csv("C:/Users/madal/Downloads/routes/routes.csv")
states <- read.csv("~/GLASGOW/Dissertation/data/States_Num_Name.csv", sep=";")

routename <- new.dissimilarity %>% select(partition, RouteName)  
transects <- unique(routename$RouteName)

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
spmatrix.t1_binary[spmatrix.t1_binary > 0] = 1

spmatrix.t2_binary <- spmatrix.t2.full
spmatrix.t2_binary[spmatrix.t2_binary > 0] = 1


#BINARY DATA FROM 2001 AND 2016 BY TRANSECTS

data_2001 <- as.data.frame(spmatrix.t1_binary)%>% 
  mutate("partition"=rownames(spmatrix.t1_binary)) %>% 
  merge(routename, by="partition") %>% 
  select(-partition) 

data_2016 <- as.data.frame(spmatrix.t2_binary)%>% 
  mutate("partition"=rownames(spmatrix.t2_binary)) %>% 
  merge(routename, by="partition")  %>% 
  select(-partition) 

new.2001 <- aggregate(. ~ RouteName, data=data_2001, sum) 
rownames(new.2001) <- transects
new.2001 <- new.2001 %>% select(-RouteName)

new.2016<- aggregate(. ~ RouteName, data=data_2016, sum)
rownames(new.2016) <- transects
new.2016 <- new.2016 %>% select(-RouteName)

new.2001[new.2001 > 0] = 1
new.2016[new.2016 > 0] = 1


# NEW HABITAT DATA BY TRANSECTS

data.transects.4km <- data.segments.4km %>% 
  merge(routename, by="partition") %>% 
  select(-partition) 

temperature <- data.transects.4km %>% select(delta.temp, RouteName)
temperature <- aggregate(delta.temp~RouteName, data=temperature, FUN = mean)


data.transects.4km <- aggregate(.~RouteName, data=data.transects.4km, sum)
data.transects.4km <- data.transects.4km %>% select(-delta.temp)
data.transects.4km <- merge(data.transects.4km,temperature, by= "RouteName")





# Jaccard ----

jaccard_by_transects <- data.frame("RouteName" = rownames(new.2001), "Jaccard dissim"=0, "Sorensen dissim"=0, "Jaccard turnover"=0)
jaccard.trans <- beta.temp(x = new.2001, y = new.2016, index.family = "jaccard")
sorensen.trans <- beta.temp(x = new.2001, y = new.2016, index.family = "sorensen")

jaccard_by_transects$Jaccard.dissim <- jaccard.trans$beta.jac     # Jaccard dissimilarity index
jaccard_by_transects$Sorensen.dissim <- sorensen.trans$beta.sor 
jaccard_by_transects$Jaccard.turnover <- jaccard.trans$beta.jtu


dissimilarity_transects <- merge(jaccard_by_transects, data.transects.4km, by="RouteName")




# BCR ----
bcr <- routes %>% select(RouteName, BCR)

dissimilarity_transects <- merge(dissimilarity_transects, bcr, by= "RouteName")

ggplot(dissimilarity_transects, aes(as.factor(BCR), Jaccard.turnover, color=as.factor(BCR))) +
  geom_boxplot()+
  ggtitle("Variation of the Jaccard turnover index by BCR ")

ggplot(dissimilarity_transects, aes(as.factor(BCR), Jaccard.turnover)) +
  geom_boxplot()+
  ggtitle("Variation of the Jaccard turnover index by BCR - transects")

#Delta Urban
ggplot(dissimilarity_transects, aes(delta.urban, Jaccard.turnover)) +
  geom_point()+
  facet_wrap(~BCR) +
  geom_smooth(method="glm") +
  ggtitle("Jaccard turnover as a function of delta urban by BCR - glm transects")

ggplot(dissimilarity_transects, aes(delta.urban, Jaccard.turnover)) +
  geom_point()+
  facet_wrap(~BCR) +
  geom_smooth(method="loess") +
  ggtitle("Jaccard turnover as a function of delta urban by BCR - loess transects")

ggplot(dissimilarity_transects, aes(delta.urban, Jaccard.turnover)) +
  geom_point()+
  geom_smooth(method="loess") +
  ggtitle("Jaccard turnover as a function of delta urban - loess")

#Delta Forest
ggplot(dissimilarity_transects, aes(delta.forest, Jaccard.turnover)) +
  geom_point()+
  facet_wrap(~BCR) +
  geom_smooth(method="glm") +
  ggtitle("Jaccard turnover as a function of delta forest by BCR - glm")

ggplot(dissimilarity_transects, aes(delta.forest, Jaccard.turnover)) +
  geom_point()+
  facet_wrap(~BCR) +
  geom_smooth(method="loess") +
  ggtitle("Jaccard turnover as a function of delta forest by BCR - loess")

ggplot(dissimilarity_transects, aes(delta.forest, Jaccard.turnover)) +
  geom_point()+
  geom_smooth(method="loess") +
  ggtitle("Jaccard turnover as a function of delta forest - loess")


#Delta temperature
ggplot(dissimilarity_transects, aes(delta.temp, Jaccard.turnover)) +
  geom_point()+
  facet_wrap(~BCR) +
  geom_smooth(method="glm") +
  ggtitle("Jaccard turnover as a function of delta temperature by BCR - glm")

ggplot(dissimilarity_transects, aes(delta.temp, Jaccard.turnover)) +
  geom_point()+
  facet_wrap(~BCR) +
  geom_smooth(method="loess") +
  ggtitle("Jaccard turnover as a function of delta temperature by BCR - loess")

ggplot(dissimilarity_transects, aes(delta.temp, Jaccard.turnover)) +
  geom_point()+
  geom_smooth(method="loess") +
  ggtitle("Jaccard turnover as a function of delta temperature - loess")

ggplot(dissimilarity_transects, aes(delta.temp, Jaccard.turnover)) +
  geom_point()+
  geom_smooth(method="glm") +
  ggtitle("Jaccard turnover as a function of delta temperature - loess")


###################################################################################
# Normalised beta diversity ----
#
###################################################################################

meta.2001 <- as.data.frame(spmatrix.t1.full) %>% 
  mutate("partition"=rownames(spmatrix.t1.full)) %>% 
  merge(routename, by="partition") %>% 
  select(-partition) 

meta.2001 <- aggregate(.~RouteName, data = meta.2001, sum)
rownames(meta.2001) <- meta.2001$RouteName
meta.2001 <- select(meta.2001, -RouteName)
meta.2001 <- metacommunity(t(meta.2001))


meta.2016 <- as.data.frame(spmatrix.t2.full) %>% 
  mutate("partition"=rownames(spmatrix.t2.full)) %>% 
  merge(routename, by="partition") %>% 
  select(-partition) 

meta.2016 <- aggregate(.~RouteName, data = meta.2016, sum)
rownames(meta.2016) <- meta.2016$RouteName
meta.2016 <- select(meta.2016, -RouteName)
meta.2016 <- metacommunity(t(meta.2016))


# Measures 2001 ----

measures_2001 <- data.frame("RouteName"=transects,
                            "Red_q0"=0,"Red_q1"=0,"Red_q2"=0,"Red_qinf"=0,
                            "Dist_q0"=0,"Dist_q1"=0,"Dist_q2"=0,"Dist_qinf"=0,
                            "Rep_q0"=0,"Rep_q1"=0,"Rep_q2"=0,"Rep_qinf"=0,
                            "Eff_q0"=0,"Eff_q1"=0,"Eff_q2"=0,"Eff_qinf"=0)



# REDUNDANCY 2001 ----

red_q0.t1 <- raw_sub_rho(meta.2001, 0)
red_q1.t1 <- raw_sub_rho(meta.2001, 1)
red_q2.t1 <- raw_sub_rho(meta.2001, 2)
red_qinf.t1 <- raw_sub_rho(meta.2001, Inf)

measures_2001$Red_q0 <- red_q0.t1$diversity
measures_2001$Red_q1 <- red_q1.t1$diversity
measures_2001$Red_q2 <- red_q2.t1$diversity
measures_2001$Red_qinf <- red_qinf.t1$diversity

#DISTINCTIVENESS 2001 ----

dist_q0.t1 <- raw_sub_beta(meta.2001, 0)
dist_q1.t1 <- raw_sub_beta(meta.2001, 1)
dist_q2.t1 <- raw_sub_beta(meta.2001, 2)
dist_qinf.t1 <- raw_sub_beta(meta.2001, Inf)

measures_2001$Dist_q0 <- dist_q0.t1$diversity
measures_2001$Dist_q1 <- dist_q1.t1$diversity
measures_2001$Dist_q2 <- dist_q2.t1$diversity
measures_2001$Dist_qinf <- dist_qinf.t1$diversity

# REPRESENTATIVENESS 2001 ----

rep_q0.t1 <- norm_sub_rho(meta.2001, 0)
rep_q1.t1 <- norm_sub_rho(meta.2001, 1)
rep_q2.t1 <- norm_sub_rho(meta.2001, 2)
rep_qinf.t1 <- norm_sub_rho(meta.2001, Inf)

measures_2001$Rep_q0 <- rep_q0.t1$diversity
measures_2001$Rep_q1 <- rep_q1.t1$diversity
measures_2001$Rep_q2 <- rep_q2.t1$diversity
measures_2001$Rep_qinf <- rep_qinf.t1$diversity

#EFFECTIVE NUMBER OF DISTINCT SUBCOMMUNITIES  2001 ----

eff_q0.t1 <- norm_sub_beta(meta.2001, 0)
eff_q1.t1 <- norm_sub_beta(meta.2001, 1)
eff_q2.t1 <- norm_sub_beta(meta.2001, 2)
eff_qinf.t1 <- norm_sub_beta(meta.2001, Inf)

measures_2001$Eff_q0 <- eff_q0.t1$diversity
measures_2001$Eff_q1 <- eff_q1.t1$diversity
measures_2001$Eff_q2 <- eff_q2.t1$diversity
measures_2001$Eff_qinf <- eff_qinf.t1$diversity


# Measures 2016 ----

measures_2016 <- data.frame("RouteName"=transects,
                            "Red_q0"=0,"Red_q1"=0,"Red_q2"=0,"Red_qinf"=0,
                            "Dist_q0"=0,"Dist_q1"=0,"Dist_q2"=0,"Dist_qinf"=0,
                            "Rep_q0"=0,"Rep_q1"=0,"Rep_q2"=0,"Rep_qinf"=0,
                            "Eff_q0"=0,"Eff_q1"=0,"Eff_q2"=0,"Eff_qinf"=0)



# REDUNDANCY 2001 ----

red_q0.t2 <- raw_sub_rho(meta.2016, 0)
red_q1.t2 <- raw_sub_rho(meta.2016, 1)
red_q2.t2 <- raw_sub_rho(meta.2016, 2)
red_qinf.t2 <- raw_sub_rho(meta.2016, Inf)

measures_2016$Red_q0 <- red_q0.t2$diversity
measures_2016$Red_q1 <- red_q1.t2$diversity
measures_2016$Red_q2 <- red_q2.t2$diversity
measures_2016$Red_qinf <- red_qinf.t2$diversity

#DISTINCTIVENESS 2001 ----

dist_q0.t1 <- raw_sub_beta(meta.2001, 0)
dist_q1.t1 <- raw_sub_beta(meta.2001, 1)
dist_q2.t1 <- raw_sub_beta(meta.2001, 2)
dist_qinf.t1 <- raw_sub_beta(meta.2001, Inf)

measures_2001$Dist_q0 <- dist_q0.t1$diversity
measures_2001$Dist_q1 <- dist_q1.t1$diversity
measures_2001$Dist_q2 <- dist_q2.t1$diversity
measures_2001$Dist_qinf <- dist_qinf.t1$diversity

# REPRESENTATIVENESS 2001 ----

rep_q0.t1 <- norm_sub_rho(meta.2001, 0)
rep_q1.t1 <- norm_sub_rho(meta.2001, 1)
rep_q2.t1 <- norm_sub_rho(meta.2001, 2)
rep_qinf.t1 <- norm_sub_rho(meta.2001, Inf)

measures_2001$Rep_q0 <- rep_q0.t1$diversity
measures_2001$Rep_q1 <- rep_q1.t1$diversity
measures_2001$Rep_q2 <- rep_q2.t1$diversity
measures_2001$Rep_qinf <- rep_qinf.t1$diversity

#EFFECTIVE NUMBER OF DISTINCT SUBCOMMUNITIES  2001 ----

eff_q0.t1 <- norm_sub_beta(meta.2001, 0)
eff_q1.t1 <- norm_sub_beta(meta.2001, 1)
eff_q2.t1 <- norm_sub_beta(meta.2001, 2)
eff_qinf.t1 <- norm_sub_beta(meta.2001, Inf)

measures_2001$Eff_q0 <- eff_q0.t1$diversity
measures_2001$Eff_q1 <- eff_q1.t1$diversity
measures_2001$Eff_q2 <- eff_q2.t1$diversity
measures_2001$Eff_qinf <- eff_qinf.t1$diversity




t_Nb.q0.t1<-norm_sub_beta(meta.2001, 0)      
t_Nb.q1.t1<-norm_sub_beta(meta.2001, 1)
t_Nb.q2.t1<-norm_sub_beta(meta.2001, 2)
t_Nb.qinf.t1<-norm_sub_beta(meta.2001, Inf)

t_Nb.q0.t2<-norm_sub_beta(meta.2016, 0)      
t_Nb.q1.t2<-norm_sub_beta(meta.2016, 1)
t_Nb.q2.t2<-norm_sub_beta(meta.2016, 2)
t_Nb.qinf.t2<-norm_sub_beta(meta.2016, Inf)


norm.beta.transects$t_Nb.q0.t1 <- t_Nb.q0.t1$diversity
norm.beta.transects$t_Nb.q1.t1 <- t_Nb.q1.t1$diversity
norm.beta.transects$t_Nb.q2.t1 <- t_Nb.q2.t1$diversity
norm.beta.transects$t_Nb.qinf.t1 <- t_Nb.qinf.t1$diversity

norm.beta.transects$t_Nb.q0.t2 <- t_Nb.q0.t2$diversity
norm.beta.transects$t_Nb.q1.t2 <- t_Nb.q1.t2$diversity
norm.beta.transects$t_Nb.q2.t2 <- t_Nb.q2.t2$diversity
norm.beta.transects$t_Nb.qinf.t2 <- t_Nb.qinf.t2$diversity

norm.beta.transects$t_delta.Nb.q0 <- t_Nb.q0.t2$diversity - t_Nb.q0.t1$diversity
norm.beta.transects$t_delta.Nb.q1 <- t_Nb.q1.t2$diversity - t_Nb.q1.t1$diversity
norm.beta.transects$t_delta.Nb.q2 <- t_Nb.q2.t2$diversity - t_Nb.q2.t1$diversity
norm.beta.transects$t_delta.Nb.qinf <- t_Nb.qinf.t2$diversity - t_Nb.qinf.t1$diversity



norm.beta.transects <- merge(norm.beta.transects, dissimilarity_transects, by="RouteName") #contains the BCR and the data 4km

norm.beta.transects<- mutate(norm.beta.transects, delta.total= (abs(delta.urban) + abs(delta.forest) +abs(delta.cropland)+abs(delta.grassland) +abs(delta.wetland))/2)




#######################################################################################
# plots ----
#
#######################################################################################

# FOREST
ggplot(norm.beta.transects, aes(forest.t1,t_Nb.q0.t1)) + 
  geom_point() +
  facet_wrap(~BCR)+
  geom_smooth(method="glm")+
  geom_rug() +
  ggtitle("Distinctiveness q=0 as a function of forest 2001 by BCR")

ggplot(norm.beta.transects, aes(forest.t2,t_Nb.q0.t2)) + 
  geom_point() +
  facet_wrap(~BCR)+
  geom_smooth(method="glm")+
  geom_rug() +
  ggtitle("Distinctiveness q=0 as a function of forest 2016 by BCR")



ggplot(norm.beta.transects, aes(as.factor(BCR),t_Nb.q0.t1)) + 
  geom_boxplot() +
  scale_y_continuous(limits = c(0,15))+
  ggtitle("Variation of normalised beta diversity q=0 in 2001 by transect")


ggplot(norm.beta.transects, aes(as.factor(BCR),t_Nb.q0.t2)) + 
  geom_boxplot() +
  scale_y_continuous(limits = c(0,15))+
  ggtitle("Variation of normalised beta diversity q=0 in 2016 by transects")



ggplot(all, aes(delta.Nb.q0,Jaccard.turnover)) +
  geom_point() +
  geom_smooth(method="glm") +
  geom_rug()+
  ggtitle("Jaccard turnover as a function of delta.q0 by segments")


ggplot(norm.beta.transects, aes(t_delta.Nb.q0,Jaccard.turnover)) +
  geom_point() +
  geom_smooth(method="glm") +
  geom_rug() +
  ggtitle("Jaccard turnover as a function of delta.q0 by transects")


ggplot(norm.beta.transects, aes(delta.total,delta.q0)) +
  geom_point() +
  scale_x_sqrt()+
  geom_smooth(method="loess") +
  geom_rug() +
  ggtitle("Transects")

ggplot(new.dissimilarity, aes(delta.total,delta.q0)) +
  geom_point() +
  scale_x_sqrt()+
  geom_smooth(method="loess") +
  geom_rug()+
  ggtitle("Segments")


ggplot(norm.beta.transects, aes(forest.t2, t_Nb.q0.t2)) +
  geom_point() +
  geom_smooth(method="loess") +
  geom_rug() + 
  ggtitle("Distinctiveness q=0 as a function of forest in 2016 ")

ggplot(norm.beta.transects, aes(forest.t2, t_Nb.q1.t2)) +
  geom_point() +
  geom_smooth(method="loess") +
  geom_rug() + 
  ggtitle("Distinctiveness q=1 as a function of forest in 2016 ")

ggplot(norm.beta.transects, aes(forest.t2, t_Nb.q2.t2)) +
  geom_point() +
  geom_smooth(method="loess") +
  geom_rug() + 
  ggtitle("Distinctiveness q=2 as a function of forest in 2016 ")

ggplot(norm.beta.transects, aes(forest.t2, t_Nb.qinf.t2)) +
  geom_point() +
  geom_smooth(method="loess") +
  geom_rug() + 
  ggtitle("Distinctiveness q=Inf as a function of forest in 2016 ")



