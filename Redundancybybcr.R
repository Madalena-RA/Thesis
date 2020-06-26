#'
#'*Analysis Beta diversities - BY BCR*
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

## 1.2 - 23 km buffers ---- 

lc.01.16.23km <- lc.01.16.23km %>%  rename("RouteName"=partition) %>% merge(bcr, by="RouteName")

#transects 2001 landscape
lc.01.23km <- lc.01.16.23km %>% select(RouteName, water.t1,icesnow.t1,urban.t1,deciduousforest.t1,evergreenforest.t1,
                                       mixedforest.t1,grassland.t1,cropland.t1,shrubland.t1,wetland.t1, delta.urban)
#transects 2016 landscape
lc.16.23km <- lc.01.16.23km %>% select(RouteName, water.t2,icesnow.t2,urban.t2,deciduousforest.t2,evergreenforest.t2,
                                       mixedforest.t2,grassland.t2,cropland.t2,shrubland.t2,wetland.t2, delta.urban)


#'##################################################################################
# Redundancy 2001 ----
#
#'##################################################################################

meta.2001 <- as.data.frame(spmatrix.t1.full) %>% 
  mutate("partition"=rownames(spmatrix.t1.full)) %>%  #creates new column with partitions
  merge(transects, by="partition") %>%                #adds transects names
  select(-partition)                                  #deletes partition column

meta.2001 <- aggregate(.~RouteName, data = meta.2001, sum) #adds observations per RouteName
rownames(meta.2001) <- meta.2001$RouteName                 #defines routenames as rownames
meta.2001 <- merge(meta.2001, bcr, by="RouteName")

meta.2016 <- as.data.frame(spmatrix.t2.full) %>% 
  mutate("partition"=rownames(spmatrix.t2.full)) %>% 
  merge(transects, by="partition") %>% 
  select(-partition) 

meta.2016 <- aggregate(.~RouteName, data = meta.2016, sum)
rownames(meta.2016) <- meta.2016$RouteName
meta.2016 <- merge(meta.2016, bcr, by="RouteName")







bcr.usa <- unique(meta.2001$BCR)

#Redundancy within bcr ---- 
#2001
redundancy.df.01 <- data.frame()

for (r in bcr.usa) {
  data <- meta.2001 %>% filter(BCR==r) %>% column_to_rownames("RouteName") %>% select(-BCR)
  meta.data <- metacommunity(t(data))
  redundancy <- raw_sub_rho(meta.data,1)
  redundancy.clean <- redundancy %>% rename("RouteName"=partition_name) %>% select(RouteName, diversity) %>% 
    merge(bcr,by="RouteName")
  redundancy.df.01 <- rbind(redundancy.df.01, redundancy.clean)
}

redundancy.df.01 <- redundancy.df.01 %>%  rename("Red_bcr.2001"=diversity) %>% select(-BCR)

#2016
redundancy.df.16 <- data.frame()

for (r in bcr.usa) {
  data <- meta.2016 %>% filter(BCR==r) %>% column_to_rownames("RouteName") %>% select(-BCR)
  meta.data <- metacommunity(t(data))
  redundancy <- raw_sub_rho(meta.data,1)
  redundancy.clean <- redundancy %>% rename("RouteName"=partition_name) %>% select(RouteName, diversity) %>% 
    merge(bcr,by="RouteName")
  redundancy.df.16 <- rbind(redundancy.df.16, redundancy.clean)
}

redundancy.df.16 <- redundancy.df.16 %>%  rename("Red_bcr.2016"=diversity) %>% select(-BCR)


#Representativeness within bcr ----
#2001
representativeness.df.01 <- data.frame()

for (r in bcr.usa) {
  data <- meta.2001 %>% filter(BCR==r) %>% column_to_rownames("RouteName") %>% select(-BCR)
  meta.data <- metacommunity(t(data))
  representativeness<- norm_sub_rho(meta.data,1)
  representativeness.clean <- representativeness %>% rename("RouteName"=partition_name) %>% select(RouteName, diversity) %>% 
    merge(bcr,by="RouteName")
  representativeness.df.01 <- rbind(representativeness.df.01, representativeness.clean)
}
representativeness.df.01 <- representativeness.df.01 %>% rename("Rep_bcr.2001"=diversity) %>% select(-BCR)

#2016
representativeness.df.16 <- data.frame()

for (r in bcr.usa) {
  data <- meta.2016 %>% filter(BCR==r) %>% column_to_rownames("RouteName") %>% select(-BCR)
  meta.data <- metacommunity(t(data))
  representativeness<- norm_sub_rho(meta.data,1)
  representativeness.clean <- representativeness %>% rename("RouteName"=partition_name) %>% select(RouteName, diversity) %>% 
    merge(bcr,by="RouteName")
  representativeness.df.16 <- rbind(representativeness.df.16, representativeness.clean)
}

representativeness.df.16 <- representativeness.df.16 %>% rename("Rep_bcr.2016"=diversity) %>% select(-BCR)

index.df <- data.frame("RouteName"=transects.names)

index.df <- index.df %>% merge(redundancy.df.01, by="RouteName") %>% 
  merge(redundancy.df.16, by="RouteName") %>% 
  merge(representativeness.df.01, by="RouteName") %>% 
  merge(representativeness.df.16, by="RouteName") %>% 
  merge(bcr, by="RouteName")


#Redundancy within the metacommunity ----- 
meta.01 <- metacommunity(t(select(meta.2001, -RouteName, -BCR)))
meta.16 <- metacommunity(t(select(meta.2016, -RouteName, -BCR)))

red_q1.2001 <- raw_sub_rho(meta.01, 1)
red_q1.2016 <- raw_sub_rho(meta.16, 1)

index.df <- index.df %>% mutate("Red_usa.2001"=red_q1.2001$diversity) %>% 
  mutate("Red_usa.2016"=red_q1.2016$diversity)


#Representativeness within the metacommunity -----
rep_q1.2001 <- norm_sub_rho(meta.01, 1)
rep_q1.2016 <- norm_sub_rho(meta.16, 1)

index.df <- index.df %>% mutate("Rep_usa.2001"=rep_q1.2001$diversity) %>% 
  mutate("Rep_usa.2016"=rep_q1.2016$diversity)




index.lc.01 <- index.df %>% select(-Red_bcr.2016,-Rep_bcr.2016,-Red_usa.2016,-Rep_usa.2016) %>% merge(lc.01.23km, by="RouteName")
index.lc.16 <- index.df %>% select(-Red_bcr.2001,-Rep_bcr.2001,-Red_usa.2001,-Rep_usa.2001) %>% merge(lc.16.23km, by="RouteName")

# PLOTS -----

#comparison between BCR and the whole metacommunity

#redundancy 2001
ggplot(index.lc.01, aes(as.factor(BCR), Red_bcr.2001)) +
  geom_boxplot()+
  ggtitle("Redundancy q=1 by BCR 2001","Redundancy within regions")+
  ylab("Redundancy")+
  xlab("BCR") +
  ylim(c(0,150))

ggplot(index.lc.01, aes(as.factor(BCR), Red_usa.2001)) +
  geom_boxplot()+
  ggtitle("Redundancy q=1 by BCR 2001","Redundancy within the whole metacommunity")+
  ylab("Redundancy")+
  xlab("BCR") +
  ylim(c(0,850))

#representativeness 2001
ggplot(index.lc.01, aes(as.factor(BCR), Rep_bcr.2001)) +
  geom_boxplot()+
  ggtitle("Representativeness q=1 by BCR 2001","Representativeness within regions")+
  ylab("Representativeness")+
  xlab("BCR") +
  ylim(c(0,1))

ggplot(index.lc.01, aes(as.factor(BCR), Rep_usa.2001)) +
  geom_boxplot()+
  ggtitle("Representativeness q=1 by BCR 2001","Representativeness within the whole metacommunity")+
  ylab("Representativeness")+
  xlab("BCR") +
  ylim(c(0,1))

#redundancy 2016
ggplot(index.lc.16, aes(as.factor(BCR), Red_bcr.2016)) +
  geom_boxplot()+
  ggtitle("Redundancy q=1 by BCR 2016","Redundancy within regions")+
  ylab("Redundancy")+
  xlab("BCR") + 
  ylim(c(0,150))

ggplot(index.lc.16, aes(as.factor(BCR), Red_usa.2016)) +
  geom_boxplot()+
  ggtitle("Redundancy q=1 by BCR 2016","Redundancy within the whole metacommunity")+
  ylab("Redundancy")+
  xlab("BCR") +
  ylim(c(0,850))

#representativeness 2016
ggplot(index.lc.16, aes(as.factor(BCR), Rep_bcr.2016)) +
  geom_boxplot()+
  ggtitle("Representativeness q=1 by BCR 2016","Representativeness within regions")+
  ylab("Representativeness")+
  xlab("BCR") +
  ylim(c(0,1))

ggplot(index.lc.16, aes(as.factor(BCR), Rep_usa.2016)) +
  geom_boxplot()+
  ggtitle("Representativeness q=1 by BCR 2016","Representativeness within the whole metacommunity")+
  ylab("Representativeness")+
  xlab("BCR") +
  ylim(c(0,1))


#Redundancy


ggplot(index.lc.01, aes(urban.t1, Red_bcr.2001, col=as.factor(BCR))) +
  geom_point() +
  geom_smooth(method="loess", color="black") +
  geom_rug() + 
  ggtitle("Redundancy q=1 as a function of urban in 2001", "Redundancy within regions") +
  ylab("Redundancy") +
  ylim(c(0,150))



