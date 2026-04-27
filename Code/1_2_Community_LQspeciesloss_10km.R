#
#  This R script is used to extract number of extinct species since late Quaternary for 410 
#     Neotropical, Afrotropical, and Indomalayan mammal communities from Rowan et al. (2020) PNAS
# 
#
#  Extract species occurrence from range maps of PHYLACINE dataset v1.2.1  (Faurby et al. 2018 Ecology) 
#      https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2443#support-information-section
#     data repository: https://megapast2future.github.io/PHYLACINE_1.2/
#
#  Output 
#    PHYLACINE_LQ_sp_gain_loss_no_ocean_no_homo_410community.RData
#


library(terra); library(sf)
library(dplyr)
setwd(".")

Buff.10km.cov <- read.csv("Data/Rowan_Comm_10kmBuff_Env_cov.csv")

Comm.pt <- st_as_sf(Buff.10km.cov, coords = c("Long", "Lat"), crs = st_crs("EPSG:4326")) 
Comm.10km <- st_buffer(Comm.pt, dist = 10000)

#I. Collect species occurrence based on range map overlapping 10-km buffer of community ========== 
## Use Mass values of PHYLACINE 1.2 to collect species >= 500 g  
### Range maps, trait values, and phylogeny all obtained from PHYLACINE v1.2.1 (https://megapast2future.github.io/PHYLACINE_1.2/)
mamm.func <- read.csv("./Trait_data.csv")
mamm.func <- mamm.func[!duplicated(mamm.func$Binomial.1.2),]
#Exclude mammals < 500g, specialized for oceanic and freshwater: Family of Otariidae, Phocidae, Odobenidae, and genus of Homo
terrestrial.large.sp.exclude.marine <- mamm.func[mamm.func$Terrestrial == 1 & mamm.func$Mass.g >= 500 & 
                                                   !mamm.func$Family.1.2 %in% c("Otariidae", "Phocidae", "Odobenidae") &
                                                   !mamm.func$Genus.1.2 == "Homo",]$Binomial.1.2
length(terrestrial.large.sp.exclude.marine) #1578 sp  

##1. stack present-natural ranges ------
nature.rast.dir <- "."
nature.present.file  <- list.files(nature.rast.dir)
nature.present.file.terr <- nature.present.file[nature.present.file %in% paste0(terrestrial.large.sp.exclude.marine, ".tif")]
length(nature.present.file.terr) #1578 sp without marine 

nature.present <- rast(list(rast(paste0(nature.rast.dir, nature.present.file.terr))))
Comm.10km.phyla <- st_transform(Comm.10km, st_crs(nature.present))

nature.comm.matrix <- data.frame()
for (i in 1:length(Buff.10km.cov$Community)) {
  print(i)
  nature.comm.matrix <- rbind(nature.comm.matrix, 
                              data.frame(Community = Comm.10km.phyla[i,]$Community_unify,
                                         extract(nature.present, Comm.10km.phyla[i,], method = "simple")))  
}

nature.comm.matrix.singe <- nature.comm.matrix[!duplicated(nature.comm.matrix$Community),]
nature.comm.matrix.t <- data.frame(t(nature.comm.matrix.singe[,-c(1:2)]))
colnames(nature.comm.matrix.t) <- nature.comm.matrix.singe$Community

nature.present.comm.SpList <- data.frame()
for ( j in 1:ncol(nature.comm.matrix.t)) {
    comm.occur <- data.frame(nature.comm.matrix.t[,j],
                             Species =  rownames(nature.comm.matrix.t))
    colnames(comm.occur) <- c("Occur", "Species")
    nature.present.comm.SpList <- rbind(nature.present.comm.SpList,
                                        data.frame(Community = colnames(nature.comm.matrix.t)[j],
                                             Species = comm.occur[comm.occur$Occur >0,]$Species))
}

##2. stack present range ----------------
current.rast.dir <- "."
current.file  <- list.files(current.rast.dir)
current.file.terr <- current.file[current.file %in% paste0(terrestrial.large.sp.exclude.marine, ".tif")]
length(current.file.terr)

current <- rast(list(rast(paste0(current.rast.dir, current.file.terr))))

current.comm.matrix <- data.frame()
for (i in 1:length(Buff.10km.cov$Community)) {
  print(i)
  current.comm.matrix <- rbind(current.comm.matrix, 
                              data.frame(Community = Comm.10km.phyla[i,]$Community,
                                         extract(current, Comm.10km.phyla[i,], method = "simple")))  
}

current.comm.matrix.singe <- current.comm.matrix[!duplicated(current.comm.matrix$Community),]
current.comm.matrix.t <- data.frame(t(current.comm.matrix.singe[,-c(1:2)]))
colnames(current.comm.matrix.t) <- current.comm.matrix.singe$Community

current.comm.SpList <- data.frame()
for ( j in 1:ncol(current.comm.matrix.t)) {
  comm.occur <- data.frame(current.comm.matrix.t[,j],
                           Species =  rownames(current.comm.matrix.t))
  colnames(comm.occur) <- c("Occur", "Species")
  current.comm.SpList <- rbind(current.comm.SpList,
                                      data.frame(Community = colnames(current.comm.matrix.t)[j],
                                                 Species = comm.occur[comm.occur$Occur >0,]$Species))
}

##3. species richness loss estimated by differences in species occurrence between present-natural communities and current communities-----
LQ_SR_anomaly  <- data.frame(Community = colnames(nature.comm.matrix.t), 
                            N.sp.nature = colSums(nature.comm.matrix.t),
                            N.sp.current = colSums(current.comm.matrix.t))
LQ_SR_anomaly$SR_LQ_anomaly  <- LQ_SR_anomaly$N.sp.nature - LQ_SR_anomaly$N.sp.current

##4. create species loss or gain in current communities relative to present-natural communities----
LQ_anomaly  <- data.frame(current.comm.matrix.t - nature.comm.matrix.t)
colnames(LQ_anomaly ) <- current.comm.matrix.singe$Community

LQ.loss.gain.comm.SpList <- data.frame()
for (j in 1:ncol(LQ_anomaly )) {
  comm.occur <- data.frame(LQ_anomaly [,j],
                           Species =  rownames(LQ_anomaly ))
  colnames(comm.occur) <- c("Occur", "Species")
  
  #no change
  if (length(comm.occur[comm.occur$Occur ==  1 | comm.occur$Occur ==  -1,]$Species) == 0) {
    LQ.loss.gain.comm.SpList <- rbind(LQ.loss.gain.comm.SpList,
                                      data.frame(Community = colnames(LQ_anomaly )[j],
                                                 Species = "",
                                                 LossGain = "no change"))
  } else if (length(comm.occur[comm.occur$Occur ==  1,]$Species) == 0) {
    #only loss
    LQ.loss.gain.comm.SpList <- rbind(LQ.loss.gain.comm.SpList,
                                      data.frame(Community = colnames(LQ_anomaly )[j],
                                                 Species = comm.occur[comm.occur$Occur == -1,]$Species,
                                                 LossGain = "loss"))
  } else {
    
    LQ.loss.gain.comm.SpList <- rbind(LQ.loss.gain.comm.SpList,
                                      data.frame(Community = colnames(LQ_anomaly )[j],
                                                 Species = comm.occur[comm.occur$Occur ==  1,]$Species,
                                                 LossGain = "gain"),
                                      data.frame(Community = colnames(LQ_anomaly )[j],
                                                 Species = comm.occur[comm.occur$Occur == -1,]$Species,
                                                 LossGain = "loss"))
  }
}

save(nature.present.comm.SpList, current.comm.SpList, 
     LQ_SR_anomaly , LQ.loss.gain.comm.SpList, 
     file = "Data/PHYLACINE_LQ_sp_gain_loss_no_ocean_no_homo_410community.RData")

