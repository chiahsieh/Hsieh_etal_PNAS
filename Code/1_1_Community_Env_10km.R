#
#  This R script is used to extract Environmental covariates for 410 mammal communities from Rowan et al. (2020) PNAS
#
#  Extracting environmental variables by 10-km radius buffer of community coordinates
#
#  1. Temperature variability since 3.3 Mya and contemporary temperature seasonality 
#      a. BIO1 over twelve time periods from PaleoClim with the resolution of 2.5 mins
#          obtain from http://www.paleoclim.org/
#
#
#  2. Productivity and Topographic complexity were extracted by Google Earth Engine
#      a. EVI values extracted from MODIS 16-day vegetation index (MOD13Q1) during 2000-2016, performed by Google Earth Engine
#            https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD13Q1
#            Code: GGE_10kmBuff_EVI_and_SRTMDEM.txt (extracted with Rowan_2020_Community_10kmBuff.shp)
#            Input: Rowan_Comm_10kmBuff_EVI_2000-2016.csv
#      b. Elevation range with SRTM DEM, performed by Google Earth Engine
#           https://developers.google.com/earth-engine/datasets/catalog/CGIAR_SRTM90_V4
#            Code: GGE_10kmBuff_EVI_and_SRTMDEM.txt (extracted with Rowan_2020_Community_10kmBuff.shp)
#            Input: Rowan_Comm_10kmBuff_EVI_2000-2016.csv
#
#  3. Habitat heterogeneity 
#      a. Shannon index of land-cover from Jung et al. 2020 Scientific Data, IUCN level 2 habitat map
#            Input: iucn_habitatclassification_composite_lvl2_ver004.tif
#               obtain from https://doi.org/10.1038/s41597-020-00599-8



rm(list=ls())
library(readxl)
library(terra); library(sf)
library(dplyr)

# I. Create Spatial data of 410 communities with buffer ==================================
## 1. Community checklists ------------------
###  obtain community checklist from Rowan et al. 2020 PNAS https://doi.org/10.1073/pnas.191048911
setwd(".")
Afro <- read_xlsx("pnas.1910489116.sd01.xlsx", sheet = "Env Data")
Afro <- Afro[,c("Community","Country", "Long", "Lat")]

Indo <- read_xlsx("pnas.1910489116.sd02.xlsx", sheet = "Env Data")
Indo <- Indo[,c("Community","Country", "Long", "Lat")]

Neo <- read_xlsx("pnas.1910489116.sd04.xlsx", sheet = "Env Data")
Neo <- Neo[,c("Community","Country", "Long", "Lat")]

Comm <- rbind(data.frame(Afro, Region = "Afrotropic"),
              data.frame(Indo, Region = "IndoMalyan"),
              data.frame(Neo, Region = "Neotropic"))

nrow(Comm) #410 communiteis

##2. Convert community data.frame into spatial points ------
Comm.pt <- st_as_sf(Comm, coords = c("Long", "Lat"), crs = st_crs("EPSG:4326")) 

##3. Create 10-km radius buffer by each community location ----
Comm.10km <- st_buffer(Comm.pt, dist = 10000)


# II. Extract climatic covariates  ================
##1. Stack BioClimatic layers ----------
### Access PaleoClim with the resolution at 2.5 mins http://www.paleoclim.org/
#### bio_1: Annual Mean Temperature [°C*10]
#### bio_2: Temperature Seasonality [standard deviation*100]

Pliocene_3.3Ma <- rast(list(rast("./PaleoClim_2_5min/1_Pliocene-M2_ca.3.3Ma/bio_1.tif"),
                            rast("./PaleoClim_2_5min/1_Pliocene-M2_ca.3.3Ma/bio_4.tif")))

Pliocene_3.205Ma <- rast(list(rast("./PaleoClim_2_5min/2_Pliocene-mid-Pliocene warm period_3.205Ma/bio_1.tif"),
                              rast("./PaleoClim_2_5min/2_Pliocene-mid-Pliocene warm period_3.205Ma/bio_4.tif")))

Pleistocene_787ka <- rast(list(rast("./PaleoClim_2_5min/3_Pleistocene-MIS19_ca.787ka/bio_1.tif"),
                               rast("./PaleoClim_2_5min/3_Pleistocene-MIS19_ca.787ka/bio_4.tif")))

Pleistocene_130ka <- rast(list(rast("./PaleoClim_2_5min/4_Pleistocene-Last Interglacial_ca.130ka/bio_1.tif"),
                               rast("./PaleoClim_2_5min/4_Pleistocene-Last Interglacial_ca.130ka/bio_4.tif")))

LGM_21ka <- rast(list(rast("./PaleoClim_2_5min/5_CHELSA-LGM_ca.21ka/bio_1.tif"),
                      rast("./PaleoClim_2_5min/5_CHELSA-LGM_ca.21ka/bio_4.tif")))

Pleistocene_14.7ka <- rast(list(rast("./PaleoClim_2_5min/6_Pleistocene-Heinrich Stadial_17.0-14.7ka/bio_1.tif"),
                                rast("./PaleoClim_2_5min/6_Pleistocene-Heinrich Stadial_17.0-14.7ka/bio_4.tif")))

Pleistocene_12.9ka <- rast(list(rast("./PaleoClim_2_5min/7_Pleistocene-Bølling-Allerød_14.7-12.9ka/bio_1.tif"),
                                rast("./PaleoClim_2_5min/7_Pleistocene-Bølling-Allerød_14.7-12.9ka/bio_4.tif")))

Pleistocene_11.7ka <- rast(list(rast("./PaleoClim_2_5min/8_Pleistocene-YoungerDryasStadial_12.9-11.7ka/bio_1.tif"),
                                rast("./PaleoClim_2_5min/8_Pleistocene-YoungerDryasStadial_12.9-11.7ka/bio_4.tif")))

Pleistocene_8.326ka <- rast(list(rast("./PaleoClim_2_5min/9_Pleistocene-early-Holocene_11.7-8.326ka/bio_1.tif"),
                                 rast("./PaleoClim_2_5min/9_Pleistocene-early-Holocene_11.7-8.326ka/bio_4.tif")))

Pleistocene_4.2ka <- rast(list(rast("./PaleoClim_2_5min/10_Pleistocene-mid-Holocene_8.326-4.2ka/bio_1.tif"),
                               rast("./PaleoClim_2_5min/10_Pleistocene-mid-Holocene_8.326-4.2ka/bio_4.tif")))

Pleistocene_0.3ka <- rast(list(rast("./PaleoClim_2_5min/11_Pleistocene-late-Holocene_4.2-0.3ka/bio_1.tif"),
                               rast("./PaleoClim_2_5min/11_Pleistocene-late-Holocene_4.2-0.3ka/bio_4.tif")))

Mondern_1979_2013 <- rast(list(rast("./PaleoClim_2_5min/12_CHELSA_cur_V1_2B_1979_2013/bio_1.tif"),
                               rast("./PaleoClim_2_5min/12_CHELSA_cur_V1_2B_1979_2013/bio_4.tif")))

PaleoClim <- rbind(data.frame(Community = Comm.10km$Community,
                              round(extract(Pliocene_3.3Ma, Comm.10km, "mean", na.rm=TRUE),0),
                              Time = "Pliocene_3.3Ma"),
                   data.frame(Community = Comm.10km$Community,
                              round(extract(Pliocene_3.205Ma, Comm.10km, "mean", na.rm=TRUE),0),
                              Time = "Pliocene_3.205Ma"),
                   data.frame(Community = Comm.10km$Community,
                              round(extract(Pleistocene_787ka, Comm.10km, "mean", na.rm=TRUE),0),
                              Time = "Pleistocene_787ka"),
                   data.frame(Community = Comm.10km$Community,
                              round(extract(Pleistocene_130ka, Comm.10km, "mean", na.rm=TRUE),0),
                              Time = "Pleistocene_130ka"),
                   data.frame(Community = Comm.10km$Community,
                              round(extract(LGM_21ka, Comm.10km, "mean", na.rm=TRUE),0),
                              Time = "LGM_21ka"),
                   data.frame(Community = Comm.10km$Community,
                              round(extract(Pleistocene_14.7ka, Comm.10km, "mean", na.rm=TRUE),0),
                              Time = "Pleistocene_14.7ka"),
                   data.frame(Community = Comm.10km$Community,
                              round(extract(Pleistocene_12.9ka, Comm.10km, "mean", na.rm=TRUE),0),
                              Time = "Pleistocene_12.9ka"),
                   data.frame(Community = Comm.10km$Community,
                              round(extract(Pleistocene_11.7ka, Comm.10km, "mean", na.rm=TRUE),0),
                              Time = "Pleistocene_11.7ka"),
                   data.frame(Community = Comm.10km$Community,
                              round(extract(Pleistocene_8.326ka, Comm.10km, "mean", na.rm=TRUE),0),
                              Time = "Pleistocene_8.326ka"),
                   data.frame(Community = Comm.10km$Community,
                              round(extract(Pleistocene_4.2ka, Comm.10km, "mean", na.rm=TRUE),0),
                              Time = "Pleistocene_4.2ka"),
                   data.frame(Community = Comm.10km$Community,
                              round(extract(Pleistocene_0.3ka, Comm.10km, "mean", na.rm=TRUE),0),
                              Time = "Pleistocene_0.3ka"),
                   data.frame(Community = Comm.10km$Community,
                              round(extract(Mondern_1979_2013, Comm.10km, "mean", na.rm=TRUE),0),
                              Time = "Mondern_1979_2013"))


##2. Temperature variability since 3.3 Mya ---------------
Temp_deep <- PaleoClim %>% group_by(Community) %>%
  summarize(Temp_deep_Mean = round(mean(bio_1/10, na.rm = T),2), Temp_deep_SD = round(sd(bio_1/10, na.rm = T),2))

##3. Temperature seasonality 1979-2013  ---------------
Temp_present <- data.frame(Community = PaleoClim[PaleoClim$Time == "Mondern_1979_2013", ]$Community,
                               present_temp = round(PaleoClim[PaleoClim$Time == "Mondern_1979_2013", ]$bio_1/10,2),
                               present_temp_sea = round(PaleoClim[PaleoClim$Time == "Mondern_1979_2013", ]$bio_4/100,2))

Temp_deep <- left_join(Temp_deep, Temp_present, join_by(Community == Community))


#III. Extract habitat diversity =================================
### IUCN habitat classification Jung et al. 2020 Scientific Data https://doi.org/10.1038/s41597-020-00599-8
### iucn_habitatclassification_composite_lvl2_ver004.tif obtained from https://doi.org/10.1038/s41597-020-00599-8
iucn.habitat <- rast("./iucn_habitatclassification_composite_lvl2_ver004.tif")

dt_total.land <- c()
for (i in 1: length(Comm.10km$Community)){
  habitat.buf.freq = data.frame(table(extract(iucn.habitat, Comm.10km[i,]))); colnames(habitat.buf.freq) <- c("ID","Habitat","N")
  #Calculate the number of habitat types
  no.habi = data.frame(No_Habitat=length(unique((habitat.buf.freq$Habitat))))
  #Calculate Shannon diversity 
  habitat.buf.freq.p = habitat.buf.freq
  habitat.buf.freq.p$plnp = (habitat.buf.freq.p$N/sum(habitat.buf.freq.p$N))*(log(habitat.buf.freq.p$N/sum(habitat.buf.freq.p$N))) #log is the natural logarithm in R
  shannon = data.frame(Shannon = round(-sum(habitat.buf.freq.p$plnp),3)) 
  dt_total.land <- rbind(dt_total.land,
                        cbind(Comm.10km[i,]$Community, no.habi, shannon))
} 
colnames(dt_total.land) <- c("Community", "No_Habitat", "Shannon")

#IV. EVI annual mean and seasonality obtained from Google Earth Engine as mean values within 10km-buffer =======
## monthly averages of Enhanced Vegetation Index stacked by months during 2000-2016 to indicate long-term productivity
evi <-read.csv("Data/Rowan_Comm_10kmBuff_EVI_2000-2016.csv", encoding="UTF-8")
evi$EVI_month.mean <- round(apply(evi[,2:13], 1, mean)/10000,2)
evi$EVI_seasonality <- round(apply(evi[,2:13], 1, sd)/apply(evi[,2:13], 1, mean),2)
## check communities with productivity below zero, which indicating low vegetation coverage and plant-related productivity
evi[evi$EVI_month.mean <0,]$Community
evi <- evi[,c("Community", "EVI_month.mean", "EVI_seasonality")]

#V.  Elevation obtained from Google Earth Engine with mean, maximum, and minimum cell values within 10km-buffer =======
elevation <-read.csv("Data/Rowan_Comm_10kmBuff_elevation.csv", encoding="UTF-8") 
View(elevation[!elevation$Community %in% Comm$Community,])
elevation$range <- elevation$max - elevation$min
elevation <- elevation[, c("Community", "mean", "range")]
colnames(elevation) <- c("Community", paste0("Elev_", colnames(elevation[,-1])))


# VI. Join all information ==================
Buff.10km.cov <- left_join(Comm, Temp_deep, join_by(Community == Community))
Buff.10km.cov <- left_join(Buff.10km.cov, dt_total.land, join_by(Community == Community))

#unify community name with special character 
Buff.10km.cov$Community_unify <- Buff.10km.cov$Community
##check community name with special character to add unified community name
Buff.10km.cov[!Buff.10km.cov$Community %in% evi$Community,]
Buff.10km.cov[!Buff.10km.cov$Community %in% evi$Community,]$Community_unify <- 
  c("Qui__ma", "Hoollongapar_","Wilpattu_", "Mamirau_")

evi$Community_unify <- evi$Community
#check community name with special character to add unified community name
evi[!evi$Community %in% Buff.10km.cov$Community,]$Community_unify
evi[!evi$Community %in% Buff.10km.cov$Community,]$Community_unify <- 
  c("Hoollongapar_","Mamirau_", "Qui__ma", "Wilpattu_" )

elevation$Community_unify <- elevation$Community
#check community name with special character to add unified community name
elevation[!elevation$Community %in% Buff.10km.cov$Community,]$Community_unify
#add unified community name in correct order 
elevation[!elevation$Community %in% Buff.10km.cov$Community,]$Community_unify <- 
  c("Hoollongapar_","Mamirau_", "Qui__ma", "Wilpattu_" )

Buff.10km.cov <- left_join(Buff.10km.cov, evi[,-1], join_by(Community_unify == Community_unify))
Buff.10km.cov <- left_join(Buff.10km.cov, elevation[,-1], join_by(Community_unify == Community_unify))

#subset environmental covariates for subsequent analyses 
Buff.10km.cov.analyzed <- Buff.10km.cov[, c("Community","Community_unify","Region", "Long", "Lat", "Temp_deep_Mean", "Temp_deep_SD",
                                            "present_temp", "present_temp_sea", "EVI_month.mean", "EVI_seasonality", 
                                            "Elev_mean", "Elev_range","Shannon")]

#three communities with negative EVI. Exclude in later step
View(Buff.10km.cov.analyzed[Buff.10km.cov.analyzed$EVI_month.mean <0,])

#Exclude communities with NA 
Buff.10km.cov.analyzed[is.na(Buff.10km.cov.analyzed$Temp_deep_SD) == T,]
Buff.10km.cov.analyzed[is.na(Buff.10km.cov.analyzed$present_temp) == T,]
Buff.10km.cov.analyzed[is.na(Buff.10km.cov.analyzed$EVI_month.mean) == T,]
Buff.10km.cov.analyzed[is.na(Buff.10km.cov.analyzed$Elev_range) == T,]
Buff.10km.cov.analyzed[is.na(Buff.10km.cov.analyzed$Shannon) == T,]

#407 of 410 communities has sufficient environmental covariates
nrow(Buff.10km.cov.analyzed)
write.csv(Buff.10km.cov.analyzed, "Data/Rowan_Comm_10kmBuff_Env_cov.csv", row.names = F)

