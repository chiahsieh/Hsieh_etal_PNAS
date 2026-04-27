#
#   Supplementary Figures S1 to S3 to show predation matrix and eleven ecological traits
#
#  
#   Input 
#     1. Community408_carni_web_doc_inf95_500g_SR3_commcollect.Rdata
#
#   Trait and Taxonomy from Soria et al. (2021) COMBINE dataset https://doi.org/10.1002/ecy.3344 to have fine-scale trait information to characterize ecological similarity
#     1. trait_data_imputed.csv
#
#

library(ggtree)
library(ape)
library(phytools)
library(phytools)
library(stringr); library(dplyr)
library(ggplot2); library(RColorBrewer)

setwd("/Users/chiahsieh/Documents/Research/Projects/Tropical_Foodweb/Data/Community/SR3_carnivora_predator_500g")

load("/Users/chiahsieh/Documents/Research/R/Tropical_foodweb/Carnivore_only_500g/PNAS_DataCodeRepo/Full_data/Community408_carni_web_doc_inf95_500g_SR3_commcollect.Rdata")
#collect food web species list 
comm.trophic.link.inf95.10p <- c()
for (i in 1:length(Tropic.carni.webs.inf95.10p)) {
  #Collect all interactions in a local food web  
  comm.trophic.link.inf95.10p <- rbind(comm.trophic.link.inf95.10p, 
                                       data.frame(Tropic.carni.webs.inf95.10p[[i]][3],
                                                  Community = Tropic.carni.webs.inf95.10p[[i]][2]$properties$Community_unify,
                                                  Region = Tropic.carni.webs.inf95.10p[[i]][2]$properties$Region))
}

realm.intx <- comm.trophic.link.inf95.10p
colnames(realm.intx) <- c("resource" , "consumer" , "trophic.links.documeted" ,"Community" , "Region")
realm.intx$resource <- gsub(" ", "_", realm.intx$resource)
realm.intx$consumer <- gsub(" ", "_", realm.intx$consumer)

region.sp <- unique(c(realm.intx$resource,realm.intx$consumer))
region.c.sp <- unique(realm.intx[, c("Region", "consumer")])
region.r.sp <- unique(realm.intx[, c("Region", "resource")])

mamm.func <- read.csv("/Users/chiahsieh/Documents/Research/Data/Trait/Mammal_Trait/COMBINE/trait_data_imputed.csv")
mamm.func$phylacine_binomial <- gsub(" ", "_", mamm.func$phylacine_binomial )
mamm.func.unique <- mamm.func[!duplicated(mamm.func$phylacine_binomial),]
mamm.func.study <- mamm.func.unique[mamm.func.unique$phylacine_binomial %in% region.sp, ]
mamm.func.study$adult_mass_g_log <- log(mamm.func.study$adult_mass_g,10)
mamm.func.study$brain_mass_g_log <- log(mamm.func.study$brain_mass_g,10)
mamm.func.study$adult_body_length_mm_log <- log(mamm.func.study$adult_body_length_mm,10)
mamm.func.study$dispersal_km_log <- log(mamm.func.study$dispersal_km,10)
mamm.func.study$foraging_stratum <- factor(mamm.func.study$foraging_stratum)
mamm.func.study$activity_cycle <- factor(mamm.func.study$activity_cycle, levels =c(1, 2, 3), labels =c("Noct. only", "Inter.", "Diur. only"))
mamm.func.study$trophic_level <- factor(mamm.func.study$trophic_level, labels =c("Herbivore", "Omnivore", "Carnivore"))
mamm.func.study$fossoriality <-ifelse(mamm.func.study$fossoriality == 2, 0, 1)

#Fill values from IUCN habitat types 
mamm.func.study[mamm.func.study$phylacine_binomial == "Cebus_albifrons",]$habitat_breadth_n <- 2 #https://www.iucnredlist.org/species/39951/191703935#habitat-ecology
mamm.func.study[mamm.func.study$phylacine_binomial == "Cebus_brunneus",]$habitat_breadth_n <- 1 #https://www.iucnredlist.org/species/81237954/17981252#habitat-ecology
mamm.func.study[mamm.func.study$phylacine_binomial == "Cebus_capucinus",]$habitat_breadth_n <- 2 #https://www.iucnredlist.org/species/81257277/191708164#habitat-ecology
mamm.func.study[mamm.func.study$phylacine_binomial == "Loxodonta_africana",]$habitat_breadth_n <- 6 #https://www.iucnredlist.org/species/181008073/223031019#habitat-ecology
mamm.func.study[mamm.func.study$phylacine_binomial == "Trachypithecus_barbei",]$habitat_breadth_n <- 1 #https://www.iucnredlist.org/species/41554/17960144#habitat-ecology

mamm.func.study.hypervolume <- mamm.func.study[,c("phylacine_binomial","adult_mass_g_log","brain_mass_g_log", "adult_body_length_mm_log", "dispersal_km_log",
                                                  "det_diet_breadth_n", "habitat_breadth_n", "hibernation_torpor", "fossoriality",
                                                  "foraging_stratum", "activity_cycle","trophic_level")]

region.c.sp.trait <- left_join(region.c.sp, mamm.func.study.hypervolume, join_by(consumer == phylacine_binomial))
region.r.sp.trait <- left_join(region.r.sp, mamm.func.study.hypervolume, join_by(resource == phylacine_binomial))

#I. Neotropic -------------------------------------------------------------
neo.pair <- unique(realm.intx[realm.intx$Region == "Neotropic",c("resource", "consumer")])
occ.splist.neo <- data.frame()
for (i in 1:length(unique(neo.pair$consumer))) {
  dt <- data.frame(resource = sort(unique(neo.pair$resource)),
                   consumer =  sort(unique(neo.pair$consumer))[i],
                   occur = 0)
  dt[dt$resource %in% neo.pair[neo.pair$consumer == sort(unique(neo.pair$consumer))[i],]$resource, ]$occur <- 1
  occ.splist.neo <- rbind(occ.splist.neo, dt)
}
occ.splist.neo$resource <- factor(occ.splist.neo$resource, levels = sort(unique(occ.splist.neo$resource),decreasing = T))
occ.splist.neo$occur <- factor(occ.splist.neo$occur, levels = c(0,1))

## 1. predation matrix ------
neo.pair.potential <- unique(realm.intx[realm.intx$Region == "Neotropic",c("resource", "consumer", "trophic.links.documeted")])
neo.pair.potential <- neo.pair.potential[neo.pair.potential$trophic.links.documeted == 0,]
neo.pair.potential$resource.consumer <- paste0(neo.pair.potential$resource, ".", neo.pair.potential$consumer)

occ.splist.neo.potential <- occ.splist.neo
occ.splist.neo.potential$occur <- as.numeric(as.character(occ.splist.neo.potential$occur))
occ.splist.neo.potential$resource.consumer <- paste0(occ.splist.neo.potential$resource, ".", occ.splist.neo.potential$consumer)
occ.splist.neo.potential[occ.splist.neo.potential$resource.consumer %in% neo.pair.potential$resource.consumer,]$occur <- 2
occ.splist.neo.potential$occur <- factor(occ.splist.neo.potential$occur, levels = c(0,1,2), labels = c("Absent", "Documented", "Potential"))

occur.neo.potential <- ggplot(occ.splist.neo.potential, aes(resource, consumer)) +  
  geom_point(aes(fill=occur), pch = 21, size = 3) + 
  scale_fill_manual(values = c("White", "Black", "Grey"), name="Predation") +
  coord_flip() + 
  theme_minimal()  +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                            axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) + 
  labs(x= "", y = "") 

## 2. body mass  ------
col.body <- c("yellow", "orange", "red", "darkred")
mass.neo <- rbind(data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Neotropic",]$resource,
                             bodymass = region.r.sp.trait[region.r.sp.trait$Region == "Neotropic",]$adult_mass_g_log,
                             Type = "body"),
                  data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Neotropic",]$resource,
                             bodymass = region.r.sp.trait[region.r.sp.trait$Region == "Neotropic",]$brain_mass_g_log,
                             Type = "brain"))
mass.neo$resource.ordered <- factor(mass.neo$resource, levels = sort(unique(mass.neo$resource),decreasing = T))
mass.neo.plot <- ggplot(mass.neo, aes(resource.ordered, 1, fill=bodymass)) +  
  facet_wrap(~Type) + 
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  geom_tile(color="white") +
  scale_fill_gradientn(colours = col.body,name = "log10-mass (g)") +
  coord_flip() +  theme_minimal() +theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")
## 3. body length  ------
col.length <- c("pink","#ffb6db","#ff6db6", "violetred")
length.neo <- rbind(data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Neotropic",]$resource,
                               length = region.r.sp.trait[region.r.sp.trait$Region == "Neotropic",]$adult_body_length_mm_log))
length.neo$resource.ordered <- factor(length.neo$resource, levels = sort(unique(length.neo$resource),decreasing = T))
length.neo.plot <- ggplot(length.neo, aes(resource.ordered, 1, fill=length)) +  
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  geom_tile(color="white") +
  scale_fill_gradientn(colours = col.length,name = "log10-body\nlength (mm)") +
  coord_flip() +  theme_minimal() +theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")
## 4. dispersal distance ------
col.dist <- c("lavender","#b66dff","mediumpurple4", "#490092")
dist.neo <- rbind(data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Neotropic",]$resource,
                             dist = region.r.sp.trait[region.r.sp.trait$Region == "Neotropic",]$dispersal_km_log))
dist.neo$resource.ordered <- factor(dist.neo$resource, levels = sort(unique(dist.neo$resource),decreasing = T))
dist.neo.plot <- ggplot(dist.neo, aes(resource.ordered, 1, fill=dist)) +  
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  geom_tile(color="white") +
  scale_fill_gradientn(colours = col.dist,name = "log10-dispersal\ndistance (km)") +
  coord_flip() +  theme_minimal() +theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")

## 5. dietary and habitat breadth ------
breadth.neo <- rbind(data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Neotropic",]$resource,
                                Breadth = region.r.sp.trait[region.r.sp.trait$Region == "Neotropic",]$habitat_breadth_n,
                                Type = "habitat"),
                     data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Neotropic",]$resource,
                                Breadth = region.r.sp.trait[region.r.sp.trait$Region == "Neotropic",]$det_diet_breadth_n, 
                                Type = "diet"))
breadth.neo$resource.ordered <- factor(breadth.neo$resource, levels = sort(unique(breadth.neo$resource),decreasing = T))
col <- c("#BCE4D8","#7BBFC9", "#52A8BC", "#52A8BC", "#3690AE", "#3480A2","#317197","#2E628C","#2C5985")
breadth.neo.plot <- ggplot(breadth.neo, aes(resource.ordered, 1, fill=Breadth)) +  
  facet_wrap(~Type) +
  geom_tile(color="white") +
  scale_fill_gradientn(colours = col, name = "Breadth") +
  coord_flip() +  theme_minimal() +theme(strip.text = element_blank(),
                                         panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")

## 6. hibernation and fossor  -----
col.bin <- c("grey80", "grey20")
hib.fossor.neo <- rbind(data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Neotropic",]$resource,
                                   bin = region.r.sp.trait[region.r.sp.trait$Region == "Neotropic",]$hibernation_torpor,
                                   Type = "hiber"),
                        data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Neotropic",]$resource,
                                   bin = region.r.sp.trait[region.r.sp.trait$Region == "Neotropic",]$fossoriality,
                                   Type = "fossor"))
hib.fossor.neo$resource.ordered <- factor(hib.fossor.neo$resource, levels = sort(unique(hib.fossor.neo$resource),decreasing = T))
hib.fossor.neo$bin <- factor(hib.fossor.neo$bin, levels = c(0,1))
bin.plot <- ggplot(hib.fossor.neo, aes(resource.ordered, 1, fill=bin)) +  
  facet_wrap(~Type) +
  geom_tile(color="white") +
  scale_fill_manual(values = c("grey80", "grey20"), name = "Binary") +
  coord_flip() +  theme_minimal() +theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")
## 7. activity pattern -----
activity.neo <- data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Neotropic",]$resource,
                           activity_cycle = region.r.sp.trait[region.r.sp.trait$Region == "Neotropic",]$activity_cycle)
activity.neo$resource.ordered <- factor(activity.neo$resource, levels = sort(activity.neo$resource,decreasing = T))
activity.neo.plot <- ggplot(activity.neo, aes(resource.ordered, 1, fill=activity_cycle)) +  
  geom_tile(color="white") +
  scale_fill_manual(values = c("dodgerblue","wheat2", "goldenrod"), name = "Activity pattern") +
  coord_flip() +  theme_minimal() +theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")
## 8. forage strata -----
forage.neo <- data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Neotropic",]$resource,
                         forage = region.r.sp.trait[region.r.sp.trait$Region == "Neotropic",]$foraging_stratum)
forage.neo$resource.ordered <- factor(forage.neo$resource, levels = sort(forage.neo$resource,decreasing = T))
forage.neo.plot <- ggplot(forage.neo, aes(resource.ordered, 1, fill= forage)) +  
  geom_tile(color="white") +
  scale_fill_manual(values = c("tan4","olivedrab3", "green4"), name = "Forage stratum") +
  coord_flip() +  theme_minimal() +theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")
## 9. trophic level -----
level.neo <- data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Neotropic",]$resource,
                        trophic_level = region.r.sp.trait[region.r.sp.trait$Region == "Neotropic",]$trophic_level)
level.neo$resource.ordered <- factor(level.neo$resource, levels = sort(level.neo$resource,decreasing = T))
level.neo.plot <- ggplot(level.neo, aes(resource, 1, fill=trophic_level)) +  
  geom_tile(color="white") +
  scale_fill_manual(values = c("green4","tan3", "red3"), name = "Trophic level") +
  coord_flip() +  theme_minimal() +theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")
## 10. Export to combine --------
trait.neo <- ggarrange(mass.neo.plot +  theme(legend.position = "none", strip.text = element_blank()), 
                       length.neo.plot + rremove("y.text") +  theme(legend.position = "none"), 
                       dist.neo.plot + rremove("y.text") +  theme(legend.position = "none"), 
                       breadth.neo.plot + rremove("y.text") +  theme(legend.position = "none", strip.text = element_blank()), 
                       bin.plot + rremove("y.text") +  theme(legend.position = "none", strip.text = element_blank()), 
                       activity.neo.plot + rremove("y.text") +  theme(legend.position = "none"), 
                       forage.neo.plot + rremove("y.text") +  theme(legend.position = "none"), 
                       level.neo.plot + rremove("y.text") +  theme(legend.position = "none"), 
                       nrow = 1, widths = c(.6,.18,.18,.3,.3,.18,.18,.18)) + bgcolor("white") + border("white")

trait.neo.legend <- ggarrange(mass.neo.plot  , 
                              length.neo.plot + rremove("y.text") , 
                              dist.neo.plot + rremove("y.text") , 
                              breadth.neo.plot + rremove("y.text") , 
                              bin.plot + rremove("y.text"), 
                              nrow = 1, widths = c(.6,.18,.18,.3,.3)) + bgcolor("white") + border("white")

trait.neo.legend <- ggarrange(mass.neo.plot  , 
                              activity.neo.plot + rremove("y.text") , 
                              forage.neo.plot + rremove("y.text") , 
                              level.neo.plot + rremove("y.text"), 
                              nrow = 1, widths = c(.6,.18,.18,.18)) + bgcolor("white") + border("white")


#II. Afrotropic -------------------------------------------------------------
afr.pair <- unique(realm.intx[realm.intx$Region == "Afrotropic",c("resource", "consumer")])
occ.splist.afr <- data.frame()
for (i in 1:length(unique(afr.pair$consumer))) {
  dt <- data.frame(resource = sort(unique(afr.pair$resource)),
                   consumer =  sort(unique(afr.pair$consumer))[i],
                   occur = 0)
  dt[dt$resource %in% afr.pair[afr.pair$consumer == sort(unique(afr.pair$consumer))[i],]$resource, ]$occur <- 1
  occ.splist.afr <- rbind(occ.splist.afr, dt)
}
occ.splist.afr$resource <- factor(occ.splist.afr$resource, levels = sort(unique(occ.splist.afr$resource),decreasing = T))
occ.splist.afr$occur <- factor(occ.splist.afr$occur, levels = c(0,1))

afr.pair.potential <- unique(realm.intx[realm.intx$Region == "Afrotropic",c("resource", "consumer", "trophic.links.documeted")])
afr.pair.potential <- afr.pair.potential[afr.pair.potential$trophic.links.documeted == 0,]
afr.pair.potential$resource.consumer <- paste0(afr.pair.potential$resource, ".", afr.pair.potential$consumer)

occ.splist.afr.potential <- occ.splist.afr
occ.splist.afr.potential$occur <- as.numeric(as.character(occ.splist.afr.potential$occur))
occ.splist.afr.potential$resource.consumer <- paste0(occ.splist.afr.potential$resource, ".", occ.splist.afr.potential$consumer)
occ.splist.afr.potential[occ.splist.afr.potential$resource.consumer %in% afr.pair.potential$resource.consumer,]$occur <- 2
occ.splist.afr.potential$occur <- factor(occ.splist.afr.potential$occur, levels = c(0,1,2))

## 1. predation matrix ------
occur.afr.potential <- ggplot(occ.splist.afr.potential, aes(resource, consumer)) +  
  geom_point(aes(fill=occur), pch = 21, size = 3) + 
  scale_fill_manual(values = c("White", "Black", "Grey"), name="Predation") +
  coord_flip() + 
  theme_minimal()  +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                            axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) + 
  labs(x= "", y = "") 

## 2. body mass  ------
col.body <- c("yellow", "orange", "red", "darkred")
mass.afr <- rbind(data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Afrotropic",]$resource,
                             bodymass = region.r.sp.trait[region.r.sp.trait$Region == "Afrotropic",]$adult_mass_g_log,
                             Type = "body"),
                  data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Afrotropic",]$resource,
                             bodymass = region.r.sp.trait[region.r.sp.trait$Region == "Afrotropic",]$brain_mass_g_log,
                             Type = "brain"))
mass.afr$resource.ordered <- factor(mass.afr$resource, levels = sort(unique(mass.afr$resource),decreasing = T))
mass.afr.plot <- ggplot(mass.afr, aes(resource.ordered, 1, fill=bodymass)) +  
  facet_wrap(~Type) + 
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  geom_tile(color="white") +
  scale_fill_gradientn(colours = col.body,name = "log10-mass (g)") +
  coord_flip() +  theme_minimal() +theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")
## 3. body length  ------
col.length <- c("pink","#ffb6db","#ff6db6", "violetred")
length.afr <- rbind(data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Afrotropic",]$resource,
                               length = region.r.sp.trait[region.r.sp.trait$Region == "Afrotropic",]$adult_body_length_mm_log))
length.afr$resource.ordered <- factor(length.afr$resource, levels = sort(unique(length.afr$resource),decreasing = T))
length.afr.plot <- ggplot(length.afr, aes(resource.ordered, 1, fill=length)) +  
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  geom_tile(color="white") +
  scale_fill_gradientn(colours = col.length,name = "log10-body\nlength (mm)") +
  coord_flip() +  theme_minimal() +theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")
## 4. dispersal distance ------
col.dist <- c("lavender","#b66dff","mediumpurple4", "#490092")
dist.afr <- rbind(data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Afrotropic",]$resource,
                             dist = region.r.sp.trait[region.r.sp.trait$Region == "Afrotropic",]$dispersal_km_log))
dist.afr$resource.ordered <- factor(dist.afr$resource, levels = sort(unique(dist.afr$resource),decreasing = T))
dist.afr.plot <- ggplot(dist.afr, aes(resource.ordered, 1, fill=dist)) +  
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  geom_tile(color="white") +
  scale_fill_gradientn(colours = col.dist,name = "log10-dispersal\ndistance (km)") +
  coord_flip() +  theme_minimal() +theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")


## 5. dietary and habitat breadth ------
breadth.afr <- rbind(data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Afrotropic",]$resource,
                                Breadth = region.r.sp.trait[region.r.sp.trait$Region == "Afrotropic",]$habitat_breadth_n,
                                Type = "habitat"),
                     data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Afrotropic",]$resource,
                                Breadth = region.r.sp.trait[region.r.sp.trait$Region == "Afrotropic",]$det_diet_breadth_n, 
                                Type = "diet"))
breadth.afr$resource.ordered <- factor(breadth.afr$resource, levels = sort(unique(breadth.afr$resource),decreasing = T))
col <- c("#BCE4D8","#7BBFC9", "#52A8BC", "#52A8BC", "#3690AE", "#3480A2","#317197","#2E628C","#2C5985")
breadth.afr.plot <- ggplot(breadth.afr, aes(resource.ordered, 1, fill=Breadth)) +  
  facet_wrap(~Type) +
  geom_tile(color="white") +
  scale_fill_gradientn(colours = col, name = "Breadth") +
  coord_flip() +  theme_minimal() +theme(strip.text = element_blank(),
                                         panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")

## 6. hibernation and fossor  -----
col.bin <- c("grey80", "grey20")
hib.fossor.afr <- rbind(data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Afrotropic",]$resource,
                                   bin = region.r.sp.trait[region.r.sp.trait$Region == "Afrotropic",]$hibernation_torpor,
                                   Type = "hiber"),
                        data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Afrotropic",]$resource,
                                   bin = region.r.sp.trait[region.r.sp.trait$Region == "Afrotropic",]$fossoriality,
                                   Type = "fossor"))
hib.fossor.afr$resource.ordered <- factor(hib.fossor.afr$resource, levels = sort(unique(hib.fossor.afr$resource),decreasing = T))
hib.fossor.afr$bin <- factor(hib.fossor.afr$bin, levels = c(0,1))
bin.plot <- ggplot(hib.fossor.afr, aes(resource.ordered, 1, fill=bin)) +  
  facet_wrap(~Type) +
  geom_tile(color="white") +
  scale_fill_manual(values = c("grey80", "grey20"), name = "Binary") +
  coord_flip() +  theme_minimal() +theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")
## 7. activity pattern -----
activity.afr <- data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Afrotropic",]$resource,
                           activity_cycle = region.r.sp.trait[region.r.sp.trait$Region == "Afrotropic",]$activity_cycle)
activity.afr$resource.ordered <- factor(activity.afr$resource, levels = sort(activity.afr$resource,decreasing = T))
activity.afr.plot <- ggplot(activity.afr, aes(resource.ordered, 1, fill=activity_cycle)) +  
  geom_tile(color="white") +
  scale_fill_manual(values = c("dodgerblue","wheat2", "goldenrod"), name = "Activity pattern") +
  coord_flip() +  theme_minimal() +theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")
## 8. forage strata -----
forage.afr <- data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Afrotropic",]$resource,
                         forage = region.r.sp.trait[region.r.sp.trait$Region == "Afrotropic",]$foraging_stratum)
forage.afr$resource.ordered <- factor(forage.afr$resource, levels = sort(forage.afr$resource,decreasing = T))
forage.afr.plot <- ggplot(forage.afr, aes(resource.ordered, 1, fill= forage)) +  
  geom_tile(color="white") +
  scale_fill_manual(values = c("tan4","olivedrab3", "green4"), name = "Forage stratum") +
  coord_flip() +  theme_minimal() +theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")
## 9. trophic level -----
level.afr <- data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Afrotropic",]$resource,
                        trophic_level = region.r.sp.trait[region.r.sp.trait$Region == "Afrotropic",]$trophic_level)
level.afr$resource.ordered <- factor(level.afr$resource, levels = sort(level.afr$resource,decreasing = T))
level.afr.plot <- ggplot(level.afr, aes(resource, 1, fill=trophic_level)) +  
  geom_tile(color="white") +
  scale_fill_manual(values = c("green4","tan3", "red3"), name = "Trophic level") +
  coord_flip() +  theme_minimal() +theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")


## 10. Export to combine --------
trait.afr <- ggarrange(mass.afr.plot +  theme(legend.position = "none", strip.text = element_blank()), 
                       length.afr.plot + rremove("y.text") +  theme(legend.position = "none"), 
                       dist.afr.plot + rremove("y.text") +  theme(legend.position = "none"), 
                       breadth.afr.plot + rremove("y.text") +  theme(legend.position = "none", strip.text = element_blank()), 
                       bin.plot + rremove("y.text") +  theme(legend.position = "none", strip.text = element_blank()), 
                       activity.afr.plot + rremove("y.text") +  theme(legend.position = "none"), 
                       forage.afr.plot + rremove("y.text") +  theme(legend.position = "none"), 
                       level.afr.plot + rremove("y.text") +  theme(legend.position = "none"), 
                       nrow = 1, widths = c(.6,.18,.18,.3,.3,.18,.18,.18)) + bgcolor("white") + border("white")

trait.afr.legend <- ggarrange(mass.afr.plot  , 
                              length.afr.plot + rremove("y.text") , 
                              dist.afr.plot + rremove("y.text") , 
                              breadth.afr.plot + rremove("y.text") , 
                              bin.plot + rremove("y.text"), 
                              nrow = 1, widths = c(.6,.18,.18,.3,.3)) + bgcolor("white") + border("white")

trait.afr.legend <- ggarrange(mass.afr.plot  , 
                              activity.afr.plot + rremove("y.text") , 
                              forage.afr.plot + rremove("y.text") , 
                              level.afr.plot + rremove("y.text"), 
                              nrow = 1, widths = c(.6,.18,.18,.18)) + bgcolor("white") + border("white")



#III. Indomalayan -------------------------------------------------------------
ind.pair <- unique(realm.intx[realm.intx$Region == "Indomalayan",c("resource", "consumer")])
occ.splist.ind <- data.frame()
for (i in 1:length(unique(ind.pair$consumer))) {
  dt <- data.frame(resource = sort(unique(ind.pair$resource)),
                   consumer =  sort(unique(ind.pair$consumer))[i],
                   occur = 0)
  dt[dt$resource %in% ind.pair[ind.pair$consumer == sort(unique(ind.pair$consumer))[i],]$resource, ]$occur <- 1
  occ.splist.ind <- rbind(occ.splist.ind, dt)
}
occ.splist.ind$resource <- factor(occ.splist.ind$resource, levels = sort(unique(occ.splist.ind$resource),decreasing = T))
occ.splist.ind$occur <- factor(occ.splist.ind$occur, levels = c(0,1))

ind.pair.potential <- unique(realm.intx[realm.intx$Region == "Indomalayan",c("resource", "consumer", "trophic.links.documeted")])
ind.pair.potential <- ind.pair.potential[ind.pair.potential$trophic.links.documeted == 0,]
ind.pair.potential$resource.consumer <- paste0(ind.pair.potential$resource, ".", ind.pair.potential$consumer)

occ.splist.ind.potential <- occ.splist.ind
occ.splist.ind.potential$occur <- as.numeric(as.character(occ.splist.ind.potential$occur))
occ.splist.ind.potential$resource.consumer <- paste0(occ.splist.ind.potential$resource, ".", occ.splist.ind.potential$consumer)
occ.splist.ind.potential[occ.splist.ind.potential$resource.consumer %in% ind.pair.potential$resource.consumer,]$occur <- 2
occ.splist.ind.potential$occur <- factor(occ.splist.ind.potential$occur, levels = c(0,1,2))

## 1. predation matrix ------
occur.ind.potential <- ggplot(occ.splist.ind.potential, aes(resource, consumer)) +  
  geom_point(aes(fill=occur), pch = 21, size = 3) + 
  scale_fill_manual(values = c("White", "Black", "Grey"), name="Predation") +
  coord_flip() + 
  theme_minimal()  +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                            axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) + 
  labs(x= "", y = "") 


## 2. body mass  ------
col.body <- c("yellow", "orange", "red", "darkred")
mass.ind <- rbind(data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Indomalayan",]$resource,
                             bodymass = region.r.sp.trait[region.r.sp.trait$Region == "Indomalayan",]$adult_mass_g_log,
                             Type = "body"),
                  data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Indomalayan",]$resource,
                             bodymass = region.r.sp.trait[region.r.sp.trait$Region == "Indomalayan",]$brain_mass_g_log,
                             Type = "brain"))
mass.ind$resource.ordered <- factor(mass.ind$resource, levels = sort(unique(mass.ind$resource),decreasing = T))
mass.ind.plot <- ggplot(mass.ind, aes(resource.ordered, 1, fill=bodymass)) +  
  facet_wrap(~Type) + 
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  geom_tile(color="white") +
  scale_fill_gradientn(colours = col.body,name = "log10-mass (g)") +
  coord_flip() +  theme_minimal() +theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")
## 3. body length  ------
col.length <- c("pink","#ffb6db","#ff6db6", "violetred")
length.ind <- rbind(data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Indomalayan",]$resource,
                               length = region.r.sp.trait[region.r.sp.trait$Region == "Indomalayan",]$adult_body_length_mm_log))
length.ind$resource.ordered <- factor(length.ind$resource, levels = sort(unique(length.ind$resource),decreasing = T))
length.ind.plot <- ggplot(length.ind, aes(resource.ordered, 1, fill=length)) +  
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  geom_tile(color="white") +
  scale_fill_gradientn(colours = col.length,name = "log10-body\nlength (mm)") +
  coord_flip() +  theme_minimal() +theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")
## 4. dispersal distance ------
col.dist <- c("lavender","#b66dff","mediumpurple4", "#490092")
dist.ind <- rbind(data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Indomalayan",]$resource,
                             dist = region.r.sp.trait[region.r.sp.trait$Region == "Indomalayan",]$dispersal_km_log))
dist.ind$resource.ordered <- factor(dist.ind$resource, levels = sort(unique(dist.ind$resource),decreasing = T))
dist.ind.plot <- ggplot(dist.ind, aes(resource.ordered, 1, fill=dist)) +  
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  geom_tile(color="white") +
  scale_fill_gradientn(colours = col.dist,name = "log10-dispersal\ndistance (km)") +
  coord_flip() +  theme_minimal() +theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")

## 5. dietary and habitat breadth ------
breadth.ind <- rbind(data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Indomalayan",]$resource,
                                Breadth = region.r.sp.trait[region.r.sp.trait$Region == "Indomalayan",]$habitat_breadth_n,
                                Type = "habitat"),
                     data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Indomalayan",]$resource,
                                Breadth = region.r.sp.trait[region.r.sp.trait$Region == "Indomalayan",]$det_diet_breadth_n, 
                                Type = "diet"))
breadth.ind$resource.ordered <- factor(breadth.ind$resource, levels = sort(unique(breadth.ind$resource),decreasing = T))
col <- c("#BCE4D8","#7BBFC9", "#52A8BC", "#52A8BC", "#3690AE", "#3480A2","#317197","#2E628C","#2C5985")
breadth.ind.plot <- ggplot(breadth.ind, aes(resource.ordered, 1, fill=Breadth)) +  
  facet_wrap(~Type) +
  geom_tile(color="white") +
  scale_fill_gradientn(colours = col, name = "Breadth") +
  coord_flip() +  theme_minimal() +theme(strip.text = element_blank(),
                                         panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")

## 6. hibernation and fossor  -----
col.bin <- c("grey80", "grey20")
hib.fossor.ind <- rbind(data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Indomalayan",]$resource,
                                   bin = region.r.sp.trait[region.r.sp.trait$Region == "Indomalayan",]$hibernation_torpor,
                                   Type = "hiber"),
                        data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Indomalayan",]$resource,
                                   bin = region.r.sp.trait[region.r.sp.trait$Region == "Indomalayan",]$fossoriality,
                                   Type = "fossor"))
hib.fossor.ind$resource.ordered <- factor(hib.fossor.ind$resource, levels = sort(unique(hib.fossor.ind$resource),decreasing = T))
hib.fossor.ind$bin <- factor(hib.fossor.ind$bin, levels = c(0,1))
bin.plot <- ggplot(hib.fossor.ind, aes(resource.ordered, 1, fill=bin)) +  
  facet_wrap(~Type) +
  geom_tile(color="white") +
  scale_fill_manual(values = c("grey80", "grey20"), name = "Binary") +
  coord_flip() +  theme_minimal() +theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")
## 7. activity pattern -----
activity.ind <- data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Indomalayan",]$resource,
                           activity_cycle = region.r.sp.trait[region.r.sp.trait$Region == "Indomalayan",]$activity_cycle)
activity.ind$resource.ordered <- factor(activity.ind$resource, levels = sort(activity.ind$resource,decreasing = T))
activity.ind.plot <- ggplot(activity.ind, aes(resource.ordered, 1, fill=activity_cycle)) +  
  geom_tile(color="white") +
  scale_fill_manual(values = c("dodgerblue","wheat2", "goldenrod"), name = "Activity pattern") +
  coord_flip() +  theme_minimal() +theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")
## 8. forage strata -----
forage.ind <- data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Indomalayan",]$resource,
                         forage = region.r.sp.trait[region.r.sp.trait$Region == "Indomalayan",]$foraging_stratum)
forage.ind$resource.ordered <- factor(forage.ind$resource, levels = sort(forage.ind$resource,decreasing = T))
forage.ind.plot <- ggplot(forage.ind, aes(resource.ordered, 1, fill= forage)) +  
  geom_tile(color="white") +
  scale_fill_manual(values = c("tan4","olivedrab3", "green4"), name = "Forage stratum") +
  coord_flip() +  theme_minimal() +theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")
## 9. trophic level -----
level.ind <- data.frame(resource = region.r.sp.trait[region.r.sp.trait$Region == "Indomalayan",]$resource,
                        trophic_level = region.r.sp.trait[region.r.sp.trait$Region == "Indomalayan",]$trophic_level)
level.ind$resource.ordered <- factor(level.ind$resource, levels = sort(level.ind$resource,decreasing = T))
level.ind.plot <- ggplot(level.ind, aes(resource, 1, fill=trophic_level)) +  
  geom_tile(color="white") +
  scale_fill_manual(values = c("green4","tan3", "red3"), name = "Trophic level") +
  coord_flip() +  theme_minimal() +theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) + labs(x= "", y = "")



## 10. Export to combine --------
trait.ind <- ggarrange(mass.ind.plot +  theme(legend.position = "none", strip.text = element_blank()), 
                       length.ind.plot + rremove("y.text") +  theme(legend.position = "none"), 
                       dist.ind.plot + rremove("y.text") +  theme(legend.position = "none"), 
                       breadth.ind.plot + rremove("y.text") +  theme(legend.position = "none", strip.text = element_blank()), 
                       bin.plot + rremove("y.text") +  theme(legend.position = "none", strip.text = element_blank()), 
                       activity.ind.plot + rremove("y.text") +  theme(legend.position = "none"), 
                       forage.ind.plot + rremove("y.text") +  theme(legend.position = "none"), 
                       level.ind.plot + rremove("y.text") +  theme(legend.position = "none"), 
                       nrow = 1, widths = c(.6,.18,.18,.3,.3,.18,.18,.18)) + bgcolor("white") + border("white")

trait.ind.legend <- ggarrange(mass.ind.plot  , 
                              length.ind.plot + rremove("y.text") , 
                              dist.ind.plot + rremove("y.text") , 
                              breadth.ind.plot + rremove("y.text") , 
                              bin.plot + rremove("y.text"), 
                              nrow = 1, widths = c(.6,.18,.18,.3,.3)) + bgcolor("white") + border("white")

trait.ind.legend <- ggarrange(mass.ind.plot  , 
                              activity.ind.plot + rremove("y.text") , 
                              forage.ind.plot + rremove("y.text") , 
                              level.ind.plot + rremove("y.text"), 
                              nrow = 1, widths = c(.6,.18,.18,.18)) + bgcolor("white") + border("white")




