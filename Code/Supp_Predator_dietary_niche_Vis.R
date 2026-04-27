#
#  Viusalize predator dietary niche spece 
#
#  Input 
#     1. Gower_11traits_PCoA_Axis13_inf95_10p_500g.csv
#     2. Predator_prey_pairs_for_dietary_niche_space.csv
#
#  Trait and Taxonomy from Soria et al. (2021) COMBINE dataset https://doi.org/10.1002/ecy.3344 to have fine-scale trait information to characterize ecological similarity
#     1. trait_data_imputed.csv
#


library(ggplot2); library(ggpubr); library(dplyr)
traits <- read.csv("Gower_11traits_PCoA_Axis13_inf95_10p_500g.csv")
meta.web <- read.csv("Predator_prey_pairs_for_dietary_niche_space.csv")
mamm.func <- read.csv("/Users/chiahsieh/Documents/Research/Data/Trait/Mammal_Trait/COMBINE/trait_data_imputed.csv")
mamm.func <- mamm.func[!duplicated(mamm.func$phylacine_binomial),]
mamm.func$adult_kg <- round(mamm.func$adult_mass_g/1000,2)
mamm.func$adult_g.log <- round(log10(mamm.func$adult_mass_g),2)

mamm.func.c <- mamm.func[,c("phylacine_binomial", "adult_kg", "adult_g.log")]
colnames(mamm.func.c) <- c("phylacine_binomial", "adult_kg.c", "adult_g.log.c")

##Dietary niche spaces of all and regions ======
all.r.sp <- unique(meta.web[,c("Region", "trophic.links.resource")])
all.r.sp <- left_join(all.r.sp, traits[,c(1:13)], join_by(trophic.links.resource == X))
all.r.sp <- left_join(all.r.sp, mamm.func[,c("phylacine_binomial", "adult_kg")], join_by(trophic.links.resource == phylacine_binomial))

all.hull <- all.r.sp[chull(all.r.sp$Axis.1, all.r.sp$Axis.2),] 

All.prey <- ggplot() + 
  geom_polygon(data = all.hull, aes(Axis.1,  Axis.2), fill = "grey80", color = NA, alpha = 0.7, lwd = 0.5) + 
  geom_point(data = all.r.sp, aes(Axis.1,  Axis.2), color = "black",alpha = 0.8) + 
  geom_point(data = all.r.sp[all.r.sp$adult_kg >= 10 &  all.r.sp$adult_kg < 45,], aes(Axis.1,  Axis.2), color = "#00DDF9",alpha = 0.8) + 
  geom_point(data = all.r.sp[all.r.sp$adult_kg >= 45,], aes(Axis.1,  Axis.2), color = "#0087CA",alpha = 0.8) +
  #geom_point(data = all.r.sp[all.r.sp$adult_kg >= 45,], aes(Axis.1,  Axis.2), color = "red",alpha = 0.8) + 
  ylim(-3.5,3) + xlim(-4.8,4.8) + theme_bw() + 
  annotate("text", x = -4, y = 2.8, label = "Global", size = 5) + 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        strip.text =  element_text(size = 14, face = "bold"), strip.background = element_blank(),
        legend.title = element_text(size = 16), legend.text = element_text(size = 14),
        panel.grid = element_blank()) +
  labs(x = "PCoA Axis 1 (40.32% of variation)", y ="PCoA Axis 2 (15.94% of variation)")


all.c.sp <- data.frame(consumer = unique(meta.web[,c("trophic.links.consumer")]))
all.c.sp <- left_join(all.c.sp, mamm.func.c, join_by(consumer == phylacine_binomial))

color.option <- "A"

###1. Neotropic -------------
neo.region.r <- all.r.sp[all.r.sp$Region == "Neotropic" & all.r.sp$trophic.links.resource %in% 
                           meta.web[meta.web$Region == "Neotropic",]$trophic.links.resource,]
neo.predator <- unique(meta.web[meta.web$Region == "Neotropic",]$trophic.links.consumer)

neo.predator.hull <- c()
for (i in 1:length(neo.predator)) {
  predator.prey <- neo.region.r[neo.region.r$trophic.links.resource %in% meta.web[meta.web$Region == "Neotropic" & meta.web$trophic.links.consumer == neo.predator[i],]$trophic.links.resource,]
  neo.predator.hull <- rbind(neo.predator.hull,
                             data.frame(consumer = neo.predator[i],
                                        number.of.resource = length(predator.prey$trophic.links.resource),
                                        predator.prey[chull(predator.prey$Axis.1, predator.prey$Axis.2),]))
}  

neo.predator.hull <- left_join(neo.predator.hull, mamm.func.c, join_by(consumer == phylacine_binomial))
neo.predator.hull$consumer.order <- factor(neo.predator.hull$consumer, levels = rev(unique(neo.predator.hull[rev(order(neo.predator.hull$adult_kg.c)),]$consumer)))

Neo.prey.metaweb.by.predator <- ggplot() + 
  geom_polygon(data = all.hull, aes(Axis.1,  Axis.2), fill = "grey80", color ="NA", alpha = 0.7, lwd = 0.5) + 
  geom_polygon(data = neo.predator.hull, aes(Axis.1,  Axis.2, group = consumer.order, fill = adult_kg.c), color = "NA", alpha = 0.2, lwd = 1) + 
  geom_point(data = neo.region.r, aes(Axis.1,  Axis.2), color = "black",alpha = 0.8) + 
  geom_point(data = neo.region.r[neo.region.r$adult_kg >= 10 &  neo.region.r$adult_kg < 45,], aes(Axis.1,  Axis.2), color ="#00DDF9",alpha = 0.8) + 
  geom_point(data = neo.region.r[neo.region.r$adult_kg >= 45,], aes(Axis.1,  Axis.2), color = "#0087CA",alpha = 0.8) + 
  ylim(-3.5,3) + xlim(-4.8,4.8) + theme_bw() + 
  annotate("text", x = -3.6, y = 2.8, label = "Neotropic", size = 5) + 
  scale_fill_viridis_c(limits = range(all.c.sp$adult_kg.c), option = color.option) + 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        strip.text =  element_text(size = 14, face = "bold"), strip.background = element_blank(),
        legend.title = element_text(size = 10), legend.text = element_text(size = 10),
        panel.grid = element_blank(), legend.position = c(0.12, 0.25), legend.background = element_rect(fill = NA), legend.key.size = unit(0.35, "cm")) +
  labs(x = "PCoA Axis 1 (40.32% of variation)", y ="PCoA Axis 2 (15.94% of variation)", fill = "Predator\nbody\nmass (kg)")

###2. Afrotropic -------------
afr.region.r <- all.r.sp[all.r.sp$Region == "Afrotropic" & all.r.sp$trophic.links.resource %in% 
                           meta.web[meta.web$Region == "Afrotropic",]$trophic.links.resource,]
afr.predator <- unique(meta.web[meta.web$Region == "Afrotropic",]$trophic.links.consumer)
length(afr.predator)

afr.predator.hull <- c()
for (i in 1:length(afr.predator)) {
  predator.prey <- afr.region.r[afr.region.r$trophic.links.resource %in% meta.web[meta.web$Region == "Afrotropic" & meta.web$trophic.links.consumer == afr.predator[i],]$trophic.links.resource,]
  afr.predator.hull <- rbind(afr.predator.hull,
                             data.frame(consumer = afr.predator[i],
                                        number.of.resource = length(predator.prey$trophic.links.resource),
                                        predator.prey[chull(predator.prey$Axis.1, predator.prey$Axis.2),]))
}  

afr.predator.hull <- left_join(afr.predator.hull, mamm.func.c, join_by(consumer == phylacine_binomial))
afr.predator.hull$consumer.order <- factor(afr.predator.hull$consumer, levels = rev(unique(afr.predator.hull[rev(order(afr.predator.hull$adult_kg.c)),]$consumer)))

Afr.prey.metaweb.by.predator <- ggplot() + 
  geom_polygon(data = all.hull, aes(Axis.1,  Axis.2), fill = "grey80", color ="NA", alpha = 0.7, lwd = 0.5) + 
  geom_polygon(data = afr.predator.hull, aes(Axis.1,  Axis.2, group = consumer.order, fill = adult_kg.c ), color = "NA", alpha = 0.2, lwd = 1) + 
  geom_point(data = afr.region.r, aes(Axis.1,  Axis.2), color = "black",alpha = 0.8) + 
  geom_point(data = afr.region.r[afr.region.r$adult_kg >= 10 &  afr.region.r$adult_kg < 45,], aes(Axis.1,  Axis.2), color ="#00DDF9",alpha = 0.8) + 
  geom_point(data = afr.region.r[afr.region.r$adult_kg >= 45,], aes(Axis.1,  Axis.2), color = "#0087CA",alpha = 0.8) + 
  ylim(-3.5,3) + xlim(-4.8,4.8) + theme_bw() + 
  annotate("text", x = -3.6, y = 2.8, label = "Afrotropic", size = 5) + 
  scale_fill_viridis_c(limits = range(all.c.sp$adult_kg.c), option = color.option) + 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        strip.text =  element_text(size = 14, face = "bold"), strip.background = element_blank(),
        legend.title = element_text(size = 16), legend.text = element_text(size = 14),
        panel.grid = element_blank(),  legend.position = "none") +
  labs(x = "PCoA Axis 1 (40.32% of variation)", y ="PCoA Axis 2 (15.94% of variation)", fill = "Predator\nbody\nmass (kg)")


###3. Indomalaya -------------
ind.region.r <- all.r.sp[all.r.sp$Region == "Indomalayan" & all.r.sp$trophic.links.resource %in% 
                           meta.web[meta.web$Region == "Indomalayan",]$trophic.links.resource,]
ind.predator <- unique(meta.web[meta.web$Region == "Indomalayan",]$trophic.links.consumer)
length(ind.predator)

ind.predator.hull <- c()
for (i in 1:length(ind.predator)) {
  predator.prey <- ind.region.r[ind.region.r$trophic.links.resource %in% meta.web[meta.web$Region == "Indomalayan" & meta.web$trophic.links.consumer == ind.predator[i],]$trophic.links.resource,]
  ind.predator.hull <- rbind(ind.predator.hull,
                             data.frame(consumer = ind.predator[i],
                                        number.of.resource = length(predator.prey$trophic.links.resource),
                                        predator.prey[chull(predator.prey$Axis.1, predator.prey$Axis.2),]))
}  

ind.predator.hull <- left_join(ind.predator.hull, mamm.func.c, join_by(consumer == phylacine_binomial))
ind.predator.hull$consumer.order <- factor(ind.predator.hull$consumer, levels = rev(unique(ind.predator.hull[rev(order(ind.predator.hull$adult_kg.c)),]$consumer)))

Ind.prey.metaweb.by.predator <- 
  ggplot() + 
  geom_polygon(data = all.hull, aes(Axis.1,  Axis.2), fill = "grey80", color ="NA", alpha = 0.7, lwd = 0.5) + 
  geom_polygon(data = ind.predator.hull, aes(Axis.1,  Axis.2, group = consumer.order, fill = adult_kg.c ), color = "NA", alpha = 0.2, lwd = 1) + 
  geom_point(data = ind.region.r, aes(Axis.1,  Axis.2), color = "black",alpha = 0.8) + 
  geom_point(data = ind.region.r[ind.region.r$adult_kg >= 10 &  ind.region.r$adult_kg < 45,], aes(Axis.1,  Axis.2), color ="#00DDF9",alpha = 0.8) + 
  geom_point(data = ind.region.r[ind.region.r$adult_kg >= 45,], aes(Axis.1,  Axis.2), color = "#0087CA",alpha = 0.8) + 
  ylim(-3.5,3) + xlim(-4.8,4.8) + theme_bw() + 
  annotate("text", x = -3.3, y = 2.8, label = "Indomalaya", size = 5) + 
  scale_fill_viridis_c(limits = range(all.c.sp$adult_kg.c), option = color.option) + 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        strip.text =  element_text(size = 14, face = "bold"), strip.background = element_blank(),
        legend.title = element_text(size = 16), legend.text = element_text(size = 14),
        panel.grid = element_blank(),  legend.position = "none") +
  labs(x = "PCoA Axis 1 (40.32% of variation)", y ="PCoA Axis 2 (15.94% of variation)", fill = "Predator\nbody\nmass (kg)")


ggarrange(All.prey+ rremove("xlab")+ rremove("ylab"), 
          Neo.prey.metaweb.by.predator + rremove("xlab") + rremove("ylab"),  
          Afr.prey.metaweb.by.predator+ rremove("xlab")+ rremove("ylab"), 
          Ind.prey.metaweb.by.predator+ rremove("xlab")+ rremove("ylab"), 
          ncol = 2, nrow = 2, labels = c(" ", " ", " ", " "), heights = c(1, 1, 1,1), widths =  c(1, 1, 1,1), align = "hv")
