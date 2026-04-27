#
# Estimate lineage history of present natural communities and contamporary food web lineage structure 
#
#
# Input
#   1. Community408_carni_web_doc_inf95_500g_SR3_commcollect.Rdata
#   2. PHYLACINE_LQ_sp_gain_loss_no_ocean_no_homo_410community.RData
#   3. Rowan_Comm_10kmBuff_Env_cov.csv
#
# Source function: DRstat from Title and Rabosky 2019 Methods in Ecology and Evolution 
#   https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13153
# 
# Phylogeny and Taxonomy from PHYLACINE
#   https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2443#support-information-section
#   data repository: https://megapast2future.github.io/PHYLACINE_1.2/
#    1. Complete_phylogeny.nex
#    2. Synonymy_table_valid_species_only.csv
#
#  Trait values PHYLACINE dataset v1.2.1  (Faurby et al. 2018 Ecology) 
#      https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2443#support-information-section
#     data repository: https://megapast2future.github.io/PHYLACINE_1.2/
#    1. Trait_data.csv
#
#  
#  Output 
#    1. Foodweb_lineage_history.csv
#
#  Supplementary Figures S4 and S16
#

setwd(".")

# I. Estimated species-specific lineage hisotry using DR metric -----------------
## Source function: DRstat from Title and Rabosky 2019 Methods in Ecology and Evolution 
## https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13153
source("./tipRateFunctions.R")
##phylogeny from PHYLACINE 1.2 (data repository: https://megapast2future.github.io/PHYLACINE_1.2/)
mammal.tree <- read.nexus("./Complete_phylogeny.nex") 

set.seed(123)
sample.tree.n <- 100
sample.tree <- sample(1:1000, size =  sample.tree.n, replace = F)

Species.tipDR <-c()
for (z in 1:sample.tree.n){
  print(z)
  #generate species' tip speciation rate with full tree (N = 5831 sp)
  DR.whole <- DRstat(mammal.tree[[sample.tree[z]]])
  DR.whole.M <- data.frame(Species = names(DR.whole), DR = DR.whole,
                           sample.tree = rep(sample.tree[z], length(DR.whole)))
  DR.whole.M$Species <- gsub("_", " ",DR.whole.M$Species)
  
  Species.tipDR <- rbind(Species.tipDR, DR.whole.M)
  
}

##average over 1000 phylogeny trees 
DR.whole.M.mean <- Species.tipDR %>% group_by(Species) %>% summarise(mean.sp.DR = mean(DR))
DR.whole.M.mean$Species <- sub(" ", "_", DR.whole.M.mean$Species)

#II. Food-web lineage structure ---------
load("Data/Community408_carni_web_doc_inf95_500g_SR3_commcollect.Rdata")
#collect food web species list 
comm.trophic.link.inf95.10p <- c()
for (i in 1:length(Tropic.carni.webs.inf95.10p)) {
  #Collect all interactions in a local food web  
  comm.trophic.link.inf95.10p <- rbind(comm.trophic.link.inf95.10p, 
                                       data.frame(Tropic.carni.webs.inf95.10p[[i]][3],
                                                  Community = Tropic.carni.webs.inf95.10p[[i]][2]$properties$Community_unify,
                                                  Region = Tropic.carni.webs.inf95.10p[[i]][2]$properties$Region))
}
comm.trophic.link.inf95.10p.splist <- comm.trophic.link.inf95.10p %>% group_by(Community) %>% reframe(Species = unique(c(trophic.links.resource, trophic.links.consumer)))
comm.trophic.link.inf95.10p.splist$Species <- sub(" ", "_", comm.trophic.link.inf95.10p.splist$Species)
comm.trophic.link.inf95.10p.splist.DR <- left_join(comm.trophic.link.inf95.10p.splist, DR.whole.M.mean, join_by(Species == Species))

foodweb.mean.DR <- comm.trophic.link.inf95.10p.splist.DR %>% group_by(Community) %>% summarise(foodweb_meanDR = round(mean(mean.sp.DR),2))


#III. Present-natural community ---------
load("Data/PHYLACINE_LQ_sp_gain_loss_no_ocean_no_homo_410community.RData")
studied.nature.present.comm.SpList <- nature.present.comm.SpList[nature.present.comm.SpList$Community %in% comm.trophic.link.inf95.10p$Community,]
length(unique(studied.nature.present.comm.SpList$Community))

studied.nature.present.comm.SpList.DR <- left_join(studied.nature.present.comm.SpList, DR.whole.M.mean, join_by(Species== Species))
natural.present.comm.mean.DR <- studied.nature.present.comm.SpList.DR %>% group_by(Community) %>% summarise(naturalpresent_comm_meanDR = round(mean(mean.sp.DR),2))

#IV. Combine ----
food.web.evo.history <- left_join(natural.present.comm.mean.DR, foodweb.mean.DR, join_by(Community == Community))
write.csv(food.web.evo.history, "Data/Foodweb_lineage_history.csv", row.names = F)



#Figure S4  ----------------------
### Phylogeny of studied species 
#library(devtools)
#devtools::install_github("YuLab-SMU/ggtree")

library(ape); library(ggtree); library(treeio); library(ggnewscale)
library(ggplot2) 
# Get mammal taxonomy  from  PHYLACINE 1.2 (data repository: https://megapast2future.github.io/PHYLACINE_1.2/)
phyla.taxon <- read.csv("./Synonymy_table_valid_species_only.csv")
region.metaweb.intx <- comm.trophic.link.inf95.10p

region.metaweb.intx$trophic.links.resource <- gsub(" ","_",region.metaweb.intx$trophic.links.resource)
region.metaweb.intx$trophic.links.consumer <- gsub(" ","_",region.metaweb.intx$trophic.links.consumer)

carni.predator.sp <- data.frame(consumer = unique(region.metaweb.intx$trophic.links.consumer))
carni.predator.sp <- left_join(carni.predator.sp, phyla.taxon[,c("Binomial.1.2", "Family.1.2")], join_by(consumer==Binomial.1.2))
predator.family <- data.frame(carni.predator.sp[,c("Family.1.2")])
rownames(predator.family) <- carni.predator.sp$consumer; colnames(predator.family) <- "family"
length(unique(rownames(predator.family))); length(unique(predator.family$family))

prey.sp <- data.frame(resource = unique(region.metaweb.intx$trophic.links.resource))
prey.sp <- left_join(prey.sp, phyla.taxon[,c("Binomial.1.2", "Order.1.2")], join_by(resource==Binomial.1.2))
prey.order <- data.frame(prey.sp[,c("Order.1.2")])
rownames(prey.order) <- prey.sp$resource; colnames(prey.order) <- "order"
length(unique(rownames(prey.order))); length(unique(prey.order$order))


predator.fam.colr <- c("Canidae" = "#000000","Felidae"="grey40", "Herpestidae"="grey60", "Hyaenidae"="grey80",
                       "Mustelidae" ="lightsteelblue4","Nandiniidae" ="lightsteelblue3","Procyonidae"="dodgerblue4",
                       "Ursidae"="lightblue","Viverridae"="dodgerblue3")

prey.order.colr <- c("Carnivora"="#000000","Cetartiodactyla"="#009292","Cingulata"="#004949","Didelphimorphia"="#ff6db6","Hyracoidea"="#ffb6db",
                     "Lagomorpha"="#b66dff",  "Perissodactyla"="#6db6ff","Pholidota"="#006ddb","Pilosa"="#490092","Primates"="#b6dbff",
                     "Proboscidea"="#920000","Rodentia"="lightsteelblue4","Tubulidentata"="#924900")

###plot tree
predator.tree <- ape::drop.tip(mammal.tree[[1]], mammal.tree[[1]]$tip.label[!mammal.tree[[1]]$tip.label %in% rownames(predator.family)])
prey.tree <- ape::drop.tip(mammal.tree[[1]], mammal.tree[[1]]$tip.label[!mammal.tree[[1]]$tip.label %in% rownames(prey.order)])

#### All predators and prey
predator.circ <- ggtree(predator.tree, layout='circular') 
pred.all <- gheatmap(predator.circ, predator.family, offset=.8, width=.1,
                     colnames_angle=95, colnames = F,color=NULL) +
  scale_fill_manual(values = predator.fam.colr, name="Family") 

prey.circ <- ggtree(prey.tree, layout='circular') 
prey.all <- gheatmap(prey.circ, prey.order, offset=.8, width=.1,
                     colnames_angle=95, colnames = F,color=NULL) +
  scale_fill_manual(values = prey.order.colr, name="Order") 

#### Neotropical predators and prey
carni.predator.sp.neo <- data.frame(consumer = unique(region.metaweb.intx[region.metaweb.intx$Region == "Neotropic",]$trophic.links.consumer))
carni.predator.sp.neo <- left_join(carni.predator.sp.neo, phyla.taxon[,c("Binomial.1.2", "Family.1.2")], join_by(consumer==Binomial.1.2))
predator.family.neo <- data.frame(carni.predator.sp.neo[,c("Family.1.2")])
rownames(predator.family.neo) <- carni.predator.sp.neo$consumer; colnames(predator.family.neo) <- "family"
length(unique(rownames(predator.family.neo))); length(unique(predator.family.neo$family))

prey.sp.neo <- data.frame(resource = unique(region.metaweb.intx[region.metaweb.intx$Region == "Neotropic",]$trophic.links.resource))
prey.sp.neo <- left_join(prey.sp.neo, phyla.taxon[,c("Binomial.1.2", "Order.1.2")], join_by(resource==Binomial.1.2))
prey.order.neo <- data.frame(prey.sp.neo[,c("Order.1.2")])
rownames(prey.order.neo) <- prey.sp.neo$resource; colnames(prey.order.neo) <- "order"
length(unique(rownames(prey.order.neo))); length(unique(prey.order.neo$order))

predator.neo.tree <- ape::drop.tip(mammal.tree[[1]], mammal.tree[[1]]$tip.label[!mammal.tree[[1]]$tip.label %in% rownames(predator.family.neo)])
prey.neo.tree <- ape::drop.tip(mammal.tree[[1]], mammal.tree[[1]]$tip.label[!mammal.tree[[1]]$tip.label %in% rownames(prey.order.neo)])

predator.neo.circ <- ggtree(predator.neo.tree, layout='circular') 
pred.neo <- gheatmap(predator.neo.circ, predator.family, offset=.8, width=.1,
                     colnames_angle=95, colnames = F,color=NULL) +
  scale_fill_manual(values = predator.fam.colr, name="Family")

prey.neo.circ <- ggtree(prey.neo.tree, layout='circular') 
prey.neo <- gheatmap(prey.neo.circ, prey.order, offset=.8, width=.1,
                     colnames_angle=95, colnames = F,color=NULL) +
  scale_fill_manual(values = prey.order.colr, name="Order") 

#### Afrotropical predators and prey
carni.predator.sp.afro <- data.frame(consumer = unique(region.metaweb.intx[region.metaweb.intx$Region == "Afrotropic",]$trophic.links.consumer))
carni.predator.sp.afro <- left_join(carni.predator.sp.afro, phyla.taxon[,c("Binomial.1.2", "Family.1.2")], join_by(consumer==Binomial.1.2))
predator.family.afro <- data.frame(carni.predator.sp.afro[,c("Family.1.2")])
rownames(predator.family.afro) <- carni.predator.sp.afro$consumer; colnames(predator.family.afro) <- "family"
length(unique(rownames(predator.family.afro))); length(unique(predator.family.afro$family))

prey.sp.afro <- data.frame(resource = unique(region.metaweb.intx[region.metaweb.intx$Region == "Afrotropic",]$trophic.links.resource))
prey.sp.afro <- left_join(prey.sp.afro, phyla.taxon[,c("Binomial.1.2", "Order.1.2")], join_by(resource==Binomial.1.2))
prey.order.afro <- data.frame(prey.sp.afro[,c("Order.1.2")])
rownames(prey.order.afro) <- prey.sp.afro$resource; colnames(prey.order.afro) <- "order"
length(unique(rownames(prey.order.afro))); length(unique(prey.order.afro$order))

predator.afro.tree <- ape::drop.tip(mammal.tree[[1]], mammal.tree[[1]]$tip.label[!mammal.tree[[1]]$tip.label %in% rownames(predator.family.afro)])
prey.afro.tree <- ape::drop.tip(mammal.tree[[1]], mammal.tree[[1]]$tip.label[!mammal.tree[[1]]$tip.label %in% rownames(prey.order.afro)])

predator.afro.circ <- ggtree(predator.afro.tree, layout='circular') 
pred.afro <- gheatmap(predator.afro.circ, predator.family, offset=.8, width=.1,
                      colnames_angle=95, colnames = F,color=NULL) +
  scale_fill_manual(values = predator.fam.colr, name="Family") 

prey.afro.circ <- ggtree(prey.afro.tree, layout='circular') 
prey.afro <- gheatmap(prey.afro.circ, prey.order, offset=.8, width=.1,
                      colnames_angle=95, colnames = F,color=NULL) +
  scale_fill_manual(values = prey.order.colr, name="Order") 

#### Indomalayan predators and prey
carni.predator.sp.ind <- data.frame(consumer = unique(region.metaweb.intx[region.metaweb.intx$Region == "Indomalayan",]$trophic.links.consumer))
carni.predator.sp.ind <- left_join(carni.predator.sp.ind, phyla.taxon[,c("Binomial.1.2", "Family.1.2")], join_by(consumer==Binomial.1.2))
predator.family.ind <- data.frame(carni.predator.sp.ind[,c("Family.1.2")])
rownames(predator.family.ind) <- carni.predator.sp.ind$consumer; colnames(predator.family.ind) <- "family"
length(unique(rownames(predator.family.ind))); length(unique(predator.family.ind$family))

prey.sp.ind <- data.frame(resource = unique(region.metaweb.intx[region.metaweb.intx$Region == "Indomalayan",]$trophic.links.resource))
prey.sp.ind <- left_join(prey.sp.ind, phyla.taxon[,c("Binomial.1.2", "Order.1.2")], join_by(resource==Binomial.1.2))
prey.order.ind <- data.frame(prey.sp.ind[,c("Order.1.2")])
rownames(prey.order.ind) <- prey.sp.ind$resource; colnames(prey.order.ind) <- "order"
length(unique(rownames(prey.order.ind))); length(unique(prey.order.ind$order))

predator.ind.tree <- ape::drop.tip(mammal.tree[[1]], mammal.tree[[1]]$tip.label[!mammal.tree[[1]]$tip.label %in% rownames(predator.family.ind)])
prey.ind.tree <- ape::drop.tip(mammal.tree[[1]], mammal.tree[[1]]$tip.label[!mammal.tree[[1]]$tip.label %in% rownames(prey.order.ind)])

predator.ind.circ <- ggtree(predator.ind.tree, layout='circular') 
pred.ind <- gheatmap(predator.ind.circ, predator.family, offset=.8, width=.1,
                     colnames_angle=95, colnames = F,color=NULL) +
  scale_fill_manual(values = predator.fam.colr, name="Family") 

prey.ind.circ <- ggtree(prey.ind.tree, layout='circular') 
prey.ind <- gheatmap(prey.ind.circ, prey.order, offset=.8, width=.1,
                     colnames_angle=95, colnames = F,color=NULL) +
  scale_fill_manual(values = prey.order.colr, name="Order") 

ggarrange(pred.all, prey.all, pred.neo, prey.neo, pred.afro, prey.afro, pred.ind, prey.ind,
          labels = c("A", "B", "C", "D", "E", "F","G", "H"), vjust = 4, ncol = 2, nrow = 4, align = "hv")



#Figure S16 ----------------
### Associations of body mass and species-specific DR metric by species 
#Mass values from  PHYLACINE 1.2 (data repository: https://megapast2future.github.io/PHYLACINE_1.2/)
mamm.func <- read.csv("./Trait_data.csv")
mamm.func <- mamm.func[!duplicated(mamm.func$Binomial.1.2),]
mamm.func$Mass.g.log <- log10(mamm.func$Mass.g)

mamm.func$Guild <- ""
mamm.func[mamm.func$Diet.Plant > 50,]$Guild <- "Herbivore"
mamm.func[mamm.func$Diet.Invertebrate > 50,]$Guild <- "Insectivore"
mamm.func[mamm.func$Diet.Vertebrate > 50,]$Guild <- "Carnivore"
mamm.func[mamm.func$Guild =="",]$Guild <- "Omnivore"

mamm.func <- mamm.func[,c("Binomial.1.2","Mass.g", "Mass.g.log","Guild","Family.1.2", "IUCN.Status.1.2")]

#Join region information 
buff.10km.env.cov <- read.csv("Data/Rowan_Comm_10kmBuff_Env_cov.csv")
buff.10km.env.cov <- buff.10km.env.cov[buff.10km.env.cov$EVI_month.mean >0,]

studied.nature.present.comm.SpList <- left_join(studied.nature.present.comm.SpList, buff.10km.env.cov[,c("Community_unify", "Region")], join_by(Community == Community_unify))
studied.nature.present.region.SpList <- unique(studied.nature.present.comm.SpList[,c( "Species" , "Region")])
studied.nature.present.region.SpList.DR <- left_join(studied.nature.present.region.SpList, mamm.func, join_by(Species == Binomial.1.2))
studied.nature.present.region.SpList.DR <- left_join(studied.nature.present.region.SpList.DR, DR.whole.M.mean, join_by(Species == Species))

Region.color <- c("Neotropic" = "#adc9bc","Afrotropic" = "#ffb453" , "Indomalaya" = "#e6aaf0")
dodge <- position_dodge(width = 10)

studied.nature.present.region.SpList.DR$Region <- factor(studied.nature.present.region.SpList.DR$Region, levels = c("Neotropic", "Afrotropic", "Indomalayan"), labels = c("Neotropic", "Afrotropic", "Indomalaya"))

c.bm.DR <- ggplot(studied.nature.present.region.SpList.DR[studied.nature.present.region.SpList.DR$Guild == "Carnivore",], aes(Mass.g.log, mean.sp.DR, color = Region, fill = Region)) + 
  xlim(2.6, 5.5) +
  geom_point(alpha = 0.8) + geom_smooth(method = "lm", alpha = 0.2) + 
  scale_color_manual(values = Region.color) + scale_fill_manual(values = Region.color) +
  theme_bw() +  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14),
                      legend.title = element_text(size = 16),legend.text = element_text(size = 14),
                      panel.grid = element_blank(),legend.position  ="none") +
  labs(x = expression('log' [10] * '(mean body mass [g])'), y = "Species-specific lineage\ndiversification history: carnivore", fill = "", color ="")

bm.DR <- ggplot(studied.nature.present.region.SpList.DR, aes(Mass.g.log, mean.sp.DR, color = Region, fill = Region)) + 
  xlim(2.6, 7.2) +  ylim(0,2) +
  geom_point(alpha = 0.8) + geom_smooth(method = "lm", alpha = 0.2) + 
  scale_color_manual(values = Region.color) + scale_fill_manual(values = Region.color) +
  theme_bw() +  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14),
                      legend.title = element_text(size = 16),legend.text = element_text(size = 14),
                      panel.grid = element_blank(),legend.position  ="none") +
  labs(x = expression('log' [10] * '(mean body mass [g])'), y = "Species-specific lineage\ndiversification history: all guilds", fill = "", color ="")

ggarrange(c.bm.DR, bm.DR, labels = c("A", "B"), ncol = 2, nrow = 1, common.legend = T)
