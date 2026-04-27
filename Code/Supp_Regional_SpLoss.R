#
# R script to summarize regional species loss
#
#  Input 
#   1. Rowan_Comm_10kmBuff_Env_cov.csv
#   2. Foodweb_metrics.csv
#   3. PHYLACINE_LQ_sp_gain_loss_no_ocean_no_homo_410community.RData
#   
#  Trait values PHYLACINE dataset v1.2.1  (Faurby et al. 2018 Ecology) 
#      https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2443#support-information-section
#     data repository: https://megapast2future.github.io/PHYLACINE_1.2/
#   1. Trait_data.csv
#
#  Supplementary Tables S3 - S4
#  Figure S5
#


food.web.metrics <- read.csv( "Foodweb_metrics.csv")
buff.10km.env.cov <- read.csv("Rowan_Comm_10kmBuff_Env_cov.csv")

load("PHYLACINE_LQ_sp_gain_loss_no_ocean_no_homo_410community.RData")
food.web.metrics.anaylzed <- buff.10km.env.cov[buff.10km.env.cov$Community_unify %in% food.web.metrics$Community,]

#Call PHYLACINE trait to set guild information for present natural communities 
mamm.func <- read.csv("/Users/chiahsieh/Documents/Research/Data/Trait/Mammal_Trait/PHYLACINE/Traits/Trait_data.csv")
mamm.func <- mamm.func[!duplicated(mamm.func$Binomial.1.2),]
mamm.func$Mass.kg <- round(mamm.func$Mass.g/1000,2)
mamm.func$Mass.g.log <- log10(mamm.func$Mass.g)
mamm.func$Guild <- ""
mamm.func[mamm.func$Diet.Plant > 50,]$Guild <- "Herbivore"
mamm.func[mamm.func$Diet.Invertebrate > 50,]$Guild <- "Insectivore"
mamm.func[mamm.func$Diet.Vertebrate > 50,]$Guild <- "Carnivore"
mamm.func[mamm.func$Guild =="",]$Guild <- "Omnivore"
mamm.func <- mamm.func[,c("Binomial.1.2",  "Mass.g.log","Mass.kg","Guild","Family.1.2", "IUCN.Status.1.2")]


natural.present.realm.sp.pool <- left_join(nature.present.comm.SpList[nature.present.comm.SpList$Community %in% food.web.metrics.anaylzed$Community_unify,], food.web.metrics.anaylzed[,c("Community_unify", "Region")], join_by(Community == Community_unify))
natural.present.realm.sp.pool.trait <- left_join(unique(natural.present.realm.sp.pool[,c( "Species","Region")]), mamm.func, join_by(Species == Binomial.1.2))

current.realm.sp.pool <- left_join(current.comm.SpList[current.comm.SpList$Community %in% food.web.metrics.anaylzed$Community_unify,], food.web.metrics.anaylzed[,c("Community_unify", "Region")], join_by(Community == Community_unify))
current.realm.sp.pool.trait <- left_join(unique(current.realm.sp.pool[,c("Species","Region")]), mamm.func, join_by(Species == Binomial.1.2))

LQ.loss.realm.sp.pool <- left_join(LQ.loss.gain.comm.SpList[LQ.loss.gain.comm.SpList$Community %in% food.web.metrics.anaylzed$Community_unify,], food.web.metrics.anaylzed[,c("Community_unify", "Region")], join_by(Community == Community_unify))
LQ.loss.realm.sp.pool.trait <- left_join(unique(LQ.loss.realm.sp.pool[LQ.loss.realm.sp.pool$LossGain == "loss", c("Species","Region")]), mamm.func, join_by(Species == Binomial.1.2))

LQ.loss.realm.sp.pool.trait.extinction <- LQ.loss.realm.sp.pool.trait[LQ.loss.realm.sp.pool.trait$IUCN.Status.1.2 %in% c("EP", "EW", "EX"),]                              
LQ.loss.realm.sp.pool.trait.range.loss <- LQ.loss.realm.sp.pool.trait[!LQ.loss.realm.sp.pool.trait$IUCN.Status.1.2 %in% c("EP", "EW", "EX"),]                                  


realm.species.pool <- rbind(data.frame(natural.present.realm.sp.pool.trait, Status = "nature.present"),
                            data.frame(current.realm.sp.pool.trait, Status = "current"),
                            data.frame(LQ.loss.realm.sp.pool.trait.extinction, Status = "extinction"),
                            data.frame(LQ.loss.realm.sp.pool.trait.range.loss, Status = "range.contraction"))

#Summary table for Table S3-----------------------
### Number of unique species lost from at least one local community due to extinction or range contractions relative to the total number of species in the natural-present realm species pool
table(realm.species.pool$Region, realm.species.pool$Status,realm.species.pool$Guild)

#Summary table for Table S4-----------------------
### Number of communities having carnivore species lost in local species pools by dietary guilds 
LQ.loss.gain.comm.SpList <- left_join(LQ.loss.gain.comm.SpList, food.web.metrics.anaylzed[,c("Community_unify", "Region")], join_by(Community == Community_unify))
LQ.loss.comm.SpList <- left_join(LQ.loss.gain.comm.SpList[LQ.loss.gain.comm.SpList$LossGain =="loss",], mamm.func, join_by(Species == Binomial.1.2))
LQ.loss.comm.SpList[is.na(LQ.loss.comm.SpList$Mass.g) == T,]
LQ.loss.carni.SpList <- LQ.loss.comm.SpList[LQ.loss.comm.SpList$Guild == "Carnivore",] %>% group_by(Region, Species, Family.1.2, IUCN.Status.1.2, Mass.g.log) %>% summarise(N.comm = length(Community))
#LQ.loss.herb.SpList <- LQ.loss.comm.SpList[LQ.loss.comm.SpList$Guild == "Herbivore",] %>% group_by(Region, Species, Family.1.2, IUCN.Status.1.2, Mass.g.log) %>% summarise(N.comm = length(Community))
#LQ.loss.insect.SpList <- LQ.loss.comm.SpList[LQ.loss.comm.SpList$Guild == "Insectivore",] %>% group_by(Region, Species, Family.1.2, IUCN.Status.1.2, Mass.g.log) %>% summarise(N.comm = length(Community))
#LQ.loss.omni.SpList <- LQ.loss.comm.SpList[LQ.loss.comm.SpList$Guild == "Omnivore",] %>% group_by(Region, Species, Family.1.2, IUCN.Status.1.2, Mass.g.log) %>% summarise(N.comm = length(Community))


#Supplementary Figure S5 -----------------------
### Realm variation in species body mass for each dietary guild in realm species pools
realm.species.pool$Status <- factor(realm.species.pool$Status, levels = c("nature.present",  "current", "extinction", "range.contraction"), labels = c("Present\nnatural   ", "Current  ","Extinction ", "Range\ncontraction  "))
realm.species.pool$Guild <- factor(realm.species.pool$Guild, levels = c("Herbivore", "Carnivore", "Insectivore", "Omnivore"))
realm.species.pool$Region <- factor(realm.species.pool$Region, levels = c("Neotropic", "Afrotropic", "Indomalayan"))

status.color <- c("Present\nnatural   " = "grey40", "Current  " = "grey60", "Extinction " = "brown1", "Range\ncontraction  " = "pink")

ggplot(realm.species.pool, aes(y = Mass.g.log, x = Region, fill = Status)) + 
  geom_boxplot(alpha = 1, width=1.5, varwidth = T, color = "grey20", linewidth = 0.3) + 
  facet_wrap(~Guild, ncol = 2) +
  ylim(c(2.5,7.8)) +
  scale_fill_manual(values = status.color) + theme_bw() +
  theme_bw() +  theme(axis.text = element_text(size = 16),
                      axis.title = element_text(size = 18),
                      legend.title = element_text(size = 0),
                      legend.text = element_text(size = 12),
                      strip.text = element_text(size = 18, face = "bold"),
                      strip.background = element_blank(),
                      legend.background = element_rect(color = "grey40", fill = NA, linewidth = 0.2),
                      legend.position = c(0.752, 0.93),
                      legend.justification = "center",
                      panel.grid = element_blank()) +
  guides(fill = guide_legend(ncol = 4)) +
  labs(x = "", y = expression(log[10]*"(body mass [g])"), fill = "", color = "")



