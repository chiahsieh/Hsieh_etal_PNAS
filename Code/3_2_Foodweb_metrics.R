#
# Measure food web metrics: predator-basal prey ratio, mean body mass of predator and prey, normalized generality and normalized vulnerability, and mean predator dietary niche breadth and overlap
#
# Input
#  1. Community408_carni_web_doc_inf95_500g_SR3_commcollect.Rdata
#  2. Region_sp_intx_500g.RData
#  3. Regional_carni_links_500g_cosine_similarity_in95_10p_all.csv
#
#  Trait and Taxonomy from Soria et al. (2021) COMBINE dataset https://doi.org/10.1002/ecy.3344 to have fine-scale trait information to characterize ecological similarity
#     1. trait_data_imputed.csv
#
# Output 
#  Predator_prey_pairs_for_dietary_niche_space.csv [Input 4_1_Predator_dietary_niche_hypervolume.py]
#  Gower_11traits_PCoA_Axis13_inf95_10p_500g.csv [Input 4_1_Predator_dietary_niche_hypervolume.py]
#  Predator_dietary_niche_overlap_hypervolume.csv
#  Predator_dietary_niche_hypervolume.csv
#
#  Final output 
#   Foodweb_metrics.csv
#
#  Supplementary Figures S6, S8, and S14
#

library(cheddar); library(dplyr)
load("Data/Community408_carni_web_doc_inf95_500g_SR3_commcollect.Rdata")

#I. Food web richness and predator-basal prey ratio -------
food.web.dt <- data.frame()
for (i in 1:length(Tropic.carni.webs.doc)){
  
  food.web.dt <- rbind(food.web.dt,
                       data.frame(Community = Tropic.carni.webs.doc[[i]]$properties$Community_unify,
                                  Region = Tropic.carni.webs.doc[[i]]$properties$Region,
                                  #Documented links
                                  N.Sp.doc = NumberOfNodes(Tropic.carni.webs.doc[[i]]),
                                  N.Sp.c.doc = length(unique(Tropic.carni.webs.doc[[i]]$trophic.links$consumer)),
                                  N.Sp.r.doc = length(unique(Tropic.carni.webs.doc[[i]]$trophic.links$resource)),
                                  N.Sp.basal.doc = length(BasalNodes(Tropic.carni.webs.doc[[i]])),
                                  Per.basal.doc = round(FractionBasalNodes(Tropic.carni.webs.doc[[i]]),2),
                                  Per.inter.doc = round(FractionIntermediateNodes(Tropic.carni.webs.doc[[i]]),2),
                                  Per.top.doc = round(FractionTopLevelNodes(Tropic.carni.webs.doc[[i]]),2),
                                  N.links.doc = NumberOfTrophicLinks(Tropic.carni.webs.doc[[i]]),
                                
                                  #Documented + potential links (cosine similarity0.95 body mass +- 10%):
                                  N.Sp.inf95.10p = NumberOfNodes(Tropic.carni.webs.inf95.10p[[i]]),
                                  N.Sp.c.inf95.10p = length(unique(Tropic.carni.webs.inf95.10p[[i]]$trophic.links$consumer)),
                                  N.Sp.r.inf95.10p = length(unique(Tropic.carni.webs.inf95.10p[[i]]$trophic.links$resource)),
                                  N.Sp.basal.inf95.10p = length(BasalNodes(Tropic.carni.webs.inf95.10p[[i]])),
                                  Per.basal.inf95.10p = round(FractionBasalNodes(Tropic.carni.webs.inf95.10p[[i]]),2),
                                  Per.inter.inf95.10p = round(FractionIntermediateNodes(Tropic.carni.webs.inf95.10p[[i]]),2),
                                  Per.top.inf95.10p = round(FractionTopLevelNodes(Tropic.carni.webs.inf95.10p[[i]]),2),
                                  N.links.inf95.10p = NumberOfTrophicLinks(Tropic.carni.webs.inf95.10p[[i]])))
                                  
}
food.web.dt$cr.ratio <- food.web.dt$N.Sp.c.inf95.10p/food.web.dt$N.Sp.basal.inf95.10p

#II. Mean body mass of predator and prey -------
load("Data/Region_sp_intx_500g.RData")
## Traits from Soria et al. (2021) COMBINE dataset https://doi.org/10.1002/ecy.3344
mamm.func <- read.csv("./trait_data_imputed.csv")
mamm.func <- mamm.func[!duplicated(mamm.func$phylacine_binomial),]
mamm.func$adult_mass_g.log <- log10(mamm.func$adult_mass_g)

#prepare body mass in g 
mamm.func.r <- mamm.func[,c("phylacine_binomial","adult_mass_g.log") ]; colnames(mamm.func.r) <- c("phylacine_binomial","adult_mass_g.log.r")
mamm.func.c <- mamm.func[,c("phylacine_binomial","adult_mass_g.log") ]; colnames(mamm.func.c) <- c("phylacine_binomial","adult_mass_g.log.c")

## collect predator and prey species per food web 
comm.intx.inf95.10p <- c()
for (i in 1:length(Tropic.carni.webs.inf95.10p)) {
  #Collect all interactions in a local food web  
  comm.intx.inf95.10p <- rbind(comm.intx.inf95.10p, 
                                       data.frame(Tropic.carni.webs.inf95.10p[[i]][3],
                                                  Community = Tropic.carni.webs.inf95.10p[[i]][2]$properties$Community_unify,
                                                  Region = Tropic.carni.webs.inf95.10p[[i]][2]$properties$Region))
}
length(unique(comm.intx.inf95.10p$Community))

comm.intx.inf95.10p.c <- left_join(unique(comm.intx.inf95.10p[,c("Region", "Community", "trophic.links.consumer")]), mamm.func.c, join_by(trophic.links.consumer==phylacine_binomial))
comm.intx.inf95.10p.r <- left_join(unique(comm.intx.inf95.10p[,c("Region", "Community", "trophic.links.resource")]), mamm.func.r, join_by(trophic.links.resource==phylacine_binomial))

comm.intx.inf95.10p.c.mean <- comm.intx.inf95.10p.c %>% group_by(Community) %>% summarise(Region = unique(Region), mass.log.c = round(mean(adult_mass_g.log.c),3))
comm.intx.inf95.10p.r.mean <- comm.intx.inf95.10p.r %>% group_by(Community) %>% summarise(Region = unique(Region), mass.log.r = round(mean(adult_mass_g.log.r),3))

hist(comm.intx.inf95.10p.c.mean$mass.log.c)

food.web.dt.summary <- left_join(food.web.dt, comm.intx.inf95.10p.c.mean, join_by(Community == Community, Region == Region))
food.web.dt.summary <- left_join(food.web.dt.summary, comm.intx.inf95.10p.r.mean, join_by(Community == Community, Region == Region))

#III. Normalized generality and normalized vulnerability  -----------------------------------------
### comparisons between binary and probabilistic links
library(performance)
##1. Prepare predator-prey pairs and non-interact species in regional species pools ---------
###all species pair based on regional species pool
metaweb.cosine.similarity.raw <- read.csv("Data/Regional_carni_links_500g_cosine_similarity_in95_10p_all.csv")
unique.cosine.pair <- unique(metaweb.cosine.similarity.raw[,c("region","consumer", "potential_resources")])
####total 432 trait-matched prey species based on >= 0.95 cosine similarity +- 10% body mass 
metaweb.cosine.similarity.potential.prey <- metaweb.cosine.similarity.raw[metaweb.cosine.similarity.raw$similarity  >= 0.95 & 
                                                                            metaweb.cosine.similarity.raw$body.ratio.prey.pairs >= 0.9 & metaweb.cosine.similarity.raw$body.ratio.prey.pairs <= 1.1,]
nrow(metaweb.cosine.similarity.potential.prey); length(unique(metaweb.cosine.similarity.potential.prey$potential_resources))  

unique.cosine.pair.potential.prey <- unique(metaweb.cosine.similarity.potential.prey[,c("region","consumer", "potential_resources")])
length(unique(unique.cosine.pair.potential.prey$potential_resources))  

###regional trophic link
metaweb.trophic.link.inf95.10p <- data.frame()
for (i in 1:length(Tropic.carni.webs.inf95.10p)) {
  #Collect all interactions in a local food web  
  metaweb.trophic.link.inf95.10p <- rbind(metaweb.trophic.link.inf95.10p, 
                                          data.frame(Tropic.carni.webs.inf95.10p[[i]][3],
                                                     Region = Tropic.carni.webs.inf95.10p[[i]][2]$properties$Region))
}
length(unique(metaweb.trophic.link.inf95.10p$trophic.links.consumer)) #Total 64 predator species 
length(unique(metaweb.trophic.link.inf95.10p$trophic.links.resource)) #Total 423 prey species (used for global dietary niche)
###check how many prey and predator species within each realm meta-web
metaweb.trophic.link.inf95.10p %>% group_by(Region) %>% summarise(N.c.sp = length(unique(trophic.links.consumer)),
                                                                  N.r.sp = length(unique(trophic.links.resource)))

metaweb.trophic.link.inf95.10p.unique <- metaweb.trophic.link.inf95.10p[!duplicated(metaweb.trophic.link.inf95.10p),]
metaweb.trophic.link.inf95.10p.unique$trophic_id <- seq(1, nrow(metaweb.trophic.link.inf95.10p.unique), by = 1)
nrow(metaweb.trophic.link.inf95.10p.unique) #1295 regional unique predator-prey
length(unique(c(metaweb.trophic.link.inf95.10p.unique$trophic.links.resource, 
                metaweb.trophic.link.inf95.10p.unique$trophic.links.consumer))) #446 unique species
#export unique predator-prey pairs to estimate predator dietnary niche breath and overlap in 
write.csv(metaweb.trophic.link.inf95.10p.unique, "Predator_prey_pairs_for_dietary_niche_space.csv",row.names = F)


##2. Prepare body mass in g to join body mass ratio between predator and prey ----------------
mamm.func.r.g <- mamm.func[,c("phylacine_binomial","adult_mass_g") ]; colnames(mamm.func.r.g) <- c("phylacine_binomial","adult_mass_g.r")
mamm.func.c.g <- mamm.func[,c("phylacine_binomial","adult_mass_g") ]; colnames(mamm.func.c.g) <- c("phylacine_binomial","adult_mass_g.c")

metaweb.cosine.similarity.potential.prey <- left_join(metaweb.cosine.similarity.potential.prey, mamm.func.c.g, join_by(consumer == phylacine_binomial))
metaweb.cosine.similarity.potential.prey <- left_join(metaweb.cosine.similarity.potential.prey, mamm.func.r.g, join_by(resource == phylacine_binomial))


##3. Using a predator species as example: Panthera pardus in Afrotropical meta-web -----------------------
# presence and absence of predator prey interaction
afr.Panthera.pardus <- data.frame(meta.web.cosine.prey  = unique(unique.cosine.pair[unique.cosine.pair$region == "afr",]$potential_resources),
                                  meta.web.predator = "Panthera pardus",
                                  potential.interaction = 0)

afr.Panthera.pardus[afr.Panthera.pardus$meta.web.cosine.prey %in% 
                      unique.cosine.pair.potential.prey[unique.cosine.pair.potential.prey$region == "afr" & unique.cosine.pair.potential.prey$consumer == "Panthera pardus",]$potential_resources,]$potential.interaction <- 1

afr.Panthera.pardus.mass <- left_join(afr.Panthera.pardus, mamm.func.r.g, join_by(meta.web.cosine.prey == phylacine_binomial))
afr.Panthera.pardus.mass <- left_join(afr.Panthera.pardus.mass, mamm.func.c.g, join_by(meta.web.predator == phylacine_binomial))
afr.Panthera.pardus.mass$log_ratio <- log(afr.Panthera.pardus.mass$adult_mass_g.c/afr.Panthera.pardus.mass$adult_mass_g.r)

## Body size model (Rohr et al. 2010 Am Nat)
model.2 <- glm(
  potential.interaction ~ log_ratio + I(log_ratio^2),
  family = binomial(link = "logit"),
  data = afr.Panthera.pardus.mass
)
coef(model.2)

afr.Panthera.pardus.mass$P_hat <- predict(model.2, type = "response")
afr.Panthera.pardus.mass$interaction.binary <- factor(afr.Panthera.pardus.mass$potential.interaction)

ggplot(afr.Panthera.pardus.mass, aes(log_ratio, P_hat, color = interaction.binary)) + geom_point() + 
  geom_rug(data = afr.Panthera.pardus.mass[afr.Panthera.pardus.mass$interaction == 1,], aes(log_ratio, P_hat, color = interaction.binary), sides="b") + theme_bw() +
  scale_color_manual(values = c("grey60", "black")) +
  labs(x = "log(predator.mass/prey.mass)", y = "Probability", title = "Afrotropic: Panthera pardus\nRegional meta-web", color ="interaction")


##4.Body size model: Rohr et al. 2010 The American Naturalist ---------------
analyszed.predator <- metaweb.trophic.link.inf95.10p.unique

###a. Neotropic ------------------------
##predator and prey list 
meta.web.predator <- unique(unique.cosine.pair[unique.cosine.pair$region == "neo",]$consumer)
meta.web.cosine.prey  <- unique(unique.cosine.pair[unique.cosine.pair$region == "neo",]$potential_resources)

neo.body.size.model.coeff.dt <- c()
neo.interaction.probability <- c()

for (i in 1:length(meta.web.predator)) {
  
  meta.web <- data.frame(meta.web.cosine.prey  = meta.web.cosine.prey,
                         meta.web.predator = meta.web.predator[i],
                         potential.interaction = 0)
  
  #potential ecological similary prey without considering co-occurrence
  meta.web[meta.web$meta.web.cosine.prey %in% 
             unique.cosine.pair.potential.prey[unique.cosine.pair.potential.prey$region == "neo" & unique.cosine.pair.potential.prey$consumer == meta.web.predator[i],]$potential_resources,]$potential.interaction <- 1
  
  meta.web.mass <- left_join(meta.web, mamm.func.r.g, join_by(meta.web.cosine.prey == phylacine_binomial))
  meta.web.mass <- left_join(meta.web.mass, mamm.func.c.g, join_by(meta.web.predator == phylacine_binomial))
  meta.web.mass$log_ratio <- log(meta.web.mass$adult_mass_g.c/meta.web.mass$adult_mass_g.r)
  
  #logit body-size model 
  body.size.model <- glm(
    potential.interaction ~ log_ratio + I(log_ratio^2),
    family = binomial(link = "logit"),
    data = meta.web.mass
  )
  
  body.size.model.summary <- summary(body.size.model) 
  #dataframe to save predator-specific body size model coefficient
  neo.body.size.model.coeff.dt <- rbind(neo.body.size.model.coeff.dt,
                                        data.frame(Region = "Neotropic",
                                                   consumer = meta.web.predator[i],
                                                   interacting.prey = length(unique((meta.web.mass[meta.web.mass$potential.interaction == 1,]$meta.web.cosine.prey))),
                                                   alpha = body.size.model.summary$coefficients[1,1],
                                                   alpha.se = body.size.model.summary$coefficients[1,2],
                                                   beta = body.size.model.summary$coefficients[2,1],
                                                   beta.se = body.size.model.summary$coefficients[2,2],
                                                   gamma = body.size.model.summary$coefficients[3,1],
                                                   gamma.se = body.size.model.summary$coefficients[3,2],
                                                   pseudo_R2 = r2(body.size.model),
                                                   converge = body.size.model$converged))
  
  meta.web.mass$P_hat <- predict(body.size.model, type = "response")
  
  predator.prey.bodysize.p <- meta.web.mass
  #dataframe to save predator-prey pair specific probabilistic link
  neo.interaction.probability <- rbind(neo.interaction.probability,
                                       data.frame(Region = "Neotropic",
                                                  predator.prey.bodysize.p))
  
  meta.web.mass$interaction.binary <- factor(meta.web.mass$potential.interaction, levels = rev(c("1", "0")))
  ### visual by x-axis with body mass ratio of predator mass to prey mass
  ggplot(meta.web.mass, aes(log_ratio, P_hat, color = interaction.binary)) + geom_point() + 
    geom_rug(data = meta.web.mass[meta.web.mass$potential.interaction == 1,], aes(log_ratio, P_hat, color = interaction.binary), sides="b") + theme_bw() +
    scale_color_manual(values = c("grey60", "black")) +
    labs(x = "log(predator.mass/prey.mass)", y = "Probability", title = paste0("Neotropic: ", meta.web.predator[i]), color ="interaction")
  ### visual by x-axis with prey body mass 
  ggplot(meta.web.mass, aes(log10(adult_mass_g.r), P_hat, color = interaction.binary)) + geom_point() + 
    geom_point(data = meta.web.mass[meta.web.mass$potential.interaction == 1,], aes(log10(adult_mass_g.r), P_hat, color = interaction.binary), color = "black") + 
    geom_rug(data = meta.web.mass[meta.web.mass$potential.interaction == 1,], aes(log10(adult_mass_g.r), P_hat, color = interaction.binary), sides="b") + theme_bw() +
    scale_color_manual(values = c("grey60", "black")) +
    labs(x = "log10(Prey mass)", y = "Probability", title = paste0("Neotropic: ", meta.web.predator[i]), color ="interaction")
  
}

###b. Afrotropic ------------------------
##predator and prey list 
meta.web.predator <- unique(unique.cosine.pair[unique.cosine.pair$region == "afr",]$consumer)
meta.web.cosine.prey  <- unique(unique.cosine.pair[unique.cosine.pair$region == "afr",]$potential_resources)

afr.body.size.model.coeff.dt <- c()
afr.interaction.probability <- c()
for (i in 1:length(meta.web.predator)) {
  print(i)
  meta.web <- data.frame(meta.web.cosine.prey  = meta.web.cosine.prey,
                         meta.web.predator = meta.web.predator[i],
                         potential.interaction = 0)
  
  #potential ecological similary prey without considering co-occurrence
  meta.web[meta.web$meta.web.cosine.prey %in% 
             unique.cosine.pair.potential.prey[unique.cosine.pair.potential.prey$region == "afr" & unique.cosine.pair.potential.prey$consumer == meta.web.predator[i],]$potential_resources,]$potential.interaction <- 1
  
  meta.web.mass <- left_join(meta.web, mamm.func.r.g, join_by(meta.web.cosine.prey == phylacine_binomial))
  meta.web.mass <- left_join(meta.web.mass, mamm.func.c.g, join_by(meta.web.predator == phylacine_binomial))
  meta.web.mass$log_ratio <- log(meta.web.mass$adult_mass_g.c/meta.web.mass$adult_mass_g.r)
  
  #logit body-size model 
  body.size.model <- glm(
    potential.interaction ~ log_ratio + I(log_ratio^2),
    family = binomial(link = "logit"),
    data = meta.web.mass
  )
  
  body.size.model.summary <- summary(body.size.model) 
  #dataframe to save predator-specific body size model coefficient
  afr.body.size.model.coeff.dt <- rbind(afr.body.size.model.coeff.dt,
                                        data.frame(Region = "Afrotropic",
                                                   consumer = meta.web.predator[i],
                                                   interacting.prey = length(unique((meta.web.mass[meta.web.mass$potential.interaction == 1,]$meta.web.cosine.prey))),
                                                   alpha = body.size.model.summary$coefficients[1,1],
                                                   alpha.se = body.size.model.summary$coefficients[1,2],
                                                   beta = body.size.model.summary$coefficients[2,1],
                                                   beta.se = body.size.model.summary$coefficients[2,2],
                                                   gamma = body.size.model.summary$coefficients[3,1],
                                                   gamma.se = body.size.model.summary$coefficients[3,2],
                                                   pseudo_R2 = r2(body.size.model),
                                                   converge = body.size.model$converged))
  
  meta.web.mass$P_hat <- predict(body.size.model, type = "response")
  predator.prey.bodysize.p <- meta.web.mass
  #dataframe to save predator-prey pair specific probabilistic link
  afr.interaction.probability <- rbind(afr.interaction.probability,
                                       data.frame(Region = "Afrotropic",
                                                  predator.prey.bodysize.p))
  
  meta.web.mass$interaction.binary <- factor(meta.web.mass$potential.interaction, levels = rev(c("1", "0")))
  ### visual by x-axis with body mass ratio of predator mass to prey mass
  ggplot(meta.web.mass, aes(log_ratio, P_hat, color = interaction.binary)) + geom_point() + 
    geom_rug(data = meta.web.mass[meta.web.mass$potential.interaction == 1,], aes(log_ratio, P_hat, color = interaction.binary), sides="b") + theme_bw() +
    scale_color_manual(values = c("grey60", "black")) +
    labs(x = "log(predator.mass/prey.mass)", y = "Probability", title = paste0("Afrotropic: ", meta.web.predator[i]), color ="interaction")
  ### visual by x-axis with prey body mass 
  ggplot(meta.web.mass, aes(log10(adult_mass_g.r), P_hat, color = interaction.binary)) + geom_point() + 
    geom_point(data = meta.web.mass[meta.web.mass$potential.interaction == 1,], aes(log10(adult_mass_g.r), P_hat, color = interaction.binary), color = "black") + 
    geom_rug(data = meta.web.mass[meta.web.mass$potential.interaction == 1,], aes(log10(adult_mass_g.r), P_hat, color = interaction.binary), sides="b") + theme_bw() +
    scale_color_manual(values = c("grey60", "black")) +
    labs(x = "log10(Prey mass)", y = "Probability", title = paste0("Afrotropic: ", meta.web.predator[i]), color ="interaction")
  
}

###c. Indomalaya ------------------------
##predator and prey list 
meta.web.predator <- unique(unique.cosine.pair[unique.cosine.pair$region == "ind",]$consumer)
meta.web.cosine.prey  <- unique(unique.cosine.pair[unique.cosine.pair$region == "ind",]$potential_resources)

ind.body.size.model.coeff.dt <- c()
ind.interaction.probability <- c()
for (i in 1:length(meta.web.predator)) {
  
  meta.web <- data.frame(meta.web.cosine.prey  = meta.web.cosine.prey,
                         meta.web.predator = meta.web.predator[i],
                         potential.interaction = 0)
  
  #potential ecological similary prey without considering co-occurrence
  meta.web[meta.web$meta.web.cosine.prey %in% 
             unique.cosine.pair.potential.prey[unique.cosine.pair.potential.prey$region == "ind" & unique.cosine.pair.potential.prey$consumer == meta.web.predator[i],]$potential_resources,]$potential.interaction <- 1
  
  meta.web.mass <- left_join(meta.web, mamm.func.r.g, join_by(meta.web.cosine.prey == phylacine_binomial))
  meta.web.mass <- left_join(meta.web.mass, mamm.func.c.g, join_by(meta.web.predator == phylacine_binomial))
  meta.web.mass$log_ratio <- log(meta.web.mass$adult_mass_g.c/meta.web.mass$adult_mass_g.r)
  
  #logit body-size model 
  body.size.model <- glm(
    potential.interaction ~ log_ratio + I(log_ratio^2),
    family = binomial(link = "logit"),
    data = meta.web.mass
  )
  
  body.size.model.summary <- summary(body.size.model) 
  #dataframe to save predator-specific body size model coefficient
  ind.body.size.model.coeff.dt <- rbind(ind.body.size.model.coeff.dt,
                                        data.frame(Region = "Indomalayan",
                                                   consumer = meta.web.predator[i],
                                                   interacting.prey = length(unique((meta.web.mass[meta.web.mass$potential.interaction == 1,]$meta.web.cosine.prey))),
                                                   alpha = body.size.model.summary$coefficients[1,1],
                                                   alpha.se = body.size.model.summary$coefficients[1,2],
                                                   beta = body.size.model.summary$coefficients[2,1],
                                                   beta.se = body.size.model.summary$coefficients[2,2],
                                                   gamma = body.size.model.summary$coefficients[3,1],
                                                   gamma.se = body.size.model.summary$coefficients[3,2],
                                                   pseudo_R2 = r2(body.size.model),
                                                   converge = body.size.model$converged))
  
  meta.web.mass$P_hat <- predict(body.size.model, type = "response")
  #P.equation <- coef(body.size.model)[1] + coef(body.size.model)[2]*meta.web.mass$log_ratio + coef(body.size.model)[3]*(meta.web.mass$log_ratio^2)
  #meta.web.mass$P_hat_self <- (exp(P.equation)/(1+exp(P.equation))) #same with predict function
  
  predator.prey.bodysize.p <- meta.web.mass
  #dataframe to save predator-prey pair specific probabilistic link
  ind.interaction.probability <- rbind(ind.interaction.probability,
                                       data.frame(Region = "Indomalayan",
                                                  predator.prey.bodysize.p))
  
  meta.web.mass$interaction.binary <- factor(meta.web.mass$potential.interaction, levels = rev(c("1", "0")))
  ### visual by x-axis with body mass ratio of predator mass to prey mass
  ggplot(meta.web.mass, aes(log_ratio, P_hat, color = interaction.binary)) + geom_point() + 
    geom_rug(data = meta.web.mass[meta.web.mass$potential.interaction == 1,], aes(log_ratio, P_hat, color = interaction.binary), sides="b") + theme_bw() +
    scale_color_manual(values = c("grey60", "black")) +
    labs(x = "log(predator.mass/prey.mass)", y = "Probability", title = paste0("Indomalayan: ", meta.web.predator[i]), color ="interaction")
  ### visual by x-axis with prey body mass 
  ggplot(meta.web.mass, aes(log10(adult_mass_g.r), P_hat, color = interaction.binary)) + geom_point() +
    geom_point(data = meta.web.mass[meta.web.mass$potential.interaction == 1,], aes(log10(adult_mass_g.r), P_hat, color = interaction.binary), color = "black") + 
    geom_rug(data = meta.web.mass[meta.web.mass$potential.interaction == 1,], aes(log10(adult_mass_g.r), P_hat, color = interaction.binary), sides="b") + theme_bw() +
    scale_color_manual(values = c("grey60", "black")) +
    labs(x = "log10(Prey mass)", y = "Probability", title = paste0("Indomalayan: ", meta.web.predator[i]), color ="interaction")
  
}

##5. Allometric niche model: Williams et al 2010. PLOS One --------------------------
### assumption that predators follow a global curve of symmetric preference with optimal mass and niche width

interaction.p <- rbind(neo.interaction.probability, afr.interaction.probability, ind.interaction.probability)
interaction.p$log10.adult_mass_g.r <- log10(interaction.p$adult_mass_g.r)
interaction.p$log10.adult_mass_g.c <- log10(interaction.p$adult_mass_g.c)

species.mass <- unique(c(interaction.p$log10.adult_mass_g.r, interaction.p$log10.adult_mass_g.c))
logmin <- min(species.mass)
logmax <- max(species.mass)

interaction.p$adult_mass_g.r_n <- (interaction.p$log10.adult_mass_g.r - logmin)/(logmax - logmin)
interaction.p$adult_mass_g.c_n <- (interaction.p$log10.adult_mass_g.c - logmin)/(logmax - logmin)

ggplot(interaction.p, aes(adult_mass_g.r_n)) + geom_histogram(fill = "grey80", position = "identity") + 
  geom_histogram(data = interaction.p, aes(adult_mass_g.c_n), fill = "grey20", position = "identity") + theme_bw() +
  labs( x = "n: mass position of prey", y = "n: mass position of predator")


region.list <- unique(interaction.p$Region)
Allometric_niche_model_p <- c()
Allometric_niche_model_coeff <- c()

for (z in 1:length(region.list)) {
  
  region.predator <- sort(unique(interaction.p[interaction.p$Region == region.list[z] ,]$meta.web.predator))
  region.predator <- region.predator[region.predator %in% unique(analyszed.predator[analyszed.predator$Region == region.list[z],]$trophic.links.consumer)]
  print(length(region.predator))
  
  for (j in 1:length(region.predator)) {
    
    region.predator.available.prey <-  interaction.p[interaction.p$Region == region.list[z] & interaction.p$meta.web.predator == region.predator[j],]
    region.predator.prey.range <-  region.predator.available.prey[region.predator.available.prey$potential.interaction == 1,]
    c.j <- mean(region.predator.prey.range$adult_mass_g.r_n)
    r.j <- sd(region.predator.prey.range$adult_mass_g.r_n)
    
    #if only one prey species, set sd = 0.0001
    if (is.na(r.j) == T) {
      
      r.j <- 0.0001
      region.predator.available.prey$allometric.niche.model.p <- exp(-(((region.predator.available.prey$adult_mass_g.r_n - c.j)/(r.j/2))^2))
      
    } else {
      
      region.predator.available.prey$allometric.niche.model.p <- exp(-(((region.predator.available.prey$adult_mass_g.r_n - c.j)/(r.j/2))^2))
      
    }
    
    region.predator.available.prey$potential.interaction <- factor(region.predator.available.prey$potential.interaction)
    ggplot(region.predator.available.prey, aes(log10.adult_mass_g.r, allometric.niche.model.p, color = potential.interaction))  + geom_point() + theme_bw() +
      geom_point(data = region.predator.available.prey[region.predator.available.prey$potential.interaction == 1,], aes(log10.adult_mass_g.r, allometric.niche.model.p), color = "black") +
      scale_color_manual(values = c("grey60", "black")) + labs(x = "log10(Prey mass)", y = "Interaction probability", color = "Potential\ninteraction",
                                                               title = paste0(region.list[z], ": ",region.predator[j]))
    
    Allometric_niche_model_p <- rbind(Allometric_niche_model_p, region.predator.available.prey)
    Allometric_niche_model_coeff <- rbind(Allometric_niche_model_coeff, data.frame(Region = region.list[z],
                                                                                   consumer = region.predator[j], 
                                                                                   consumer.diet.position = c.j, 
                                                                                   consumer.diet.range = r.j))
  }
}

##6. Measure normalized generality and vulnerability with binary link and probabilistic links from two models -------
compare.link.p <- Allometric_niche_model_p #P_hat: Body size model vs. allometric.niche.model.p: Allometric niche model
compare.link.p$pair <- paste0(compare.link.p$meta.web.cosine.prey, ".",compare.link.p$meta.web.predator)
compare.link.p.unique <- unique(compare.link.p[,c("Region", "meta.web.cosine.prey", "meta.web.predator", "P_hat","allometric.niche.model.p","pair")])

###a. predator generality and prey vulnerability, and linkage density to normalize measures--------------------
predator.specific.generality <-c()
prey.specific.vulnerability <-c()
web.link.d <- c()

for (i in 1:length(Tropic.carni.webs.inf95.10p)) {
  
  links <- Tropic.carni.webs.inf95.10p[[i]]$trophic.links
  links$pair <- paste0(links$resource, ".",links$consumer)
  
  links.compare.link.p.unique <- left_join(links, compare.link.p.unique[compare.link.p.unique$Region == Tropic.carni.webs.inf95.10p[[i]]$properties$Region, c("Region","pair","P_hat", "allometric.niche.model.p")], join_by(pair == pair))
  links.compare.link.p.unique <- unique(links.compare.link.p.unique)
  
  #food web node richness
  web.S <- length(unique(c(links.compare.link.p.unique$consumer, links.compare.link.p.unique$resource)))
  #food web linkage density based on binary links and probabilistic links
  web.link.d <- rbind(web.link.d,
                      data.frame(Community = Tropic.carni.webs.inf95.10p[[i]]$properties$Community_unify,
                                 link.d = round(length(links.compare.link.p.unique$resource)/web.S,2),
                                 link.d.weighted.body.size.model.p = round(sum(links.compare.link.p.unique$P_hat)/web.S,2),
                                 link.d.weighted.allometric.model.p = round(sum(links.compare.link.p.unique$allometric.niche.model.p)/web.S,2)))
  #predator specific generality (number of prey links) based on binary links and probabilistic links
  predator.specific.generality <- rbind(predator.specific.generality,
                                        data.frame(Region = Tropic.carni.webs.inf95.10p[[i]]$properties$Region, Community = Tropic.carni.webs.inf95.10p[[i]]$properties$Community_unify,
                                                   links.compare.link.p.unique %>% group_by(consumer) %>% summarise(generality.sp = length(unique(resource)),
                                                                                                                   generality.weighted.body.size.model.p = sum(P_hat),
                                                                                                                   generality.weighted.allometric.model.p = sum(allometric.niche.model.p))))
  #prey specific .vulnerability (number of predator links) based on binary links and probabilistic links
  prey.specific.vulnerability <- rbind(prey.specific.vulnerability,
                                       data.frame(Region = Tropic.carni.webs.inf95.10p[[i]]$properties$Region, Community = Tropic.carni.webs.inf95.10p[[i]]$properties$Community_unify,
                                                  links.compare.link.p.unique %>% group_by(resource) %>% summarise(vulnerability.sp = length(unique(consumer)),
                                                                                                                  vulnerability.weighted.body.size.model.p = sum(P_hat),
                                                                                                                  vulnerability.weighted.allometric.model.p = sum(allometric.niche.model.p))))
}

predator.specific.generality.link.d <- left_join(predator.specific.generality, web.link.d, join_by(Community == Community))
prey.specific.vulnerability.link.d <- left_join(prey.specific.vulnerability, web.link.d, join_by(Community == Community))

web.generality <- predator.specific.generality.link.d %>% group_by(Region, Community) %>% summarise(generality.sp = round(mean(generality.sp),1),
                                                                                             generality.weighted.body.size.model.p = round(mean(generality.weighted.body.size.model.p),1),
                                                                                             generality.weighted.allometric.model.p = round(mean(generality.weighted.allometric.model.p),1),
                                                                                             
                                                                                             generality.sp.normalized = round(mean(generality.sp/link.d),1),
                                                                                             generality.weighted.body.size.model.p.normalized  = round(mean(generality.weighted.body.size.model.p/link.d.weighted.body.size.model.p),1),
                                                                                             generality.weighted.allometric.model.p.normalized  = round(mean(generality.weighted.allometric.model.p/link.d.weighted.allometric.model.p),1))


web.vulnerability <- prey.specific.vulnerability.link.d %>% group_by(Region, Community) %>% summarise(vulnerability.sp = round(mean(vulnerability.sp),1),
                                                                                               vulnerability.weighted.body.size.model.p = round(mean(vulnerability.weighted.body.size.model.p),1),
                                                                                               vulnerability.weighted.allometric.model.p = round(mean(vulnerability.weighted.allometric.model.p),1),
                                                                                               
                                                                                               vulnerability.sp.normalized = round(mean(vulnerability.sp/link.d),1),
                                                                                               vulnerability.weighted.body.size.model.p.normalized = round(mean(vulnerability.weighted.body.size.model.p/link.d.weighted.body.size.model.p),1),
                                                                                               vulnerability.weighted.allometric.model.p.normalized = round(mean(vulnerability.weighted.allometric.model.p/link.d.weighted.allometric.model.p),1))

web.generality$Region <- factor(web.generality$Region, levels =c("Neotropic", "Afrotropic", "Indomalayan"))
web.vulnerability$Region <- factor(web.vulnerability$Region, levels =c("Neotropic", "Afrotropic", "Indomalayan"))


food.web.dt.summary <- left_join(food.web.dt.summary, web.link.d, join_by(Community == Community))
food.web.dt.summary <- left_join(food.web.dt.summary, web.generality, join_by(Community == Community, Region == Region))
food.web.dt.summary <- left_join(food.web.dt.summary, web.vulnerability, join_by(Community == Community, Region == Region))


#
# Prepare niche axes to estimate predator dietary niche breadth by prey traits and measure mean predator specific dietary niche breadth and overlap within food web
#
# 
# Input: 
#   1. Region_sp_intx_500g.RData
#   2. cosine_similarity_500g_inferred_regional_carni_links_all.csv
#
# Output: 
#   1. Gower_11traits_PCoA_Axis13_inf95_10p_500g.csv
#

library(dplyr); library(ggplot2); library(ggpubr)
library(mFD)

load("Data/Community408_carni_web_doc_inf95_500g_SR3_commcollect.Rdata")

#IV. Generate prey species in realm meta-web per predator to estimate geometric hypervolumne =====================================
##1. Collect realm meta-web---------------
#call community collection list
metaweb.trophic.link.inf95.10p <- data.frame()
for (i in 1:length(Tropic.carni.webs.inf95.10p)) {
  #Collect all interactions in a local food web  
  metaweb.trophic.link.inf95.10p <- rbind(metaweb.trophic.link.inf95.10p, 
                                          data.frame(Tropic.carni.webs.inf95.10p[[i]][3],
                                                     Region = Tropic.carni.webs.inf95.10p[[i]][2]$properties$Region))
}
length(unique(metaweb.trophic.link.inf95.10p$trophic.links.consumer)) #Total 64 predator species 
length(unique(metaweb.trophic.link.inf95.10p$trophic.links.resource)) #Total 423 prey species (used for global dietary niche)
#Check how many prey and predator species within each realm meta-web
metaweb.trophic.link.inf95.10p %>% group_by(Region) %>% summarise(N.c.sp = length(unique(trophic.links.consumer)),
                                                                  N.r.sp = length(unique(trophic.links.resource)))

#Collect unique regional trophic links
metaweb.trophic.link.inf95.10p.unique <- metaweb.trophic.link.inf95.10p[!duplicated(metaweb.trophic.link.inf95.10p),]
metaweb.trophic.link.inf95.10p.unique$trophic_id <- seq(1, nrow(metaweb.trophic.link.inf95.10p.unique), by = 1)
nrow(metaweb.trophic.link.inf95.10p.unique) #1295 regional unique predator-prey
length(unique(c(metaweb.trophic.link.inf95.10p.unique$trophic.links.resource, 
                metaweb.trophic.link.inf95.10p.unique$trophic.links.consumer))) #446 unique species

### Check consumer that only has one resource species in meta-web----
metaweb.consumer.one.resource.inf95.10p <- metaweb.trophic.link.inf95.10p.unique %>% group_by(trophic.links.consumer, Region) %>%
  summarize(N.Sp.r = length(trophic.links.resource))
nrow(metaweb.consumer.one.resource.inf95.10p[metaweb.consumer.one.resource.inf95.10p$N.Sp.r == 1,]) #9 consumer sp. has only on prey
metaweb.consumer.one.resource.inf95.10p[metaweb.consumer.one.resource.inf95.10p$N.Sp.r > 100,]$trophic.links.consumer

#V. Join functional traits to total prey within three realms to construct total predator dietary niche space ==================================

mamm.func <- read.csv("./trait_data_imputed.csv")
mamm.func <- mamm.func[!duplicated(mamm.func$phylacine_binomial),]
#check scientific name match
metaweb.trophic.link.inf95.10p.unique[!metaweb.trophic.link.inf95.10p.unique$trophic.links.resource %in% mamm.func$phylacine_binomial,]$trophic.links.resource #all match

#check missing values in targeted traits
mamm.func[is.na(mamm.func$adult_mass_g)==T & mamm.func$phylacine_binomial %in% metaweb.trophic.link.inf95.10p.unique$trophic.links.resource,]$phylacine_binomial 
mamm.func[is.na(mamm.func$det_diet_breadth_n)==T & mamm.func$phylacine_binomial %in% metaweb.trophic.link.inf95.10p.unique$trophic.links.resource,]$phylacine_binomial 
mamm.func[is.na(mamm.func$trophic_level)==T & mamm.func$phylacine_binomial %in% metaweb.trophic.link.inf95.10p.unique$trophic.links.resource,]$phylacine_binomial 
mamm.func[is.na(mamm.func$activity_cycle)==T & mamm.func$phylacine_binomial %in% metaweb.trophic.link.inf95.10p.unique$trophic.links.resource,]$phylacine_binomial 
mamm.func[is.na(mamm.func$foraging_stratum)==T & mamm.func$phylacine_binomial %in% metaweb.trophic.link.inf95.10p.unique$trophic.links.resource,]$phylacine_binomial 
mamm.func[is.na(mamm.func$habitat_breadth_n)==T & mamm.func$phylacine_binomial %in% metaweb.trophic.link.inf95.10p.unique$trophic.links.resource,]$phylacine_binomial 
unique(mamm.func[is.na(mamm.func$habitat_breadth_n)==T & mamm.func$phylacine_binomial %in% metaweb.trophic.link.inf95.10p.unique$trophic.links.resource,]$phylacine_binomial)
#fill-in three species without habitat breadth based on IUCN categories
#Fill values from IUCN habitat types 
mamm.func[mamm.func$phylacine_binomial == "Cebus albifrons",]$habitat_breadth_n <- 2 #https://www.iucnredlist.org/species/39951/191703935#habitat-ecology
mamm.func[mamm.func$phylacine_binomial == "Cebus brunneus",]$habitat_breadth_n <- 1 #https://www.iucnredlist.org/species/81237954/17981252#habitat-ecology
mamm.func[mamm.func$phylacine_binomial == "Cebus capucinus",]$habitat_breadth_n <- 2 #https://www.iucnredlist.org/species/81257277/191708164#habitat-ecology
mamm.func[mamm.func$phylacine_binomial == "Loxodonta africana",]$habitat_breadth_n <- 6 #https://www.iucnredlist.org/species/181008073/223031019#habitat-ecology
mamm.func[mamm.func$phylacine_binomial == "Trachypithecus barbei",]$habitat_breadth_n <- 1 #https://www.iucnredlist.org/species/41554/17960144#habitat-ecology

#collect resource species trait
mamm.func <- mamm.func[, c("phylacine_binomial","adult_mass_g","brain_mass_g", "adult_body_length_mm", "dispersal_km",
                           "det_diet_breadth_n", "habitat_breadth_n", "hibernation_torpor", "fossoriality",
                           "foraging_stratum", "activity_cycle","trophic_level")]
mamm.func.unique <- mamm.func[!duplicated(mamm.func$phylacine_binomial),]

mamm.func.unique$adult_mass_g_log <- log(mamm.func.unique$adult_mass_g,10)
mamm.func.unique$brain_mass_g_log <- log(mamm.func.unique$brain_mass_g,10)
mamm.func.unique$adult_body_length_mm_log <- log(mamm.func.unique$adult_body_length_mm,10)
mamm.func.unique$dispersal_km_log <- log(mamm.func.unique$dispersal_km,10)
mamm.func.unique$foraging_stratum <- factor(mamm.func.unique$foraging_stratum)
mamm.func.unique$activity_cycle <- factor(mamm.func.unique$activity_cycle)
mamm.func.unique$trophic_level <- factor(mamm.func.unique$trophic_level)
###Join functional traits to prey species
metaweb.trophic.link.r.inf95.10p.func <- left_join(metaweb.trophic.link.inf95.10p.unique, mamm.func.unique,
                                                   join_by(trophic.links.resource == phylacine_binomial))

#VI. Measure Gower dissimilarity of total prey species within three realms =========
##1. Calculate pairwise Gower dissimilarity matrix ----
metaweb.trophic.link.r.inf95.10p.func.gower <- metaweb.trophic.link.r.inf95.10p.func[!duplicated(metaweb.trophic.link.r.inf95.10p.func$trophic.links.resource),
                                                                                     c("adult_mass_g_log", "brain_mass_g_log", "adult_body_length_mm_log",
                                                                                       "dispersal_km_log", "det_diet_breadth_n", "habitat_breadth_n",
                                                                                       "hibernation_torpor", "fossoriality", "foraging_stratum", "activity_cycle", "trophic_level")]
nrow(metaweb.trophic.link.r.inf95.10p.func.gower) #423 prey species
rownames(metaweb.trophic.link.r.inf95.10p.func.gower) <- metaweb.trophic.link.r.inf95.10p.func[!duplicated(metaweb.trophic.link.r.inf95.10p.func$trophic.links.resource),]$trophic.links.resource

trait_type <- data.frame(trait_name = colnames(metaweb.trophic.link.r.inf95.10p.func.gower),
                         trait_type =c("Q","Q","Q","Q","Q","Q","Q","Q", "N","N","N"))

metaweb.trophic.link.r.inf95.10p.func.gower.dist <- as.matrix(funct.dist(metaweb.trophic.link.r.inf95.10p.func.gower, trait_type, 
                                                                         metric        = "gower",
                                                                         scale_euclid  = "scale_center",
                                                                         ordinal_var   = "classic",
                                                                         weight_type   = "equal",
                                                                         stop_if_NA    = TRUE))
metaweb.trophic.link.r.inf95.10p.func.gower.dist <- metaweb.trophic.link.r.inf95.10p.func.gower.dist*10

##2. Ordination by PCoA =====
library(ape)
pcoa.inf <- pcoa(metaweb.trophic.link.r.inf95.10p.func.gower.dist)

barplot(pcoa.inf$values$Eigenvalues)
(pcoa.inf$values$Relative_eig[1]+pcoa.inf$values$Relative_eig[2]+pcoa.inf$values$Relative_eig[3]+pcoa.inf$values$Relative_eig[4]+pcoa.inf$values$Relative_eig[5] + pcoa.inf$values$Relative_eig[6] 
  +pcoa.inf$values$Relative_eig[7] +pcoa.inf$values$Relative_eig[8] + pcoa.inf$values$Relative_eig[9] 
  +pcoa.inf$values$Relative_eig[10] + pcoa.inf$values$Relative_eig[11]+ pcoa.inf$values$Relative_eig[12] +
    pcoa.inf$values$Relative_eig[13]
)/sum(pcoa.inf$values$Relative_eig[pcoa.inf$values$Relative_eig > 0]) # 95%; axis.1 = 40.32%, axis.2 = 15.94%

## Plot accumulated explanation rate (Figure S14)
acc.relative.eig <- data.frame(axis = paste0("Axis_",seq(1, length(pcoa.inf$values$Relative_eig), by = 1)), relative_eig_value = pcoa.inf$values$Relative_eig)
acc.relative.eig <- acc.relative.eig[acc.relative.eig$relative_eig_value > 0,]
acc.relative.eig$acc_axis <- seq(from = 1, to = length(acc.relative.eig$axis), by = 1)
acc.relative.eig$acc_explain_variation <- 0
acc.relative.eig$acc_explain_variation[1] <- acc.relative.eig$relative_eig_value[1]/sum(acc.relative.eig$relative_eig_value)

## Supplementary Figure S14 ----------
### Association of cumulative explained variation with the number of PCoA component based on the pairwise Gower-distance of eleven ecological traits among the 423 prey species. 
for (i in 2:length(acc.relative.eig$axis)) {
  acc.relative.eig$acc_explain_variation[i] <- acc.relative.eig$acc_explain_variation[i-1] + (acc.relative.eig$relative_eig_value[i]/sum(acc.relative.eig$relative_eig_value))
}
ggplot(acc.relative.eig[acc.relative.eig$relative_eig_value > 0,], aes(acc_axis, acc_explain_variation*100)) + 
  geom_vline(xintercept = 13, lty=2, col = "grey60") + geom_hline(yintercept = 95.5, lty=2, col = "grey60") + geom_point() +   
  theme_bw() +  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14),
                      legend.title = element_text(size = 16),legend.text = element_text(size = 14),
                      panel.grid = element_blank(),legend.position  ="none") + scale_x_continuous(expand = c(0.08, 0)) +
  labs(x = "Number of PCoA component", y = "Cumulative explained variation (%)")

##Extract PCoA traits to construct multidimensional trait space
traits.inf <- data.frame(pcoa.inf$vectors[,1:13])
traits.inf$Species <- rownames(traits.inf)
write.csv(traits.inf, "Data/Gower_11traits_PCoA_Axis13_inf95_10p_500g.csv")


#VII. Community-mean predator geometric resource hypervolume ====================================
##1. Mean predator dietary niche breadth within food web --------
predator.regional.r.alpha.geometric.inf <- read.csv("/Data/Predator_dietary_niche_hypervolume.csv")
hypervolume.alpha.comm.inf <- data.frame()
for (i in 1: length(Tropic.carni.webs.inf95.10p)) {
  if(Tropic.carni.webs.inf95.10p[[i]]$properties$Region == "Afrotropic") {
    comm.con.sp <- unique(Tropic.carni.webs.inf95.10p[[i]]$trophic.links$consumer)
    comm.alpha <- predator.regional.r.alpha.geometric.inf[predator.regional.r.alpha.geometric.inf$region == "Afrotropic" &
                                                            predator.regional.r.alpha.geometric.inf$consumer %in% comm.con.sp,] 
    hypervolume.alpha.comm.inf <- rbind(hypervolume.alpha.comm.inf, data.frame(Community = Tropic.carni.webs.inf95.10p[[i]]$properties$Community_unify,
                                                                               Region = Tropic.carni.webs.inf95.10p[[i]]$properties$Region,
                                                                               N.Sp.c.hv.alpha.inf = length(comm.alpha$consumer),
                                                                               resource.hv.alpha.m.inf = round(mean(comm.alpha$hypervolume),2)))
    
  } else if (Tropic.carni.webs.inf95.10p[[i]]$properties$Region == "Indomalayan") {
    comm.con.sp <- unique(Tropic.carni.webs.inf95.10p[[i]]$trophic.links$consumer)
    comm.alpha <- predator.regional.r.alpha.geometric.inf[predator.regional.r.alpha.geometric.inf$region == "Indomalayan" &
                                                            predator.regional.r.alpha.geometric.inf$consumer %in% comm.con.sp,] 
    hypervolume.alpha.comm.inf <- rbind(hypervolume.alpha.comm.inf, data.frame(Community = Tropic.carni.webs.inf95.10p[[i]]$properties$Community_unify,
                                                                               Region = Tropic.carni.webs.inf95.10p[[i]]$properties$Region,
                                                                               N.Sp.c.hv.alpha.inf = length(comm.alpha$consumer),
                                                                               resource.hv.alpha.m.inf = round(mean(comm.alpha$hypervolume),2)))
  } else {
    comm.con.sp <- unique(Tropic.carni.webs.inf95.10p[[i]]$trophic.links$consumer)
    comm.alpha <-  predator.regional.r.alpha.geometric.inf[predator.regional.r.alpha.geometric.inf$region == "Neotropic" &
                                                             predator.regional.r.alpha.geometric.inf$consumer %in% comm.con.sp,] 
    hypervolume.alpha.comm.inf <- rbind(hypervolume.alpha.comm.inf, data.frame(Community = Tropic.carni.webs.inf95.10p[[i]]$properties$Community_unify,
                                                                               Region = Tropic.carni.webs.inf95.10p[[i]]$properties$Region,
                                                                               N.Sp.c.hv.alpha.inf = length(comm.alpha$consumer),
                                                                               resource.hv.alpha.m.inf = round(mean(comm.alpha$hypervolume),2)))
  }
}
nrow(hypervolume.alpha.comm.inf[hypervolume.alpha.comm.inf$N.Sp.c.hv.alpha.inf > 2,])

##2. Mean predator dietary niche overlap within food web --------
predator.regional.r.beta.geometric.inf <- read.csv("Data/Predator_dietary_niche_overlap_hypervolume.csv")
hypervolume.beta.comm.inf <- data.frame()
for (i in 1: length(Tropic.carni.webs.inf95.10p)) {
  if(Tropic.carni.webs.inf95.10p[[i]]$properties$Region == "Afrotropic") {
    comm.con.sp <- unique(Tropic.carni.webs.inf95.10p[[i]]$trophic.links$consumer)
    comm.beta <- predator.regional.r.beta.geometric.inf[predator.regional.r.beta.geometric.inf$region == "Afrotropic" &
                                                          predator.regional.r.beta.geometric.inf$consumer1 %in% comm.con.sp & 
                                                          predator.regional.r.beta.geometric.inf$consumer2 %in% comm.con.sp,] 
    hypervolume.beta.comm.inf <- rbind(hypervolume.beta.comm.inf, data.frame(Community = Tropic.carni.webs.inf95.10p[[i]]$properties$Community_unify,
                                                                             Region = Tropic.carni.webs.inf95.10p[[i]]$properties$Region,
                                                                             N.Sp.c.hv.beta.inf = length(comm.con.sp),
                                                                             resource.hv.beta.m.inf = round(mean(comm.beta$overlap_hypervolume),2)))
    
  } else if (Tropic.carni.webs.inf95.10p[[i]]$properties$Region == "Indomalayan") {
    comm.con.sp <- unique(Tropic.carni.webs.inf95.10p[[i]]$trophic.links$consumer)
    comm.beta <- predator.regional.r.beta.geometric.inf[predator.regional.r.beta.geometric.inf$region == "Indomalayan" &
                                                          predator.regional.r.beta.geometric.inf$consumer1 %in% comm.con.sp & 
                                                          predator.regional.r.beta.geometric.inf$consumer2 %in% comm.con.sp,] 
    hypervolume.beta.comm.inf <- rbind(hypervolume.beta.comm.inf, data.frame(Community = Tropic.carni.webs.inf95.10p[[i]]$properties$Community_unify,
                                                                             Region = Tropic.carni.webs.inf95.10p[[i]]$properties$Region,
                                                                             N.Sp.c.hv.beta.inf = length(comm.con.sp),
                                                                             resource.hv.beta.m.inf = round(mean(comm.beta$overlap_hypervolume),2)))
  } else {
    comm.con.sp <- unique(Tropic.carni.webs.inf95.10p[[i]]$trophic.links$consumer)
    comm.beta <-  predator.regional.r.beta.geometric.inf[predator.regional.r.beta.geometric.inf$region == "Neotropic" &
                                                           predator.regional.r.beta.geometric.inf$consumer1 %in% comm.con.sp & 
                                                           predator.regional.r.beta.geometric.inf$consumer2 %in% comm.con.sp,] 
    hypervolume.beta.comm.inf <- rbind(hypervolume.beta.comm.inf, data.frame(Community = Tropic.carni.webs.inf95.10p[[i]]$properties$Community_unify,
                                                                             Region = Tropic.carni.webs.inf95.10p[[i]]$properties$Region,
                                                                             N.Sp.c.hv.beta.inf = length(comm.con.sp),
                                                                             resource.hv.beta.m.inf = round(mean(comm.beta$overlap_hypervolume),2)))
  }
}

#Combine hypervolume to food web metrics 
food.web.dt.summary <- left_join(food.web.dt.summary, hypervolume.alpha.comm.inf[,-3], join_by(Community==Community, Region == Region))
food.web.dt.summary <- left_join(food.web.dt.summary, hypervolume.beta.comm.inf[,-3], join_by(Community==Community, Region == Region))
write.csv(food.web.dt.summary, "Data/Foodweb_metrics.csv", row.names = F)


#VII. Supplementary Figures related to food web metrics   --------
##1. Figure S6 --------
### Density distribution of proportion of the total predator species and basal prey species across the 389 mammal food webs by realm 
alpha.dgree = 0.9
basal.color <- "grey80";inter.color <- "grey50"; top.color <- "grey20" 
food.web.dt.summary$Per.predator.inf95.10p <- food.web.dt.summary$Per.inter.inf95.10p + food.web.dt.summary$Per.top.inf95.10p
food.web.dt.summary$Region <- factor(food.web.dt.summary$Region, levels = c("Neotropic", "Afrotropic", "Indomalayan"))

ggplot() +
  xlim(c(0.05,1)) + ylim(c(0,15)) + facet_wrap(~Region, ncol = 1) +
  geom_density(data = food.web.dt.summary, 
               aes(Per.basal.inf95.10p), fill = basal.color, color = "NA", alpha = alpha.dgree) +
  geom_density(data = food.web.dt.summary, 
               aes(Per.predator.inf95.10p), fill = top.color, color = "NA", alpha = alpha.dgree) + theme_bw() + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14),
        strip.text =  element_text(size = 14, face = "bold"), strip.background = element_blank(),
        legend.title = element_text(size = 16), legend.text = element_text(size = 14),
        panel.grid = element_blank()) +
  labs(y = "density", x = "Proportion of species richness in food web")


##2. Figure S8 --------
### Density distribution of proportion of the total predator species and basal prey species across the 389 mammal food webs by realm 
total.regional.predaator.dietary.breadth <- read.csv("Data/Regional_predator_dietary_hypervolume.csv")
colnames(total.regional.predaator.dietary.breadth) <- c("X", "region", "region.hypervolume", "region.n.species")

predator.dietary.breadth <- read.csv("Predator_dietary_niche_hypervolume.csv")
predator.dietary.breadth <- left_join(predator.dietary.breadth, mamm.func[,c("phylacine_binomial", "adult_mass_g")], join_by(consumer == phylacine_binomial))
predator.dietary.breadth <- left_join(predator.dietary.breadth, total.regional.predaator.dietary.breadth[,-1], join_by(region == region))

predator.dietary.breadth$region.hypervolume.log <- log10(predator.dietary.breadth$region.hypervolume)
predator.dietary.breadth$hypervolume.c.log <- log10(predator.dietary.breadth$hypervolume)
predator.dietary.breadth$adult_mass_g.log <- log10(predator.dietary.breadth$adult_mass_g)

predator.dietary.breadth$region <- factor(predator.dietary.breadth$region, levels = c("Neotropic", "Afrotropic", "Indomalayan"), labels = c("Neotropic", "Afrotropic", "Indomalaya"))
ggplot(predator.dietary.breadth, aes(n.species, hypervolume.c.log, fill = adult_mass_g.log)) + 
  facet_wrap(~region, scales = "free") +  
  ylim(2.5,5.45) +
  geom_abline(data = predator.dietary.breadth, aes(intercept = region.hypervolume.log, slope = 0), lty = 1, color = "red2", alpha = 0.8, linewidth = 0.5) + 
  geom_point(shape = 21, color = "NA", size = 2.5, alpha = 0.8) + scale_fill_viridis_c(option = "viridis") + theme_bw() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14, color = "grey20"),
        legend.text = element_text(size = 12), legend.title = element_text(size = 14, color = "grey20"),
        panel.grid = element_blank(), strip.background = element_blank(), strip.text = element_text(size = 16, color = "black")) +
  labs(x = "Number of prey species", y = expression(log[10]*"(dietary niche breadth)"), fill = expression(log[10]*"(body mass [g])"))

