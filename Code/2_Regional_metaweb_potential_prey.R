#
# R script to identify potential predator-prey pairs in regional meta webs, with
#  potential prey of a focal predator identified by its trait similarity with documented prey (trait-matching prey)
#  
# Input: 
#  Documented predator-prey pairs and community checklist >= 500g of the 408 communities that consisted of predator-prey interactions
#     1. Region_sp_intx_500g.RData 
#
#  Trait and Taxonomy from Soria et al. (2021) COMBINE dataset https://doi.org/10.1002/ecy.3344 to have fine-scale trait information to characterize ecological similarity
#     1. trait_data_imputed.csv
#
#
# Output: cosine similarity of predator-prey pairs of documented prey and potential prey in each regional meta-web; 
#          potential prey based on criteria of cosine similarity >= 0.95 (in95) and 10% body mass differences (10p) relative to documented prey, 
#    1. Regional_carni_links_500g_cosine_similarity_in95_10p_all.csv
#  


#Call documented predator-prey interactions in regional meta webs and species list across 408 communities
load("Data/Region_sp_intx_500g.RData")

#I. Collect species ecological traits to quantify ecological similarity between species pairs ==============
## ecological traits were obtained from Soria et al. 2011  https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecy.3344
## ecological trait data (trait_data_imputed.csv) deposit location: 
#    https://figshare.com/articles/dataset/COMBINE_a_Coalesced_Mammal_Database_of_Intrinsic_and_extrinsic_traits/13028255
mamm.func <- read.csv("./trait_data_imputed.csv")

mamm.func.analyzed <- mamm.func[mamm.func$phylacine_binomial %in% df_region.match$Species,]
mamm.func.analyzed <- mamm.func.analyzed[!duplicated(mamm.func.analyzed$phylacine_binomial),]
## fill values from IUCN habitat types 
mamm.func.analyzed[mamm.func.analyzed$phylacine_binomial == "Cebus aequatorialis",]$habitat_breadth_n <- 1 #https://www.iucnredlist.org/species/4081/191702052#habitat-ecology
mamm.func.analyzed[mamm.func.analyzed$phylacine_binomial == "Cebus albifrons",]$habitat_breadth_n <- 2 #https://www.iucnredlist.org/species/39951/191703935#habitat-ecology
mamm.func.analyzed[mamm.func.analyzed$phylacine_binomial == "Cebus brunneus",]$habitat_breadth_n <- 1 #https://www.iucnredlist.org/species/81237954/17981252#habitat-ecology
mamm.func.analyzed[mamm.func.analyzed$phylacine_binomial == "Cebus capucinus",]$habitat_breadth_n <- 2 #https://www.iucnredlist.org/species/81257277/191708164#habitat-ecology
mamm.func.analyzed[mamm.func.analyzed$phylacine_binomial == "Hylobates muelleri",]$habitat_breadth_n <- 1 #https://www.iucnredlist.org/species/39888/17990934#habitat-ecology
mamm.func.analyzed[mamm.func.analyzed$phylacine_binomial == "Loxodonta africana",]$habitat_breadth_n <- 6 #https://www.iucnredlist.org/species/181008073/223031019#habitat-ecology
mamm.func.analyzed[mamm.func.analyzed$phylacine_binomial == "Nomascus gabriellae",]$habitat_breadth_n <- 1 #https://www.iucnredlist.org/species/128073282/17968950#habitat-ecology
mamm.func.analyzed[mamm.func.analyzed$phylacine_binomial == "Trachypithecus barbei",]$habitat_breadth_n <- 1 #https://www.iucnredlist.org/species/41554/17960144#habitat-ecology

nrow(mamm.func.analyzed) #686 species >= 500g have sufficient trait | Total 741 species in 408 communities

##1. Collect body mass to calcuate body mass ratio of prey pairs to identify potential prey relative to documeted prey ------
mamm.func$adult_kg <- round(mamm.func$adult_mass_g/1000,2)
##prepare body mass information for consumer 
mamm.func.c <- mamm.func[,c("phylacine_binomial", "order", "adult_kg")]; colnames(mamm.func.c) <- c("phylacine_binomial", "order_c", "adult_kg_c")
##prepare body mass information for documented resource & potential resource 
mamm.func.r <- mamm.func[,c("phylacine_binomial", "order", "adult_kg")]; colnames(mamm.func.r) <- c("phylacine_binomial", "order_r", "adult_kg_r")
mamm.func.p.r <- mamm.func[,c("phylacine_binomial", "order", "adult_kg")]; colnames(mamm.func.p.r) <- c("phylacine_binomial", "order_p_r", "adult_kg_p_r")

##2. Transform categorical trait for cosine-similarity -----------------------
mamm.func.analyzed$fossorialityG <- 0
mamm.func.analyzed[mamm.func.analyzed$fossoriality == 1,]$fossorialityG <- 1

mamm.func.analyzed$forage_G <- 0
mamm.func.analyzed$forage_S <- 0
mamm.func.analyzed$forage_Ar <- 0
mamm.func.analyzed$forage_A <- 0
mamm.func.analyzed[mamm.func.analyzed$foraging_stratum == "G",]$forage_G <- 1
mamm.func.analyzed[mamm.func.analyzed$foraging_stratum == "S",]$forage_S <- 1
mamm.func.analyzed[mamm.func.analyzed$foraging_stratum == "Ar",]$forage_Ar <- 1

mamm.func.analyzed$activity_noct <- 0
mamm.func.analyzed$activity_diur <- 0
mamm.func.analyzed$activity_inter <- 0
mamm.func.analyzed[mamm.func.analyzed$activity_cycle == 1,]$activity_noct <- 1
mamm.func.analyzed[mamm.func.analyzed$activity_cycle == 3,]$activity_diur <- 1
mamm.func.analyzed[mamm.func.analyzed$activity_cycle == 2,]$activity_inter <- 1

mamm.func.analyzed$trophic_level_1 <- 0
mamm.func.analyzed$trophic_level_2 <- 0
mamm.func.analyzed$trophic_level_3 <- 0
mamm.func.analyzed[mamm.func.analyzed$trophic_level == 1,]$trophic_level_1 <- 1
mamm.func.analyzed[mamm.func.analyzed$trophic_level == 3,]$trophic_level_2 <- 1
mamm.func.analyzed[mamm.func.analyzed$trophic_level == 2,]$trophic_level_3 <- 1

mamm.func.analyzed.feature <- mamm.func.analyzed[,c("phylacine_binomial","adult_mass_g" ,"brain_mass_g","adult_body_length_mm",    
                                                    "max_longevity_d","female_maturity_d","age_first_reproduction_d", 
                                                    "gestation_length_d","litter_size_n","litters_per_year_n","interbirth_interval_d" ,
                                                    "weaning_age_d" ,"generation_length_d","dispersal_km","hibernation_torpor","fossorialityG" ,
                                                    "dphy_invertebrate", "dphy_vertebrate","dphy_plant", "det_diet_breadth_n", "habitat_breadth_n", 
                                                    "forage_G","forage_S","forage_Ar" ,"forage_A","activity_noct","activity_diur" ,"activity_inter" ,
                                                    "trophic_level_1","trophic_level_2","trophic_level_3")]

df_sp_feature.match.cat.trans <- mamm.func.analyzed.feature[,-1]
rownames(df_sp_feature.match.cat.trans) <- mamm.func.analyzed.feature$phylacine_binomial


#II. Calculate cosine similarity between  predator-i - documented prey-j and predator-i - species-k in regional meta web =====================
## cosine similarity values range from 0 to 1. Self pairs or prey pairs with idential traits will have similarity values = 1
## Call function for cosine similarity
library(lsa)
get_predicted_resource_by_diff <- function(region, n=10, threshold=0, similarity="cosine") #based on vector in multi-trait space
{
  #Trait data intput directory
  load("Data/Region_sp_intx_500g.RData")
  
  #define region trait space 
  df_region_intx.match$Region <- factor(df_region_intx.match$Region, levels = c("Neotropic", "Afrotropic", "Indomalayan"),
                                        labels = c("neo", "afr", "ind"))
  # get regional intx (species pair would like to calculate cosine similarity)
  df_region_sp <- df_region.match[df_region.match$Region == region,]
  df_community_intx <- df_region_intx.match[df_region_intx.match$Region %in% region, ]
  
  # regional species: consumer, resource, and isolated species in regional meta web
  region_consumers <- unique(df_community_intx$consumer)
  region_resources <- unique(df_community_intx$resource)
  region_sp <- unique(df_region_sp$Species)
  isolated_sp <- setdiff(region_sp, unique(c(region_consumers, region_resources)))
  
  # Normalization on whole scale
  # continuous variables
  continuous_var = c("adult_mass_g", "brain_mass_g", "adult_body_length_mm", "max_longevity_d", 
                     "female_maturity_d", "age_first_reproduction_d", "gestation_length_d", "litter_size_n", 
                     "litters_per_year_n", "interbirth_interval_d", "weaning_age_d", "generation_length_d", 
                     "dispersal_km", "det_diet_breadth_n", "habitat_breadth_n")
  # seperate the continuous and rest of variables
  # contiuous
  df_continuous_feature <- df_sp_feature.match.cat.trans[, continuous_var]
  df_continuous_feature.scale <- scale(df_continuous_feature)
  df_continuous_feature.scale <- df_continuous_feature.scale[region_sp, ]
  
  # discret
  df_rest_feature <- df_sp_feature.match.cat.trans[, !(colnames(df_sp_feature.match.cat.trans) %in% continuous_var)]
  percentage_features <- c("dphy_invertebrate", "dphy_vertebrate", "dphy_plant")
  df_rest_feature[, percentage_features] <- df_rest_feature[, percentage_features] / 100
  df_rest_feature <- df_rest_feature[region_sp, ]
  
  # consumer diff matrix
  df_cr_intx <- data.frame()
  for (consumer in region_consumers)
  {
    df_consumer_resource <- df_community_intx[df_community_intx$consumer == consumer,]
    sp_resources <- unique(df_consumer_resource$resource)
    
    # normalize the continous variables
    mat_continuous_feature_norm <- apply(as.matrix(df_continuous_feature.scale), 1, function(x) {
      
      return((as.matrix(df_continuous_feature.scale[consumer, ] - x)))
    })
    df_continuous_feature_norm <- data.frame(t(mat_continuous_feature_norm))
    colnames(df_continuous_feature_norm) <- colnames(df_continuous_feature.scale)
    
    # substract the rest of variables
    mat_rest_feature_norm <- apply(as.matrix(df_rest_feature), 1, function(x) {
      return(as.matrix(df_rest_feature[consumer, ]) - x)
    })
    df_rest_feature_norm <- data.frame(t(mat_rest_feature_norm))
    colnames(df_rest_feature_norm) <- colnames(df_rest_feature)
    df_rest_feature_norm <- df_rest_feature_norm[rownames(df_continuous_feature_norm), ]
    
    # combine continous and rest variables
    df_consumer_feature_norm <- cbind(df_continuous_feature_norm, df_rest_feature_norm)
    df_consumer_feature_norm <- df_consumer_feature_norm[!(rownames(df_consumer_feature_norm) == consumer), ]
    
    # compute the feature similarity
    if(similarity == "gower")
    {
      mat_region_sp_simi <- 1 - as.matrix(daisy(as.matrix(df_consumer_feature_norm), metric="gower", stand=FALSE))
    } else {
      mat_region_sp_simi <- cosine(t(as.matrix(df_consumer_feature_norm)))
    }
    df_region_sp_simi <- as.data.frame(mat_region_sp_simi)
    colnames(df_region_sp_simi) <- rownames(df_consumer_feature_norm)
    rownames(df_region_sp_simi) <- rownames(df_consumer_feature_norm)
    
    
    for (sp_resource in sp_resources)
    {
      df_sp_resource_simi <- df_region_sp_simi[order(df_region_sp_simi[, sp_resource], decreasing=TRUE), ]
      
      # choose top n species with the highest similarity then filter by threshold
      #df_sp_resource_simi <- df_sp_resource_simi[1:(n+1), ]
      df_sp_resource_simi <- df_sp_resource_simi[1:n, ]
      df_sp_resource_simi <- df_sp_resource_simi[df_sp_resource_simi[, sp_resource] >= threshold, ]
      
      # if nothing left, then next sp_resource
      if(nrow(df_sp_resource_simi) == 0) {next}
      
      # documented
      documented <- ifelse(rownames(df_sp_resource_simi) %in% sp_resources, 1, 0)
      
      # save results
      df_cr_intx <- rbind(df_cr_intx, data.frame("consumer"=consumer,
                                                 "resource"=sp_resource,
                                                 "potential_resources"=rownames(df_sp_resource_simi),
                                                 "similarity"=df_sp_resource_simi[,sp_resource],
                                                 "documented"=documented,
                                                 check.names=FALSE))
    }
  }
  
  return (df_cr_intx)
}

##1. Neotropic meta web with the body mass ratio ---------------------
inferred.links.neo <- get_predicted_resource_by_diff("neo", n = (length(unique(df_region.match[df_region.match$Region == "neo",]$Species))-2), 
                                                     threshold=0, similarity="cosine")
inferred.links.neo.links <- left_join(inferred.links.neo, mamm.func.c, join_by(consumer == phylacine_binomial))
inferred.links.neo.links <- left_join(inferred.links.neo.links, mamm.func.r, join_by(resource == phylacine_binomial))
inferred.links.neo.links <- left_join(inferred.links.neo.links, mamm.func.p.r, join_by(potential_resources == phylacine_binomial))
inferred.links.neo.links$body.ratio.prey.pairs <- round(inferred.links.neo.links$adult_kg_p_r/inferred.links.neo.links$adult_kg_r,2)
inferred.links.neo.links$region <- "neo"
inferred.links.neo.links.filtered <- inferred.links.neo.links[inferred.links.neo.links$similarity  >= 0.95 & 
                                                                inferred.links.neo.links$body.ratio.prey.pairs >= 0.9 & inferred.links.neo.links$body.ratio.prey.pairs <= 1.1,]
length(unique(inferred.links.neo.links.filtered$potential_resources)) #113 trait-matched prey species, including both documented and potential prey

##2. Afrotropic meta web with the body mass ratio ---------------------
inferred.links.afr <- get_predicted_resource_by_diff("afr", n = (length(unique(df_region.match[df_region.match$Region == "afr",]$Species))-2), 
                                                     threshold=0, similarity="cosine")
inferred.links.afr.links <- left_join(inferred.links.afr, mamm.func.c, join_by(consumer == phylacine_binomial))
inferred.links.afr.links <- left_join(inferred.links.afr.links, mamm.func.r, join_by(resource == phylacine_binomial))
inferred.links.afr.links <- left_join(inferred.links.afr.links, mamm.func.p.r, join_by(potential_resources == phylacine_binomial))
inferred.links.afr.links$body.ratio.prey.pairs <- round(inferred.links.afr.links$adult_kg_p_r/inferred.links.afr.links$adult_kg_r,2)
inferred.links.afr.links$region <- "afr"
##subset prey pairs with cosine similarity >= 0.95 & body mass ratio 10% 
inferred.links.afr.links.filtered <- inferred.links.afr.links[inferred.links.afr.links$similarity >= 0.95 & 
                                                                inferred.links.afr.links$body.ratio.prey.pairs >= 0.9 & inferred.links.afr.links$body.ratio.prey.pairs <= 1.1,]
length(unique(inferred.links.afr.links.filtered$potential_resources)) #204 trait-matched prey species, including both documented and potential prey 

##3. Indomaylayan meta web with the body mass ratio ---------------------
inferred.links.ind <- get_predicted_resource_by_diff("ind", n = (length(unique(df_region.match[df_region.match$Region == "ind",]$Species))-2), 
                                                     threshold=0, similarity="cosine")

inferred.links.ind.links <- left_join(inferred.links.ind, mamm.func.c, join_by(consumer == phylacine_binomial))
inferred.links.ind.links <- left_join(inferred.links.ind.links, mamm.func.r, join_by(resource == phylacine_binomial))
inferred.links.ind.links <- left_join(inferred.links.ind.links, mamm.func.p.r, join_by(potential_resources == phylacine_binomial))
inferred.links.ind.links$body.ratio.prey.pairs <- round(inferred.links.ind.links$adult_kg_p_r/inferred.links.ind.links$adult_kg_r,2)
inferred.links.ind.links$region <- "ind"
inferred.links.ind.links.filtered <- inferred.links.ind.links[inferred.links.ind.links$similarity  >= 0.95 & 
                                                                inferred.links.ind.links$body.ratio.prey.pairs >= 0.9 & inferred.links.ind.links$body.ratio.prey.pairs <= 1.1,]
length(unique(inferred.links.ind.links$potential_resources)) #227 trait-matched prey species, including both documented and potential prey

##4. Combine all links (raw and filtered by similarity and body mass ratio, not unique) =======
inferred.region.links <- rbind(inferred.links.afr.links, inferred.links.neo.links, inferred.links.ind.links) 
nrow(inferred.region.links); length(unique(inferred.region.links$potential_resources))

inferred.region.links.filtered <- rbind(inferred.links.afr.links.filtered, inferred.links.neo.links.filtered, inferred.links.ind.links.filtered) 
nrow(inferred.region.links.filtered); length(unique(inferred.region.links.filtered$potential_resources)) #total 432 trait-matched prey species based on >= 0.95 cosince sim +- 10% body mass

write.csv(inferred.region.links[,c("region", "consumer", "order_c", "resource", "order_r", "potential_resources" , "order_p_r", 
                                   "similarity","documented","body.ratio.prey.pairs")], 
          "Regional_carni_links_500g_cosine_similarity_in95_10p_all.csv", row.names = F)


#Supp: Evaluate cosine similarity threshold (Figure S18) ------------------
library(ggplot2); library(ggpubr)
length(unique(inferred.region.links$potential_resources)) #total 686 species fed in trait space
inferred.region.links.filtered <- inferred.region.links[ inferred.region.links$body.ratio.prey >= 0.9 & inferred.region.links$body.ratio.prey <= 1.1,]
inferred.region.links.filtered$similarity.round <- round(inferred.region.links.filtered$similarity,2)
inferred.region.links.filtered$region <- factor(inferred.region.links.filtered$region, levels = c("neo", "afr", "ind"),
                                         labels = c("Neotropic", "Afrotropic", "Indomalayan"))
cosine.similarity.threshold <- seq(0.5, 1, by = 0.01)

cosine.sim <- c()
for (i in 1:length(cosine.similarity.threshold)) {
  cosine.sim <- rbind(cosine.sim,
                      data.frame(inferred.region.links.filtered[ inferred.region.links.filtered$similarity.round >= cosine.similarity.threshold[i], ] %>% 
                                   group_by(consumer, region) %>% summarise(N.links = length(unique(potential_resources))),
                                 cosine.sim = cosine.similarity.threshold[i]))
}
cosine.sim$consumer <- factor(cosine.sim$consumer)
cosine.sim %>% group_by(region) %>% summarise(N.c.sp = length(unique(consumer)))
cosine.sim <- left_join(cosine.sim, unique(inferred.region.links[,c("consumer", "adult_kg_c")]), join_by(consumer == consumer))
cosine.sim$adult_g_c.log <- log10(cosine.sim$adult_kg_c*1000)
hist(cosine.sim$adult_g_c.log )

realm.sp.pool <- df_region.match %>% group_by(Region) %>% summarise(realm.N.sp = length(unique(Species)))
realm.sp.pool$Region <- factor(realm.sp.pool$Region, levels = c("neo", "afr", "ind"), labels = c("Neotropic", "Afrotropic", "Indomalayan"))
cosine.sim <- left_join(cosine.sim, realm.sp.pool, join_by(region == Region))
cosine.sim$region <- factor(cosine.sim$region, levels = c("Neotropic", "Afrotropic", "Indomalayan"))

ggplot(cosine.sim, aes(cosine.sim, N.links, color = adult_g_c.log, group = adult_g_c.log)) + 
  facet_wrap(~region, ncol = 3) + geom_point() + scale_color_viridis_c() + theme_bw() + ylim(0,280) +
  geom_vline(xintercept = 0.95, lty = 2, color = "grey60") + 
  geom_abline(data = cosine.sim, aes(intercept = realm.N.sp, slope = 0), lty = 1, color = "red3") + 
  theme(strip.text = element_text(size = 14), strip.background = element_blank(), axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
  labs(x = "Cosine similarity", y = "Number of prey species", color = expression(atop(log[10]*"(body", paste("mass [g])"))))


#Supp: Evaluate cosine similarity vs. functional distance and phylogenetic distance (Figures S19-20) ------------------
mamm.func <- read.csv("./trait_data_imputed.csv")
mamm.func.analyzed <- mamm.func[mamm.func$phylacine_binomial %in% df_region.match$Species,]
mamm.func.analyzed <- mamm.func.analyzed[!duplicated(mamm.func.analyzed$phylacine_binomial),]
nrow(mamm.func.analyzed) #686 species
mamm.func.analyzed <- mamm.func.analyzed[is.na(mamm.func.analyzed$adult_mass_g) == F,]
mamm.func.analyzed[is.na(mamm.func.analyzed$brain_mass_g) == T,]
mamm.func.analyzed[is.na(mamm.func.analyzed$adult_body_length_mm) == T,]
mamm.func.analyzed[is.na(mamm.func.analyzed$max_longevity_d) == T,]
mamm.func.analyzed[is.na(mamm.func.analyzed$female_maturity_d) == T,]
mamm.func.analyzed[is.na(mamm.func.analyzed$age_first_reproduction_d) == T,]
mamm.func.analyzed[is.na(mamm.func.analyzed$gestation_length_d) == T,]
mamm.func.analyzed[is.na(mamm.func.analyzed$litter_size_n) == T,]
mamm.func.analyzed[is.na(mamm.func.analyzed$litters_per_year_n) == T,]
mamm.func.analyzed[is.na(mamm.func.analyzed$interbirth_interval_d) == T,]
mamm.func.analyzed[is.na(mamm.func.analyzed$weaning_age_d) == T,]
mamm.func.analyzed[is.na(mamm.func.analyzed$generation_length_d) == T,]
mamm.func.analyzed[is.na(mamm.func.analyzed$dispersal_km) == T,]
mamm.func.analyzed[is.na(mamm.func.analyzed$hibernation_torpor) == T,]
mamm.func.analyzed[is.na(mamm.func.analyzed$fossoriality) == T,]
mamm.func.analyzed[is.na(mamm.func.analyzed$dphy_vertebrate) == T,] 
mamm.func.analyzed[is.na(mamm.func.analyzed$det_diet_breadth_n) == T,] 
#8 species without habitat breadth information: Fill values from IUCN habitat types  
mamm.func.analyzed[is.na(mamm.func.analyzed$habitat_breadth_n) == T,]
mamm.func.analyzed[is.na(mamm.func.analyzed$habitat_breadth_n) == T,]$phylacine_binomial
#Fill values from IUCN habitat types 
mamm.func.analyzed[mamm.func.analyzed$phylacine_binomial == "Cebus aequatorialis",]$habitat_breadth_n <- 1 #https://www.iucnredlist.org/species/4081/191702052#habitat-ecology
mamm.func.analyzed[mamm.func.analyzed$phylacine_binomial == "Cebus albifrons",]$habitat_breadth_n <- 2 #https://www.iucnredlist.org/species/39951/191703935#habitat-ecology
mamm.func.analyzed[mamm.func.analyzed$phylacine_binomial == "Cebus brunneus",]$habitat_breadth_n <- 1 #https://www.iucnredlist.org/species/81237954/17981252#habitat-ecology
mamm.func.analyzed[mamm.func.analyzed$phylacine_binomial == "Cebus capucinus",]$habitat_breadth_n <- 2 #https://www.iucnredlist.org/species/81257277/191708164#habitat-ecology
mamm.func.analyzed[mamm.func.analyzed$phylacine_binomial == "Hylobates muelleri",]$habitat_breadth_n <- 1 #https://www.iucnredlist.org/species/39888/17990934#habitat-ecology
mamm.func.analyzed[mamm.func.analyzed$phylacine_binomial == "Loxodonta africana",]$habitat_breadth_n <- 6 #https://www.iucnredlist.org/species/181008073/223031019#habitat-ecology
mamm.func.analyzed[mamm.func.analyzed$phylacine_binomial == "Nomascus gabriellae",]$habitat_breadth_n <- 1 #https://www.iucnredlist.org/species/128073282/17968950#habitat-ecology
mamm.func.analyzed[mamm.func.analyzed$phylacine_binomial == "Trachypithecus barbei",]$habitat_breadth_n <- 1 #https://www.iucnredlist.org/species/41554/17960144#habitat-ecology
nrow(mamm.func.analyzed) #686 species >= 500g have sufficient trait 

#subset trait dataframe with required traits
mamm.func.analyzed.for.gower <- mamm.func.analyzed[,c("phylacine_binomial","adult_mass_g" ,"brain_mass_g","adult_body_length_mm",    
                                                      "max_longevity_d","female_maturity_d","age_first_reproduction_d", 
                                                      "gestation_length_d","litter_size_n","litters_per_year_n","interbirth_interval_d" ,
                                                      "weaning_age_d" ,"generation_length_d","dispersal_km","hibernation_torpor","fossoriality" ,
                                                      "dphy_invertebrate", "dphy_vertebrate","dphy_plant", "det_diet_breadth_n", "habitat_breadth_n", 
                                                      "foraging_stratum", "activity_cycle", "trophic_level")]
#Organize into trait matrix 
rownames(mamm.func.analyzed.for.gower) <- mamm.func.analyzed.for.gower$phylacine_binomial
df_region.match.for.gower <- df_region.match[df_region.match$Species %in% mamm.func.analyzed.for.gower$phylacine_binomial, ]
df_sp_feature.match.for.gower <- mamm.func.analyzed.for.gower[,-1]
df_region_intx.match.for.gower <- df_region_intx.match[df_region_intx.match$resource %in% mamm.func.analyzed.for.gower$phylacine_binomial &
                                                         df_region_intx.match$consumer %in% mamm.func.analyzed.for.gower$phylacine_binomial, ]

length(unique(df_region.match.for.gower$Species)) #686 species 
length(unique(c(df_region_intx.match.for.gower$resource, df_region_intx.match.for.gower$consumer)))  #350 linked species

##1. Functional distance based on Gower --------------
sp_feature_gower <- df_sp_feature.match.for.gower
rownames(sp_feature_gower) <- gsub(" ", "_", rownames(sp_feature_gower))

continuous_var = c("adult_mass_g", "brain_mass_g", "adult_body_length_mm", "max_longevity_d", 
                   "female_maturity_d", "age_first_reproduction_d", "gestation_length_d", "litter_size_n", 
                   "litters_per_year_n", "interbirth_interval_d", "weaning_age_d", "generation_length_d", "dispersal_km",
                   "det_diet_breadth_n", "habitat_breadth_n", 
                   "hibernation_torpor", "fossoriality")
cat_var <- c( "foraging_stratum" ,"activity_cycle" ,"trophic_level")
fuzzy_var <- c("dphy_invertebrate", "dphy_vertebrate", "dphy_plant")

df_cat_feature <- sp_feature_gower[, cat_var]
df_cat_feature <- df_cat_feature %>% 
  mutate(across(where(is.numeric), as.factor)) %>% 
  mutate(across(where(is.character), as.factor))

df_fuzzy_feature <- sp_feature_gower[, fuzzy_var]
df_fuzzy_feature <- df_fuzzy_feature %>% 
  mutate(across(where(is.numeric), as.integer)) 

sp_feature_gower_processed <- cbind(sp_feature_gower[,continuous_var], df_cat_feature, df_fuzzy_feature)
trait_type <- data.frame(trait_name = colnames(sp_feature_gower_processed),
                         trait_type =c(rep("Q",17), rep("N",3),rep("F",3)),
                         fuzzy_name = c(rep(NA,20),c(rep("dphy_diet", 3))))

sp_feature_gower.dist <- as.matrix(funct.dist(sp_feature_gower_processed, trait_type, 
                                              metric        = "gower",
                                              scale_euclid  = "scale_center",
                                              ordinal_var   = "classic",
                                              weight_type   = "equal",
                                              stop_if_NA    = TRUE))
dim(sp_feature_gower.dist)

##2. Phylogenetic distance  --------------
mammal.tree <- read.nexus("./Complete_phylogeny.nex") #1-k credible set downloaded from http://vertlife.org/data/birds/
observed.pruned.tree <- ape::drop.tip(mammal.tree[[1]], mammal.tree[[1]]$tip.label[!mammal.tree[[1]]$tip.label %in% rownames(sp_feature_gower)]) 
sp_phylo.dist <- cophenetic.phylo(observed.pruned.tree)
dim(sp_phylo.dist)

##3.  Phylogenetic-Functional distance --------------
mamm.func <- read.csv("./trait_data_imputed.csv")
mamm.func <- mamm.func[!duplicated(mamm.func$phylacine_binomial),]
mamm.func$adult_kg <- round(mamm.func$adult_mass_g/1000,2)
mamm.func.r <- mamm.func[,c("phylacine_binomial", "order", "adult_kg")]; colnames(mamm.func.r) <- c("phylacine_binomial", "order_r", "adult_kg_r")
mamm.func.p.r <- mamm.func[,c("phylacine_binomial", "order", "adult_kg")]; colnames(mamm.func.p.r) <- c("phylacine_binomial", "order_p_r", "adult_kg_p_r")

###FD
FDist.global.pair <- setNames(melt(sp_feature_gower.dist), c('Species', 'Compare_species', 'FDist'))
FDist.global.pair$Species <- gsub("_", " ", FDist.global.pair$Species)
FDist.global.pair$Compare_species <- gsub("_", " ", FDist.global.pair$Compare_species)
FDist.global.pair$Similarity_FD <- 1- FDist.global.pair$FDist

FDist.global.pair <- left_join(FDist.global.pair, mamm.func.r, join_by(Species ==phylacine_binomial))
FDist.global.pair <- left_join(FDist.global.pair, mamm.func.p.r, join_by(Compare_species ==phylacine_binomial))
FDist.global.pair$body.ratio.prey <- round(FDist.global.pair$adult_kg_p_r/FDist.global.pair$adult_kg_r,2)
FDist.global.pair.filtered <- FDist.global.pair[FDist.global.pair$body.ratio.prey >= 0.9 & FDist.global.pair$body.ratio.prey <= 1.1,]
FDist.global.pair.filtered$prey.pairs <- paste0(FDist.global.pair.filtered$Species, ".", FDist.global.pair.filtered$Compare_species)

###PFD
FPDist <- function(PDist, FDist, a, p, ord = TRUE){
  PDist <- as.matrix(PDist, rownames.force = TRUE)
  FDist <- as.matrix(FDist, rownames.force = TRUE)
  if(
    is.null(rownames(FDist)) ||
    is.null(colnames(FDist)) ||
    is.null(rownames(PDist)) ||
    is.null(colnames(PDist))
  ) stop('distance matrices must have row and column names')
  if(
    any(rownames(FDist) != colnames(FDist)) ||
    any(rownames(PDist) != colnames(PDist))
  ) stop('row and column names must match for distance matrices')
  if(ord){
    FDist <- FDist[order(rownames(FDist)), order(rownames(FDist))]
    PDist <- PDist[order(rownames(PDist)), order(rownames(PDist))]
  }
  if(
    any(rownames(FDist) != rownames(PDist)) ||
    any(colnames(FDist) != colnames(PDist))
  ) stop('FDist and PDist must have same row and column names')
  
  FDist <- FDist/max(FDist)
  PDist <- PDist/max(PDist)
  
  ((a*(PDist^p)) + ((1-a)*(FDist^p)))^(1/p)
}


rownames(sp_phylo.dist) == rownames(sp_feature_gower.dist)
FPDist.global=FPDist(PDist=sp_phylo.dist,FDist=sp_feature_gower.dist,a=0.5, p=2)
FPDist.global.pair <- setNames(melt(FPDist.global), c('Species', 'Compare_species', 'PFDist'))
FPDist.global.pair$Species <- gsub("_", " ", FPDist.global.pair$Species)
FPDist.global.pair$Compare_species <- gsub("_", " ", FPDist.global.pair$Compare_species)
FPDist.global.pair$Similarity_FPD <- 1- FPDist.global.pair$PFDist

FPDist.global.pair <- left_join(FPDist.global.pair, mamm.func.r, join_by(Species ==phylacine_binomial))
FPDist.global.pair <- left_join(FPDist.global.pair, mamm.func.p.r, join_by(Compare_species ==phylacine_binomial))
FPDist.global.pair$body.ratio.prey <- round(FPDist.global.pair$adult_kg_p_r/FPDist.global.pair$adult_kg_r,2)
FPDist.global.pair.filtered <- FPDist.global.pair[FPDist.global.pair$body.ratio.prey >= 0.9 & FPDist.global.pair$body.ratio.prey <= 1.1,]
FPDist.global.pair.filtered$prey.pairs <- paste0(FPDist.global.pair.filtered$Species, ".", FPDist.global.pair.filtered$Compare_species)

FD.FPD.compare <- left_join(FDist.global.pair.filtered[,c("prey.pairs", "Similarity_FD")],
                            FPDist.global.pair.filtered[,c("prey.pairs", "Similarity_FPD")], join_by(prey.pairs == prey.pairs))

ggplot(FD.FPD.compare, aes(Similarity_FD, Similarity_FPD)) + 
  geom_point() + geom_abline(slope = 1) + theme_bw() + 
  labs(x = "Functional similarity", y = "Functional-phylogenetic similarity")

##4. Cosine similarity vs. FD ------------
cosine.sim.raw.filtered <- inferred.region.links[inferred.region.links$similarity >= 0.5 & inferred.region.links$body.ratio.prey >= 0.9 & inferred.region.links$body.ratio.prey <= 1.1,]
cosine.sim.raw.filtered$prey.pairs <- paste0(cosine.sim.raw.filtered$resource, ".", cosine.sim.raw.filtered$potential_resources)

similarity.compare <- left_join(cosine.sim.raw.filtered,FD.FPD.compare, join_by(prey.pairs == prey.pairs))
similarity.compare$region <- factor(similarity.compare$region, levels= c("neo", "afr", "ind"), labels = c("Neotropic", "Afrotropic", "Indomalayan"))

Region.color <- c("Neotropic" = "#adc9bc","Afrotropic" = "#ffb453" , "Indomalayan" = "#e6aaf0")
ggplot(similarity.compare, aes(similarity, Similarity_FD, color = region)) + 
  geom_vline(xintercept = 0.95, lty = 2, color = "grey40") + geom_hline(yintercept = 0.95, lty = 2, color = "grey40") + 
  geom_point(alpha = 0.6) + geom_abline(slope = 1) + theme_bw() +  
  scale_color_manual(values = Region.color) + 
  theme(panel.spacing=unit(1,"lines"), strip.text = element_text(size = 14), axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
  labs(x = "Cosine similarity", y = "Functional similarity", color = "Realm")

