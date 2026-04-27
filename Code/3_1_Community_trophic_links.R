#
# R script to augment local food webs with identified potential prey based on cosine similarity approach 
#
# 
# Input: 
#   1. Region_sp_intx_500g.RData
#   2. Rowan_Comm_10kmBuff_Env_cov.csv
#   3. Regional_carni_links_500g_cosine_similarity_in95_10p_all.csv
#
# Output:
#  Community408_carni_web_doc_inf95_500g_SR3_commcollect.Rdata
#   community collection of mammal food webs with documented trophic links and those augmented with potential trait-matched prey () for subsequent food web structure measures


library(dplyr); library(ggplot2); library(ggpubr)
library(cheddar)

load("Data/Region_sp_intx_500g.RData")

#I. Combine documented Global meta-web and Afrotropical meta-web as total documented local webs =======
## Call community checklists with species >= 500g and sufficient trait values as potential prey
community.sp.list.408comm.500g <- community.sp.list.408comm.500g[community.sp.list.408comm.500g$Species %in% df_region.match$Species,]
length(unique(community.sp.list.408comm.500g$Species)) #686 species

Afro.extra.doc.carni.links <- df_region_intx.match[df_region_intx.match$Global.web == 0 & df_region_intx.match$African.web == 1,]
Afro.extra.doc.carni.links$intx
unique(Afro.extra.doc.carni.links$c.order)

comm.carni.web.doc.link <- list()
extra.documented.links <- c()
for (i in 1:length(comm.web)){
  print(i)
  #Subset carnivoran predator prey links in each web
  intx <- comm.web[[i]]
  carni.intx <- intx[intx$c.order == "Carnivora",]
  #Subset community checklist
  comm.sp.list <- community.sp.list.408comm.500g[community.sp.list.408comm.500g$Community == names(comm.web[i]),]$Species
  #Subset extra documented carnivora predator and links in Afrotropical meta-web by community checklist
  comm.sp.with.extra.carni.link <- Afro.extra.doc.carni.links[Afro.extra.doc.carni.links$consumer %in% comm.sp.list,]
  comm.sp.with.extra.carni.link <- comm.sp.with.extra.carni.link[comm.sp.with.extra.carni.link$resource %in%  comm.sp.list,]
  
  if (nrow(comm.sp.with.extra.carni.link)==0) {
    comm.carni.web.doc.link[[i]] <- carni.intx[,colnames(carni.intx) %in% c("resource", "consumer", "c.order")]
    
  } else {
    comm.carni.web.doc.link[[i]] <- rbind(carni.intx[,colnames(carni.intx) %in% c("resource", "consumer", "c.order")],
                                         comm.sp.with.extra.carni.link[,colnames(comm.sp.with.extra.carni.link) %in% c("resource", "consumer", "c.order")])    
    extra.documented.links <- rbind(extra.documented.links, data.frame(Community = names(comm.web[i]), comm.sp.with.extra.carni.link))
  }
}
names(comm.carni.web.doc.link) <- names(comm.web)
#save(comm.carni.web.doc.link, file = "community408_web_carni_documented_500g_with_afro_links.RData")
#write.csv(extra.documented.links, "regional_408comm_doc_carni_links_500g_extra_afro_carni_links.csv", row.names = F)

#II. Collect species list and documented local webs with Carnivoran predator and prey richness greater than 2 and sufficient environmental covariates ====================
comm.web.carni.doc.sp <- data.frame()
comm.web.carni.doc.sp.SR3 <- data.frame()
for (i in 1:length(comm.carni.web.doc.link)){
  print(i)
  #extract links from communities 
  comm.web.carni.doc.sp <- rbind(comm.web.carni.doc.sp, data.frame(Community = names(comm.carni.web.doc.link[i]),
                                                       Species = unique(c(comm.carni.web.doc.link[[i]]$resource, comm.carni.web.doc.link[[i]]$consumer))))
  #extract links from communities with carnivora predator sp. > 2  & prey species > 2 to estimate dietary niche space
  if (length(unique(comm.carni.web.doc.link[[i]]$resource)) <3) next
  if (length(unique(comm.carni.web.doc.link[[i]][comm.carni.web.doc.link[[i]]$c.order == "Carnivora",]$consumer)) <3) next
  comm.web.carni.doc.sp.SR3 <- rbind(comm.web.carni.doc.sp.SR3, data.frame(Community = names(comm.carni.web.doc.link[i]),
                                                               Species = unique(c(comm.carni.web.doc.link[[i]]$resource, comm.carni.web.doc.link[[i]]$consumer))))
}
length(unique(comm.web.carni.doc.sp$Community)) #408 communities
length(unique(comm.web.carni.doc.sp$Species)) #350 consumer + resource species
length(unique(comm.web.carni.doc.sp.SR3$Community)) #392/512 communities meet criteria 
length(unique(comm.web.carni.doc.sp.SR3$Species)) #348 consumer + resource species

#Join region information 
Buff.10km.cov <- read.csv("Data/Rowan_Comm_10kmBuff_Env_cov.csv")
comm.web.carni.doc.sp <- left_join(comm.web.carni.doc.sp, Buff.10km.cov[, c("Community_unify", "Region")], join_by(Community == Community_unify))
comm.web.carni.doc.sp.SR3 <- left_join(comm.web.carni.doc.sp.SR3, Buff.10km.cov[, c("Community_unify", "Region")],join_by(Community == Community_unify))
#write.csv(comm.web.carni.doc.sp, "regional_408comm_doc_carni_links_splist_500g_add_afro_carni_links.csv", row.names = F)
#write.csv(comm.web.carni.doc.sp.SR3, "regional_408comm_doc_carni_links_splist_500g_add_afro_carni_links_SR3.csv", row.names = F)

comm.carni.web.doc.link.SR3 <- comm.carni.web.doc.link[which(names(comm.carni.web.doc.link) %in% unique(comm.web.carni.doc.sp.SR3$Community))]
length(comm.carni.web.doc.link.SR3) #392/512 communities meet criteria 
#save(comm.carni.web.doc.link, file = "community408_web_doc_carni_500g_with_afro_links.RData")
#save(comm.carni.web.doc.link.SR3, file = "community408_web_doc_carni_500g_with_afro_links_SR3.RData")

#III. Augment potential ecologically similar prey in local webs  ===============
##Choose food webs with sufficient environmental covariates and EVI >0 as finalized food webs for following analyses
Buff.10km.cov.sufficient.value <- Buff.10km.cov[Buff.10km.cov$EVI_month.mean > 0, ]
Buff.10km.cov.sufficient.value.web <- Buff.10km.cov.sufficient.value[Buff.10km.cov.sufficient.value$Community_unify %in% names(comm.carni.web.doc.link.SR3),]
length(Buff.10km.cov$Community_unify); length(Buff.10km.cov.sufficient.value.web$Community_unify)

#389 food webs have sufficient environmental covariates, EVI >0, and carnivora predator sp. > 2
comm.carni.web.doc.link.SR3.analyzed <- comm.carni.web.doc.link.SR3[which(names(comm.carni.web.doc.link.SR3) %in% unique(Buff.10km.cov.sufficient.value.web$Community_unify))]
length(comm.carni.web.doc.link.SR3.analyzed)

#Call cosine similarity for identified potential prey based on cosine similarity >= 0.95 & body mass as 10% differences: 0.9-1.1
region.potential.prey <- read.csv("Data/Regional_carni_links_500g_cosine_similarity_in95_10p_all.csv")
region.potential.prey.inf95.10p <- region.potential.prey[region.potential.prey$similarity  >= 0.95 & 
                                                           region.potential.prey$body.ratio.prey.pairs >= 0.9 & region.potential.prey$body.ratio.prey.pairs <= 1.1,]
#total 432 trait-matched prey species based on >= 0.95 cosince sim +- 10% body mass 
nrow(region.potential.prey.inf95.10p); length(unique(region.potential.prey.inf95.10p$potential_resources)) 


comm.carni.web.SR3.inf95.10p <- c()
for (i in 1:length(comm.carni.web.doc.link.SR3.analyzed)) {
  #subset potential prey by community species list 
  comm.sp.list <- community.sp.list.408comm.500g[community.sp.list.408comm.500g$Community  == names(comm.carni.web.doc.link.SR3.analyzed[i]),]
  comm.potential.prey <- region.potential.prey.inf95.10p[region.potential.prey.inf95.10p$consumer %in% comm.sp.list$Species &
                                                     region.potential.prey.inf95.10p$potential_resources %in% comm.sp.list$Species &
                                                     region.potential.prey.inf95.10p$region %in% comm.sp.list$Region, ]
  #subset potential prey by community species list 
  comm.potential.prey$intx <- paste0(comm.potential.prey$potential_resources, ".",comm.potential.prey$consumer)
  comm.potential.prey <- comm.potential.prey[,c("potential_resources", "consumer", "intx")]
  colnames(comm.potential.prey) <- c("resource", "consumer", "intx")
  
  #collect documented carnivora links in community
  comm.carni.links.doc.SR3 <- comm.carni.web.doc.link.SR3.analyzed[[i]]
  comm.carni.links.doc.SR3 <- comm.carni.links.doc.SR3[comm.carni.links.doc.SR3$c.order == "Carnivora",]
  comm.carni.links.doc.SR3$intx <- paste0(comm.carni.links.doc.SR3$resource, ".",comm.carni.links.doc.SR3$consumer)
  
  #join additional carnivora links to documented carni webs with predator SR >= 3
  comm.carni.links.all.SR3 <- rbind(comm.carni.links.doc.SR3[,c("resource", "consumer", "intx")],
                                comm.potential.prey)
  comm.carni.links.all.SR3 <- comm.carni.links.all.SR3[!duplicated(comm.carni.links.all.SR3$intx),]
  
  #add information of whether the links is documented 
  comm.carni.links.all.SR3$documeted <- 0
  comm.carni.links.all.SR3[comm.carni.links.all.SR3$intx %in% comm.carni.links.doc.SR3$intx, ]$documeted <- 1
  comm.carni.web.SR3.inf95.10p[[i]] <- comm.carni.links.all.SR3[,-c(3)]
}
names(comm.carni.web.SR3.inf95.10p) <- names(comm.carni.web.doc.link.SR3.analyzed)



# IV. Convert local web list into community collection as cheddar package format to calculate food web metrics ================
Buff.10km.cov.sufficient.value.web$title <- Buff.10km.cov.sufficient.value.web$Community_unify
length(comm.carni.web.SR3.inf95.10p) #389 food webs
length(Buff.10km.cov.sufficient.value.web$title)

# load webs with carnivora predator >= 3 with both documented  links
comm.carni.web.doc.link.SR3 <- comm.carni.web.doc.link.SR3[which(names(comm.carni.web.doc.link.SR3) %in% unique(Buff.10km.cov.sufficient.value.web$Community_unify))]
length(comm.carni.web.doc.link.SR3) #389 food webs

##1.Collect community-level node-properties and environmental variables in the cheddar format ---- 
Carni.web.doc <- list()           
for(i in 1:length(comm.carni.web.doc.link.SR3)){
  properties <- as.list(Buff.10km.cov.sufficient.value.web[Buff.10km.cov.sufficient.value.web$Community_unify == names(comm.carni.web.doc.link.SR3)[i],])
  links <- comm.carni.web.doc.link.SR3[[i]]
  species.node <- data.frame(Species = unique(c(links$consumer,links$resource)))
  nodes <- data.frame(node=species.node$Species, species.node)
  community.inferred <- Community(nodes=nodes, properties=properties, trophic.links=links)
  Carni.web.doc[[i]] <- community.inferred
}

### tropical predator prey web with documented links & 
###     inferred links with cosine similarity >= 0.95 & 10% of body mass
Carni.web.inf95.10p <- list()           
for(i in 1:length(comm.carni.web.SR3.inf95.10p)){
  properties <- as.list(Buff.10km.cov.sufficient.value.web[Buff.10km.cov.sufficient.value.web$Community_unify == names(comm.carni.web.SR3.inf95.10p)[i],])
  links <- comm.carni.web.SR3.inf95.10p[[i]]
  species.node <- data.frame(Species = unique(c(links$consumer,links$resource)))
  nodes <- data.frame(node=species.node$Species, species.node)
  community.inferred <- Community(nodes=nodes, properties=properties, trophic.links=links)
  Carni.web.inf95.10p[[i]] <- community.inferred
}

##3. Create community collection for cheddar function -----
Tropic.carni.webs.doc <- CommunityCollection(Carni.web.doc)
Tropic.carni.webs.inf95.10p <- CommunityCollection(Carni.web.inf95.10p)


save(#community collection of mammal food webs with documented trophic links
    Tropic.carni.webs.doc,
     #community collection of mammal food webs with documented and potential trophic links
     Tropic.carni.webs.inf95.10p, 
     file = "Data/Community408_carni_web_doc_inf95_500g_SR3_commcollect.Rdata")

