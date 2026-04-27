#  Spatial regression analyses
#
#  Input 
#   1. Rowan_Comm_10kmBuff_Env_cov.csv
#   2. PHYLACINE_LQ_sp_gain_loss_no_ocean_no_homo_410community.RData
#   3. Foodweb_lineage_history.csv
#   4. Foodweb_metrics.csv
#
#  Supplementary code to assist perform model average with spatial regression model
#   1. MuMIn_sarlm_model.avg.R
#   2. Sarlm_MuMIn_addon.R
#
#
#   Supplementary Figures S11-13, S15, and S17
#
#   Vectorized boundary Olson et al. 2001 BioSciences   https://doi.org/10.1641/0006-3568(2001)051[0933:TEOTWA]2.0.CO;2
#     polygon obtained from Data basin https://databasin.org/datasets/68635d7c77f1475f9b6c1d1dbe0a4c4c/
#     manually isolate Madagascar (Malagacy) from Afrotropics for visualization
#      1. wwf_terr_realm_Malagasy.shp (orginal file name: wwf_terr_realm.shp)
#
#

library(ggplot2); library(ggpubr); library(dotwhisker) #plotting library
library(dplyr)
library(psych) #correlation plot
library(MuMIn) #All possible combinations
library(sf)    #spatial data operation 
library(spatialreg) #spatial regression models
library(spdep)

setwd("./Data/")

#I. Combine all metrics and covariates ===========
##environmental covaraites 
buff.10km.env.cov <- read.csv("Rowan_Comm_10kmBuff_Env_cov.csv")
buff.10km.env.cov <- buff.10km.env.cov[buff.10km.env.cov$EVI_month.mean >0,]

##species loss 
load("PHYLACINE_LQ_sp_gain_loss_no_ocean_no_homo_410community.RData")

##present natural community and food web lineage history
evo.history <- read.csv("Foodweb_lineage_history.csv")

## food web metrics 
food.web.metrics <- read.csv( "Foodweb_metrics.csv")

##combine all information 
food.web.metircs.env.evo.hist <- left_join(buff.10km.env.cov[buff.10km.env.cov$Community_unify %in% food.web.metrics$Community,], LQ_SR_anomaly, join_by(Community_unify == Community))
food.web.metircs.env.evo.hist <- left_join(food.web.metircs.env.evo.hist, evo.history, join_by(Community_unify == Community))
food.web.metircs.env.evo.hist <- left_join(food.web.metircs.env.evo.hist, food.web.metrics, join_by(Community_unify == Community, Region == Region))

# II. Spatial regression model ==========================================
## 1. Continental differences in food web metrics  ---------------
nrow(food.web.metircs.env.evo.hist);table(food.web.metircs.env.evo.hist$Region) #389 food webs

food.web.metrics.analyzed.pt <- st_as_sf(food.web.metircs.env.evo.hist, coords = c("Long", "Lat"), crs = "EPSG:4326")
food.web.metrics.analyzed.pt$Region <- factor(food.web.metrics.analyzed.pt$Region, levels = c("Neotropic", "Afrotropic", "Indomalayan"))
### Generate spatial weight of defined neighbor for spatial statistics 
site.nb <- knn2nb(knearneigh(food.web.metrics.analyzed.pt, k=1, longlat = TRUE), row.names = food.web.metrics.analyzed.pt$Community)
### Generate spatial weight based on neighbor relationship 
site.nb.w <- nb2listw(site.nb, style="W", zero.policy=T)

###a. Food web species richness --------
SR.lm.raw <- spatialreg::errorsarlm(N.Sp.inf95.10p ~ Region, data = food.web.metrics.analyzed.pt, listw=site.nb.w)
LQ.SR.lm.raw <- spatialreg::errorsarlm(SR_LQ_anomaly ~ Region, data = food.web.metrics.analyzed.pt, listw=site.nb.w)

SR.lm <- spatialreg::errorsarlm(N.Sp.inf95.10p ~ Region-1, data = food.web.metrics.analyzed.pt, listw=site.nb.w) #to get marginal means for each group
LQ.SR.lm <- spatialreg::errorsarlm(SR_LQ_anomaly ~ Region-1, data = food.web.metrics.analyzed.pt, listw=site.nb.w)

SR.confidence <- rbind(data.frame(Est.mean = SR.lm$coefficients, lci = confint(SR.lm)[-1,1], uci = confint(SR.lm)[-1,2] , Status = "Present", Region = names(SR.lm$coefficients)),
                       data.frame(Est.mean = LQ.SR.lm$coefficients, lci = confint(LQ.SR.lm)[-1,1], uci = confint(LQ.SR.lm)[-1,2] , Status = "LQ anomaly", Region = names(SR.lm$coefficients)))
SR.confidence$Region <- factor(SR.confidence$Region, levels = c("RegionNeotropic", "RegionAfrotropic", "RegionIndomalayan"), labels = c("Neotropic", "Afrotropic", "Indomalaya"))
SR.confidence$Status <- factor(SR.confidence$Status, levels = c("Present", "LQ anomaly"), labels = c("Present", "Loss"))

###Visualize by confidence intervals 
status.col <- c("Present" = "forestgreen", "Loss" = "goldenrod3")
SR.est.v <- ggplot(SR.confidence, aes(x = Region, y =  Est.mean, color= Status)) +
  geom_point(size = 3, position =position_dodge(width=0.3)) + ylim(0,48) +
  geom_linerange(aes(x = Region, ymin = lci, ymax = uci), lwd = 1, position =position_dodge(width=0.3)) +
  scale_color_manual(values = status.col, labels = c("Present", "Lost")) + 
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.spacing=unit(1,"lines")) +
  theme(axis.text = element_text(size = 12), axis.title  = element_text(size = 12),
        legend.title = element_blank(), legend.text = element_text(size = 10), legend.position = c(0.82,0.82),
        legend.background = element_rect(fill = "NA", color = "grey50"), 
        strip.text.x = element_text(size = 14, face="bold"), strip.text.y = element_text(color = "white"), strip.background = element_rect(fill =NA, color = NA)) +
  labs(x="", y= "Estimated number of\nspecies", color = "")

###b. Predator/Basal prey ratio--------
food.web.metrics.analyzed.pt$cr.ratio <- food.web.metrics.analyzed.pt$N.Sp.c.inf95.10p/food.web.metrics.analyzed.pt$N.Sp.basal.inf95.10p

cr.ratio.lm.raw <- spatialreg::errorsarlm(cr.ratio ~ Region, data = food.web.metrics.analyzed.pt, listw=site.nb.w)
cr.ratio.lm <- spatialreg::errorsarlm(cr.ratio ~ Region-1, data = food.web.metrics.analyzed.pt, listw=site.nb.w)

cr.ratio.confidence <- data.frame(Est.mean = cr.ratio.lm$coefficients, lci = confint(cr.ratio.lm)[-1,1], uci = confint(cr.ratio.lm)[-1,2], Region = names(SR.lm$coefficients))
cr.ratio.confidence$Region <- factor(cr.ratio.confidence$Region, levels = c("RegionNeotropic", "RegionAfrotropic", "RegionIndomalayan"), labels = c("Neotropic", "Afrotropic", "Indomalaya"))

cr.ratio.est.v <- ggplot(cr.ratio.confidence, aes(x = Region, y =  Est.mean)) +
  geom_point(size = 3, position =position_dodge(width=0.3)) + ylim(0.25, 0.55) +
  geom_linerange(aes(x = Region, ymin = lci, ymax = uci), lwd = 1, position =position_dodge(width=0.3)) +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.spacing=unit(1,"lines")) +
  theme(axis.text = element_text(size = 12), axis.title  = element_text(size = 12),
        legend.title = element_blank(), legend.text = element_text(size = 10), legend.position.inside = c(0.79,0.82),
        legend.background = element_rect(fill = "NA", color = "grey50"), 
        strip.text.x = element_text(size = 14, face="bold"), strip.text.y = element_text(color = "white"), strip.background = element_rect(fill =NA, color = NA)) +
  labs(x="", y= "Estimated preadtor-basal\nprey richness ratio", color = "")


###c. predator and prey mass -----------------
c.bm.lm.raw <- spatialreg::errorsarlm(mass.log.c ~ Region, data = food.web.metrics.analyzed.pt, listw=site.nb.w)
r.bm.lm.raw <- spatialreg::errorsarlm(mass.log.r ~ Region, data = food.web.metrics.analyzed.pt, listw=site.nb.w)

c.bm.lm <- spatialreg::errorsarlm(mass.log.c ~ Region-1, data = food.web.metrics.analyzed.pt, listw=site.nb.w)
r.bm.lm <- spatialreg::errorsarlm(mass.log.r ~ Region-1, data = food.web.metrics.analyzed.pt, listw=site.nb.w)

BM.confidence <- rbind(data.frame(Est.mean = c.bm.lm$coefficients, lci = confint(c.bm.lm)[-1,1], uci = confint(c.bm.lm)[-1,2], Status = "Predator", Region = names(SR.lm$coefficients)),
                       data.frame(Est.mean = r.bm.lm$coefficients, lci = confint(r.bm.lm)[-1,1], uci = confint(r.bm.lm)[-1,2], Status = "Prey", Region = names(SR.lm$coefficients)))
BM.confidence$Region <- factor(BM.confidence$Region, levels = c("RegionNeotropic", "RegionAfrotropic", "RegionIndomalayan"), labels = c("Neotropic", "Afrotropic", "Indomalaya"))
BM.confidence$Status <- factor(BM.confidence$Status, levels = rev(c("Prey", "Predator")))
bm.col <- c("Predator" = "black", "Prey" = "grey60")


BM.predator.est.v <- ggplot(BM.confidence[BM.confidence$Status == "Predator",], aes(x = Region, y =  Est.mean)) +
  geom_point(size = 3, position =position_dodge(width=0.3)) + ylim(3.9, 4.24) +
  geom_linerange(aes(x = Region, ymin = lci, ymax = uci), lwd = 1, position =position_dodge(width=0.3)) +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.spacing=unit(1,"lines")) +
  theme(axis.text = element_text(size = 14),axis.title  = element_text(size = 12),
        strip.text.x = element_text(size = 14, face="bold"), strip.text.y = element_text(color = "white"), strip.background = element_rect(fill =NA, color = NA)) +
  labs(x="", y= expression(atop("Estimated ", paste(log[10],"(predator body mass [g])")))) 

BM.prey.est.v <- ggplot(BM.confidence[BM.confidence$Status == "Prey",], aes(x = Region, y =  Est.mean)) +
  geom_point(size = 3, position =position_dodge(width=0.3)) + ylim(3.65, 4.3) +
  geom_linerange(aes(x = Region, ymin = lci, ymax = uci), lwd = 1, position =position_dodge(width=0.3)) +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.spacing=unit(1,"lines")) +
  theme(axis.text = element_text(size = 14),axis.title  = element_text(size = 12),
        strip.text.x = element_text(size = 14, face="bold"), strip.text.y = element_text(color = "white"), strip.background = element_rect(fill =NA, color = NA)) +
  labs(x="", y= expression(atop("Estimated ", paste(log[10],"(prey body mass [g])")))) 


###d. normalized predator generality and prey vulnerability ---------
generality.predator.normalized.lm.raw <- spatialreg::errorsarlm(generality.sp.normalized ~ Region , data = food.web.metrics.analyzed.pt, listw=site.nb.w)
vulnerability.prey.normalized.lm.raw <- spatialreg::errorsarlm(vulnerability.sp.normalized ~ Region , data = food.web.metrics.analyzed.pt, listw=site.nb.w)

generality.predator.normalized.lm <- spatialreg::errorsarlm(generality.sp.normalized ~ Region -1, data = food.web.metrics.analyzed.pt, listw=site.nb.w)
vulnerability.prey.normalized.lm <- spatialreg::errorsarlm(vulnerability.sp.normalized ~ Region -1, data = food.web.metrics.analyzed.pt, listw=site.nb.w)


web.gene.vul.normalized.confidence <- rbind(data.frame(Est.mean = generality.predator.normalized.lm$coefficients, lci = confint(generality.predator.normalized.lm)[-1,1], uci = confint(generality.predator.normalized.lm)[-1,2], 
                                                       Region = names(generality.predator.normalized.lm$coefficients), 
                                                       Model = "Generality"),
                                            data.frame(Est.mean = vulnerability.prey.normalized.lm$coefficients, lci = confint(vulnerability.prey.normalized.lm)[-1,1], uci = confint(vulnerability.prey.normalized.lm)[-1,2], 
                                                       Region = names(vulnerability.prey.normalized.lm$coefficients), 
                                                       Model = "Vulnerability"))

web.gene.vul.normalized.confidence$Region <- factor(web.gene.vul.normalized.confidence$Region, levels = c("RegionNeotropic", "RegionAfrotropic", "RegionIndomalayan"), labels = c("Neotropic", "Afrotropic", "Indomalaya"))

web.vul.est.v <- 
  ggplot(web.gene.vul.normalized.confidence[web.gene.vul.normalized.confidence$Model == "Vulnerability",], aes(x = Region, y =  Est.mean)) +
  scale_color_manual(values = c("grey20", "grey60")) + 
  ylim(c(1,1.26)) +
  geom_point(size = 3, position =position_dodge(width=0.3)) + 
  geom_linerange(aes(x = Region, ymin = lci, ymax = uci), lwd = 1, position =position_dodge(width=0.3)) +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.spacing=unit(1,"lines")) +
  theme(axis.text = element_text(size = 13), axis.title  = element_text(size = 12),
        strip.text.x = element_text(size = 14, face="bold"), strip.text.y = element_text(color = "white"), strip.background = element_rect(fill =NA, color = NA)) +
  labs(x = "", y = "Estimated mean \nnormalized prey vulnerability")

web.gene.est.v <- 
  ggplot(web.gene.vul.normalized.confidence[web.gene.vul.normalized.confidence$Model == "Generality",], aes(x = Region, y =  Est.mean)) +
  geom_point(size = 3, position =position_dodge(width=0.3)) + 
  ylim(c(3.0,4.9)) +
  geom_linerange(aes(x = Region, ymin = lci, ymax = uci), lwd = 1, position =position_dodge(width=0.3)) +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.spacing=unit(1,"lines")) +
  theme(axis.text = element_text(size = 13), axis.title  = element_text(size = 12),
        strip.text.x = element_text(size = 14, face="bold"), strip.text.y = element_text(color = "white"), strip.background = element_rect(fill =NA, color = NA)) +
  labs(x = "", y = "Estimated mean \nnormalized predator generality")

#Figure 1 
ggarrange(NULL, SR.est.v+ rremove("x.text"), NULL, cr.ratio.est.v + rremove("x.text"), NULL,
                               NULL, BM.predator.est.v+ rremove("x.text"), NULL, BM.prey.est.v+ rremove("x.text"), NULL,  
                               NULL, web.vul.est.v, NULL, web.gene.est.v, NULL,  
                               nrow =3, ncol = 5, 
                               widths = c(0.02, 0.95, 0.04, 0.95, 0.1, 
                                          0.02, 0.95, 0.04, 0.95, 0.1,
                                          0.02, 0.95, 0.04, 0.95, 0.1), 
                               heights = c(0.6, 0.6, 0.6, 0.6, 0.6, 1.2, 1.2, 1.2, 1.2, 1.2), 
                               labels = c("","B", "", "C", "",  "", "D", "", "E","", "", "F", "", "G",""), hjust = 0.5,vjust = 0.5, align = "v")


###supp. normalized predator generality and prey vulnerability between binary and probabilistic links (Figure S7) ---------
generality.predator.normalized.body.size.model.lm <- spatialreg::errorsarlm(generality.weighted.body.size.model.p.normalized ~ Region -1, data = food.web.metrics.analyzed.pt, listw=site.nb.w)
vulnerability.prey.normalized.body.size.model.lm <- spatialreg::errorsarlm(vulnerability.weighted.body.size.model.p.normalized ~ Region -1, data = food.web.metrics.analyzed.pt, listw=site.nb.w)

generality.predator.normalized.allometric.model.lm <- spatialreg::errorsarlm(generality.weighted.allometric.model.p.normalized ~ Region -1, data = food.web.metrics.analyzed.pt, listw=site.nb.w)
vulnerability.prey.normalized.allometric.model.lm <- spatialreg::errorsarlm(vulnerability.weighted.allometric.model.p.normalized ~ Region -1, data = food.web.metrics.analyzed.pt, listw=site.nb.w)


web.gene.vul.confidence.sensitivity <- rbind(data.frame(Est.mean = generality.predator.normalized.lm$coefficients, lci = confint(generality.predator.normalized.lm)[-1,1], uci = confint(generality.predator.normalized.lm)[-1,2], 
                                                        Region = names(generality.predator.normalized.lm$coefficients), 
                                                        Model = "Generality",
                                                        TL.weighted = "Binary link"),
                                             data.frame(Est.mean = vulnerability.prey.normalized.lm$coefficients, lci = confint(vulnerability.prey.normalized.lm)[-1,1], uci = confint(vulnerability.prey.normalized.lm)[-1,2], 
                                                        Region = names(vulnerability.prey.normalized.lm$coefficients), 
                                                        Model = "Vulnerability",
                                                        TL.weighted = "Binary link"),
                                             ### by weighted interaction with body size model    
                                             data.frame(Est.mean = generality.predator.normalized.body.size.model.lm$coefficients, lci = confint(generality.predator.normalized.body.size.model.lm)[-1,1], uci = confint(generality.predator.normalized.body.size.model.lm)[-1,2], 
                                                        Region = names(generality.predator.normalized.body.size.model.lm$coefficients),
                                                        Model = "Generality",
                                                        TL.weighted = "Probabilistic link:\nbody-size model"),
                                             data.frame(Est.mean = vulnerability.prey.normalized.body.size.model.lm$coefficients, lci = confint(vulnerability.prey.normalized.body.size.model.lm)[-1,1], uci = confint(vulnerability.prey.normalized.body.size.model.lm)[-1,2], 
                                                        Region = names(vulnerability.prey.normalized.body.size.model.lm$coefficients),
                                                        Model = "Vulnerability",
                                                        TL.weighted = "Probabilistic link:\nbody-size model"),
                                             ### by weighted interaction by allometric model
                                             data.frame(Est.mean = generality.predator.normalized.allometric.model.lm$coefficients, lci = confint(generality.predator.normalized.allometric.model.lm)[-1,1], uci = confint(generality.predator.normalized.allometric.model.lm)[-1,2], 
                                                        Region = names(generality.predator.normalized.allometric.model.lm$coefficients),
                                                        Model = "Generality",
                                                        TL.weighted = "Probabilistic link:\nallometric niche model"),
                                             data.frame(Est.mean = vulnerability.prey.normalized.allometric.model.lm$coefficients, lci = confint(vulnerability.prey.normalized.allometric.model.lm)[-1,1], uci = confint(vulnerability.prey.normalized.allometric.model.lm)[-1,2], 
                                                        Region = names(vulnerability.prey.normalized.allometric.model.lm$coefficients),
                                                        Model = "Vulnerability",
                                                        TL.weighted = "Probabilistic link:\nallometric niche model"))

web.gene.vul.confidence.sensitivity <- web.gene.vul.confidence.sensitivity %>%  mutate_if(is.numeric, round, 2)
web.gene.vul.confidence.sensitivity$Region <- factor(web.gene.vul.confidence.sensitivity$Region, levels = c("RegionNeotropic", "RegionAfrotropic", "RegionIndomalayan"), labels = c("Neotropic", "Afrotropic", "Indomalaya"))

web.gene.sensitivity.est.v <- 
  ggplot(web.gene.vul.confidence.sensitivity[web.gene.vul.confidence.sensitivity$Model == "Generality",], aes(x = Region, y =  Est.mean, shape = TL.weighted)) +
  geom_point(size = 3, position =position_dodge(width=0.3)) + 
  ylim(c(3.0,4.9)) +
  geom_linerange(aes(x = Region, ymin = lci, ymax = uci), lwd = 1, position =position_dodge(width=0.3)) +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.spacing=unit(1,"lines")) +
  theme(axis.text = element_text(size = 14), axis.title  = element_text(size = 14),
        legend.title = element_text(size = 14),legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 14, face="bold"), strip.text.y = element_text(color = "white"), strip.background = element_rect(fill =NA, color = NA)) +
  labs(x = "", y = "Estimated mean normalized\npredator generality", shape = "")

web.vul.sensitivity.est.v <- 
  ggplot(web.gene.vul.confidence.sensitivity[web.gene.vul.confidence.sensitivity$Model == "Vulnerability",], aes(x = Region, y =  Est.mean, shape = TL.weighted)) +
  scale_color_manual(values = c("grey20", "grey60")) + 
  ylim(c(1,1.26)) +
  geom_point(size = 3, position =position_dodge(width=0.3)) + 
  geom_linerange(aes(x = Region, ymin = lci, ymax = uci), lwd = 1, position =position_dodge(width=0.3)) +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.spacing=unit(1,"lines")) +
  theme(axis.text = element_text(size = 14), axis.title  = element_text(size = 14),
        legend.title = element_text(size = 14),legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 14, face="bold"), strip.text.y = element_text(color = "white"), strip.background = element_rect(fill =NA, color = NA)) +
  labs(x = "", y = "Estimated mean normalized\nprey vulnerability", shape = "")

ggarrange(web.gene.sensitivity.est.v, web.vul.sensitivity.est.v, ncol = 2, nrow = 1, labels = c("A","B"), 
          font.label = list(family = "Arial", size = 16, color = "black", face = "bold"), common.legend = T, legend = "top")



###e. dietary niche breadth ---------
food.web.metrics.analyzed.pt$resource.hv.alpha.m.inf.log <- log10(food.web.metrics.analyzed.pt$resource.hv.alpha.m.inf)

hyperV.lm.raw <- spatialreg::errorsarlm(resource.hv.alpha.m.inf.log ~ Region, data = food.web.metrics.analyzed.pt, listw=site.nb.w)
hyperV.lm <- spatialreg::errorsarlm(resource.hv.alpha.m.inf.log ~ Region-1, data = food.web.metrics.analyzed.pt, listw=site.nb.w)

hyperV.confidence <- data.frame(Est.mean = hyperV.lm$coefficients, lci = confint(hyperV.lm)[-1,1], uci = confint(hyperV.lm)[-1,2], Region = names(hyperV.lm$coefficients))
hyperV.confidence$Region <- factor(hyperV.confidence$Region, levels = c("RegionNeotropic", "RegionAfrotropic", "RegionIndomalayan"), labels = c("Neotropic", "Afrotropic", "Indomalaya"))

hyperV.est.v <- 
  ggplot(hyperV.confidence, aes(x = Region, y =  Est.mean)) +
  geom_point(size = 3, position =position_dodge(width=0.3)) + ylim(4.42, 4.82) +
  geom_linerange(aes(x = Region, ymin = lci, ymax = uci), lwd = 1, position =position_dodge(width=0.3)) +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.spacing=unit(1,"lines")) +
  theme(axis.text = element_text(size = 12), axis.title  = element_text(size = 12),
        legend.title = element_blank(), legend.text = element_text(size = 10), legend.position.inside = c(0.79,0.82),
        legend.background = element_rect(fill = "NA", color = "grey50"), 
        strip.text.x = element_text(size = 14, face="bold"), strip.text.y = element_text(color = "white"), strip.background = element_rect(fill =NA, color = NA)) +
  labs(x="", y= "Estimated of \nlog10(mean dietary niche)\nbreadth per predator", color = "") 

###f. dietary niche overlap ---------
mean(food.web.metrics.analyzed.pt$resource.hv.beta.m.inf); sd(food.web.metrics.analyzed.pt$resource.hv.beta.m.inf)
hyperV.beta.lm.raw <- spatialreg::errorsarlm(resource.hv.beta.m.inf ~ Region, data = food.web.metrics.analyzed.pt, listw=site.nb.w)
hyperV.beta.lm <- spatialreg::errorsarlm(resource.hv.beta.m.inf ~ Region-1, data = food.web.metrics.analyzed.pt, listw=site.nb.w)

hyperV.beta.confidence <- data.frame(Est.mean = hyperV.beta.lm$coefficients, lci = confint(hyperV.beta.lm)[-1,1], uci = confint(hyperV.beta.lm)[-1,2], Region = names(hyperV.beta.lm$coefficients))
hyperV.beta.confidence$Region <- factor(hyperV.beta.confidence$Region, levels = c("RegionNeotropic", "RegionAfrotropic", "RegionIndomalayan"), labels = c("Neotropic", "Afrotropic", "Indomalaya"))

hyperV.beta.est.v <- ggplot(hyperV.beta.confidence, aes(x = Region, y =  Est.mean)) +
  geom_point(size = 3, position =position_dodge(width=0.3)) + ylim(c(0.12,0.28)) +
  geom_linerange(aes(x = Region, ymin = lci, ymax = uci), lwd = 1, position =position_dodge(width=0.3)) +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.spacing=unit(1,"lines")) +
  theme(axis.text = element_text(size = 12), axis.title  = element_text(size = 12),
        legend.title = element_blank(), legend.text = element_text(size = 10), legend.position.inside = c(0.79,0.82),
        legend.background = element_rect(fill = "NA", color = "grey50"), 
        strip.text.x = element_text(size = 14, face="bold"), strip.text.y = element_text(color = "white"), strip.background = element_rect(fill =NA, color = NA)) +
  labs(x="", y= "Estimated mean of\ndietary niche overlap\namong coexisting predators")

#Figure 2
ggarrange(NULL, hyperV.est.v+ rremove("x.text"), NULL, hyperV.beta.est.v , NULL, nrow =5, ncol = 1, 
          heights = c(0.05, 0.9, 0.02, 1, 0.05), labels = c(""," ", "", " ", ""), vjust = 0.5, hjust = 0.5, align = "v")



##2 Regional specific model ----------------------------------
## call to conduct model average
food.web.metircs.env.evo.hist$Elev_range.log <- log10(food.web.metircs.env.evo.hist$Elev_range)
food.web.metircs.env.evo.hist$resource.hv.alpha.m.inf.log <- log10(food.web.metircs.env.evo.hist$resource.hv.alpha.m.inf)

saving.dir <- "."
options('na.action' = "na.fail")

lm.variable <- c("SR_LQ_anomaly", "foodweb_meanDR","naturalpresent_comm_meanDR","Temp_deep_SD", "present_temp_sea", "EVI_month.mean", "EVI_seasonality", 
                 "Elev_mean", "Elev_range", "Elev_range.log","Shannon")

food.web.metrics <- c("resource.hv.alpha.m.inf.log", "resource.hv.beta.m.inf")

source("MuMIn_sarlm_model.avg.R")
#manually override functions listed in this Rscript: Sarlm_MuMIn_addon.R

summ.dt <- c()
summ.sw.dt <- c()

for (i in 1:length(unique(food.web.metircs.env.evo.hist$Region))) {
  print(i)
  Region <- sort(unique(food.web.metircs.env.evo.hist$Region))[i]
  
  Region.food.web.model <- food.web.metircs.env.evo.hist[food.web.metircs.env.evo.hist$Region == Region, ]
  
  Region.comm.pt <- st_as_sf(data.frame(Community =  Region.food.web.model$Community,
                                        long = Region.food.web.model$Long, lat = Region.food.web.model$Lat),
                             coords = c("long", "lat"), crs = st_crs(4326))
  
  ### Generate spatial weight of defined neighbor for spatial statistics 
  site.nb <- knn2nb(knearneigh(Region.comm.pt, k=1, longlat = TRUE), row.names = Region.comm.pt$Community)
  ### Generate spatial weight based on neighbor relationship 
  site.nb.w <- nb2listw(site.nb, style="W", zero.policy=T)
  
  #Subset for linear regression model variables
  Region.food.web.model.scale <- as.data.frame(scale(Region.food.web.model[,lm.variable]))
  
  #Region.food.web.model.scale$N.Sp <- Region.food.web.model$N.Sp
  Region.food.web.model.scale$Region <- Region.food.web.model$Region
  Region.food.web.model.scale$Community <- Region.food.web.model$Community
  Region.food.web.model.scale <- cbind(Region.food.web.model.scale, Region.food.web.model[,food.web.metrics])
  Region.food.web.model.scale[which(is.na(Region.food.web.model.scale) == T),]
  
  ###Hyper alpha community mean 
  RHV.err.full.inf <- spatialreg::errorsarlm(resource.hv.alpha.m.inf.log ~SR_LQ_anomaly+ naturalpresent_comm_meanDR + foodweb_meanDR + Temp_deep_SD + 
                                               present_temp_sea +  EVI_seasonality  + Elev_range.log + Shannon,
                                             data = Region.food.web.model.scale, listw=site.nb.w)
  RHV.err.full.inf.all <- MuMIn::dredge(RHV.err.full.inf, rank = "AIC",beta = T, evaluate = T)
  RHV.err.full.inf.95p <- model.avg.defult(get.models(RHV.err.full.inf.all, cumsum(weight) <= .95), beta=TRUE, fit=TRUE)
  
  ###Hyper beta community mean 
  RHV.beta.err.full.inf <- spatialreg::errorsarlm(resource.hv.beta.m.inf ~SR_LQ_anomaly+  naturalpresent_comm_meanDR + foodweb_meanDR + Temp_deep_SD + 
                                                    present_temp_sea + EVI_seasonality  + Elev_range.log + Shannon,
                                                  data = Region.food.web.model.scale, listw=site.nb.w)
  RHV.beta.err.full.inf.all <- MuMIn::dredge(RHV.beta.err.full.inf, rank = "AIC",beta = T, evaluate = T)
  RHV.beta.err.full.inf.95p <- model.avg.defult(get.models(RHV.beta.err.full.inf.all, cumsum(weight) <= .95), beta=TRUE, fit=TRUE)

  summ.dt <- rbind(summ.dt,
                   data.frame(Region = Region,
                              model = rep("RHV", nrow(summary(RHV.err.full.inf.95p)$coefmat.subset[-1,])),
                              var = row.names(summary(RHV.err.full.inf.95p)$coefmat.subset[-1,]),
                              summary(RHV.err.full.inf.95p)$coefmat.subset[-1,], 
                              lci = confint(RHV.err.full.inf.95p)[-1,1],
                              uci = confint(RHV.err.full.inf.95p)[-1,2]),
                   data.frame(Region = Region,
                              model = rep("RHV_overlap", nrow(summary(RHV.beta.err.full.inf.95p)$coefmat.subset[-1,])),
                              var = row.names(summary(RHV.beta.err.full.inf.95p)$coefmat.subset[-1,]),
                              summary(RHV.beta.err.full.inf.95p)$coefmat.subset[-1,], 
                              lci = confint(RHV.beta.err.full.inf.95p)[-1,1],
                              uci = confint(RHV.beta.err.full.inf.95p)[-1,2]))
  
  summ.sw.dt <- rbind(summ.sw.dt,
                        data.frame(Region = Region,
                                   Model = "RHV",
                                   var = names(RHV.err.full.inf.95p$sw),
                                   sw = round(RHV.err.full.inf.95p$sw,2)),
                        data.frame(Region = Region,
                                   Model = "RHV_overlap",
                                   var = names(RHV.beta.err.full.inf.95p$sw),
                                   sw = round(RHV.beta.err.full.inf.95p$sw,2)))
}

sarlm.95p <- summ.dt
sarlm.95p <- sarlm.95p[sarlm.95p$model %in% c("RHV", "RHV_overlap"), ]
sarlm.95p$Region <- factor(sarlm.95p$Region, levels = c("Neotropic", "Afrotropic","Indomalayan"), labels = c("Neotropic", "Afrotropic","Indomalaya"))

sarlm.95p <- sarlm.95p %>%  mutate_if(is.numeric, round, 4)
sarlm.95p$var <- factor(sarlm.95p$var, levels =c("naturalpresent_comm_meanDR", "Temp_deep_SD", "Elev_range.log", "SR_LQ_anomaly", "foodweb_meanDR","present_temp_sea","EVI_seasonality",  "Shannon"))
sarlm.95p$model <- factor(sarlm.95p$model, levels = c("RHV", "RHV_overlap"))
sarlm.95p.order <- sarlm.95p %>% arrange(factor(var)) %>% arrange(factor(Region)) %>% arrange(factor(model)) 

sarlm.95p$sig <- rep("non-significant", nrow(sarlm.95p))
sarlm.95p[sarlm.95p$Pr...z..<0.05 & sarlm.95p$Estimate > 0,]$sig <- "generalize"
sarlm.95p[sarlm.95p$Pr...z..<0.05 & sarlm.95p$Estimate < 0,]$sig <- "specialize"

foodweb.labs <- c("A. Mean dietary niche breadth", "B. Mean dietary niche overlap")
names(foodweb.labs) <- c("RHV", "RHV_overlap")

Region.color <- c("Neotropic" = "#adc9bc","Afrotropic" = "#ffb453" , "Indomalaya" = "#e6aaf0")
eff.col <- c("non-significant" = "grey80", "specialize" = "red3", "generalize" = "dodgerblue3")

### Figure 3 -------
ggplot(sarlm.95p, aes(y = var, x =  Estimate, col = sig)) +
  scale_color_manual(values = eff.col) +   scale_fill_manual(values = Region.color) +
  facet_grid(Region~model, labeller = labeller(model = foodweb.labs)) +
  geom_rect(aes(fill = Region), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.01, show.legend = F) +
  geom_vline(xintercept = 0, lty = 2, col = "grey60") +
  geom_linerange(aes(y = var, xmin = lci, xmax = uci), lwd = 1, position =position_dodge(width=0.65)) +
  geom_point(size = 3) +  theme_bw() +
  xlim(-0.8, 0.8) +
  scale_y_discrete(labels = c("Present natural\ncommunity lineage history","Temperature variability\nsince 3.3 Mya",  "Elevation range", 
                              "Species losses\nsince late Quaternary",  "Food-web lineage structure", "Temperature seasonality", "Productivity seasonality", "Habitat diversity")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing=unit(1,"lines"),
        axis.line = element_line(colour = "black"), legend.position = "bottom") +
  theme(axis.text.x = element_text(size = 12, angle= 50,  hjust = 1, vjust = 1), 
        axis.text.y  = element_text(size = 12),
        strip.text.x = element_text(size = 14, face="bold"), strip.text.y = element_text(size = 0),
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12), legend.position = "none",
        strip.background = element_rect(fill =NA, color = NA)) +
  labs(y="", x= "Model averaged standardized coefficient", color = "") + coord_flip()

### Figure S17 Spatial regression model: AIC weights --------
all.model.sw <- summ.sw.dt
all.model.sw$var <- factor(all.model.sw$var,  levels =c("naturalpresent_comm_meanDR", "Temp_deep_SD", "Elev_range.log", "SR_LQ_anomaly", "foodweb_meanDR","present_temp_sea","EVI_seasonality",  "Shannon"),
                           c("Present natural\ncommunity lineage history","Temperature variability\nsince 3.3 Mya",  "Elevation range", 
                             "Species losses\nsince late Quaternary",  "Food-web lineage structure", "Temperature seasonality", "Productivity seasonality", "Habitat diversity"))

all.model.sw$Model <- factor(all.model.sw$Model, levels = c("RHV", "RHV_overlap"),
                             labels = c("Mean dietary niche breadth", "Mean dietary niche overlap"))
all.model.sw$Region <- factor(all.model.sw$Region, levels = c("Neotropic","Afrotropic", "Indomalayan"), labels = c("Neotropic", "Afrotropic","Indomalaya"))
Region.color <- data.frame(Region = c("Neotropic", "Afrotropic","Indomalaya"), Region.color = c("#adc9bc","#ffb453","#e6aaf0"))

all.model.sw <- left_join(all.model.sw, Region.color, join_by(Region == Region))
all.model.sw$Region <- factor(all.model.sw$Region, levels =  c("Neotropic", "Afrotropic","Indomalaya"))
#write.csv(all.model.sw, "summary.model.variate.AICweight.csv", row.names = F)

Region.color <- c("Neotropic" = "#adc9bc","Afrotropic" = "#ffb453" , "Indomalaya" = "#e6aaf0")
model.col <- c("Mean dietary niche breadth" = "grey40","Mean dietary niche overlap" = "grey60")

ggplot(data = all.model.sw, aes(x= var, y = sw, fill = Model)) +
  ylim(c(0,1.14)) +
  facet_wrap(~Region, ncol = 1)+
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            fill = all.model.sw$Region.color, alpha = 0.01, show.legend = F) +
  geom_hline(yintercept = 0.5, lty = 2, color = "grey40") +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=sw), position=position_dodge(width=0.9), vjust=-0.25, color = "grey40") + 
  scale_fill_manual(values = model.col) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.grid = element_line("white"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.ticks  = element_line(size =1)) +
  theme(axis.text.x  = element_text(size = 16, angle= 60,  hjust = 1),
        axis.text.y  = element_text(size = 16),
        axis.title.y  = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18, face = "bold"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(face="bold"),
        panel.spacing.y = unit(1, "lines"))+
  theme(legend.position = "top") +
  labs(x = "", y = "Summed AIC model weight", fill = "") 



#Figure S9. Spatial maps of predatory dietary niche breadth and overlap ------------
### realm map obtained from Data basin https://databasin.org/datasets/68635d7c77f1475f9b6c1d1dbe0a4c4c/
### original file wwf_terr_realm.shp
Realm <- st_read("./wwf_terr_realm_Malagasy.shp")
unique(Realm$REALM)
Realm <- Realm[Realm$REALM %in% c("NT","NA", "AT", "MG", "PA", "IM", "OC"), ]
Realm.study <- Realm[Realm$REALM %in% c("NT", "AT", "IM"), ]

Region.col <- c("NT" = "#adc9bc","AT" = "#ffb453" , "IM" = "#e6aaf0")
region.bound <- data.frame(x1 = c(-130,-130), x2 = c(140, 140), y1 = c(-23.43,-35), y2 = c(23.43, 35),
                           label = c("Tropics","Subtropics"))

study.site <- ggplot() + 
  geom_rect(data = region.bound, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), fill = c("#B9DDF1", "#B9DDF1"), alpha = 0.4) +
  #geom_text(data = region.bound, aes(x = -98, y = -30, label = label)) +
  geom_text(aes(x = -105, y = -19, label = "Tropics"), color = "dodgerblue4", alpha = 0.8) +
  geom_text(aes(x = -101, y = -30.5, label = "Subtropics"), color = "dodgerblue", alpha = 0.6) +
  geom_sf(data = Realm, aes(fill = REALM), fill =  "grey90", alpha = 1, color = "grey80", show.legend = F) + 
  geom_sf(data = Realm.study, aes(fill = REALM), alpha = 1, color = "NA", show.legend = F) + 
  geom_sf(data = food.web.metrics.analyzed.pt, aes(color = N.Sp.inf95.10p), size = 1.2) +
  scale_fill_manual(values = Region.col) + scale_color_viridis_c() + theme_bw() + 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14, color = "grey20"),
        legend.text = element_text(size = 12), legend.title = element_text(size = 14, color = "grey20"),
        panel.grid = element_blank()) +
  coord_sf(xlim=c(-105, 120),ylim=c(-55,35), crs = "EPSG:4326") +
  labs(x = "Longitude", y = "Latitude",fill = "Region", color = "Species\nrichness")

RHV.map <- ggplot() + 
  geom_sf(data = Realm, aes(fill = REALM), fill =  "grey90", alpha = 1, color = "grey80", show.legend = F) + 
  geom_sf(data = Realm.study, aes(fill = REALM), alpha = 1, color = "NA", show.legend = F) + 
  geom_sf(data = food.web.metrics.analyzed.pt, aes(color = log10(resource.hv.alpha.m.inf)), size = 1.2) +
  scale_fill_manual(values = Region.col) + scale_color_viridis_c() + theme_bw() + 
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14, color = "grey20"),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, color = "grey20"),
        panel.grid = element_blank()) +
  coord_sf(xlim=c(-105, 120),ylim=c(-55,35), crs = "EPSG:4326") +
  labs(x = "Longitude", y = "Latitude",fill = "Region", color = "Mean dietary niche breadth\nper predator")

overlap.map <- ggplot() + 
  geom_sf(data = Realm, aes(fill = REALM), fill =  "grey90", alpha = 1, color = "grey80", show.legend = F) + 
  geom_sf(data = Realm.study, aes(fill = REALM), alpha = 1, color = "NA", show.legend = F) + 
  geom_sf(data = food.web.metrics.analyzed.pt, aes(color = resource.hv.beta.m.inf), size = 1.2) +
  scale_fill_manual(values = Region.col) + scale_color_viridis_c() + theme_bw() + 
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14, color = "grey20"),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, color = "grey20"),
        panel.grid = element_blank()) +
  coord_sf(xlim=c(-105, 120),ylim=c(-55,35), crs = "EPSG:4326") +
  labs(x = "Longitude", y = "Latitude",fill = "Region", color = "Mean dietary niche overlap\namong coexisting predator")

LQ_Rich_Anom.map <- ggplot() + 
  geom_sf(data = Realm, aes(fill = REALM), fill =  "grey90", alpha = 1, color = "grey80", show.legend = F) + 
  geom_sf(data = Realm.study, aes(fill = REALM), alpha = 1, color = "NA", show.legend = F) + 
  geom_sf(data = food.web.metrics.analyzed.pt, aes(color = SR_LQ_anomaly), size = 1.2) +
  scale_fill_manual(values = Region.col) + scale_color_viridis_c() + theme_bw() + 
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14, color = "grey20"),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, color = "grey20"),
        panel.grid = element_blank()) +
  coord_sf(xlim=c(-105, 120),ylim=c(-55,35), crs = "EPSG:4326") +
  labs(x = "Longitude", y = "Latitude",fill = "Region", color = "Number of speciess losses \nsince late Quaternary")

ggarrange(RHV.map, overlap.map, LQ_Rich_Anom.map, nrow = 3, ncol = 1, align = "hv", labels = c("A", "B", "C"))



#Figure S10. Boxplots of Regional variations in Env. ---------------
library(reshape2)
env.evo.plot <- melt(food.web.metircs.env.evo.hist[,c("Community", "Region", "naturalpresent_comm_meanDR", "Temp_deep_SD","Elev_range.log", 
                                                  "SR_LQ_anomaly",  "foodweb_meanDR", "present_temp_sea", "EVI_seasonality", "Shannon")], id.vars=c("Community", "Region"))
env.evo.plot$Region <- factor(env.evo.plot$Region, levels = c("Neotropic", "Afrotropic", "Indomalayan"), labels = c("Neotropic", "Afrotropic", "Indomalaya"))
Region.color <- c("Neotropic" = "#adc9bc","Afrotropic" = "#ffb453" , "Indomalaya" = "#e6aaf0")

col.name <- c( "naturalpresent_comm_meanDR" = "Present natural\ncommunity lineage history",              
               "Temp_deep_SD" = "Temperature variability\nsince 3.3 Mya","Elev_range.log" = "log10(elevation range)",  
               "SR_LQ_anomaly" = "Species losses\nsince late Quaternary", 
               "foodweb_meanDR" = "Food-web lineage structure",
               "present_temp_sea" =  "Temperature seasonality", "EVI_seasonality" = "Productivity seasonality", "Shannon" =  "Habitat diversity")
ggplot(data = env.evo.plot) +
  facet_wrap(~variable, scales = "free", labeller = as_labeller(col.name), ncol =2) +
  geom_boxplot(data = env.evo.plot, aes(x = Region, y = value, fill=Region)) + 
  scale_fill_manual(values = Region.color) + 
  theme(strip.text = element_text(size = 12, face = "bold"), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),panel.border = element_rect(linetype = "solid", fill = NA),
        strip.background = element_blank(),
        legend.position = c("none"), legend.title = element_blank(),legend.key = element_blank(), 
        legend.text = element_text(size = 14)) + scale_y_continuous(expand = expansion(mult = c(0.3))) +
  labs(x = "", y = "")


#Figure S11 Correlations among the environmental predictors --------
library(GGally)

predictor <- c("naturalpresent_comm_meanDR", "Temp_deep_SD", "Elev_range.log", "SR_LQ_anomaly", "foodweb_meanDR", "present_temp_sea","EVI_seasonality",  "Shannon")
pred.label <- c("Community lineage", "Temp.deep.SD", "Elev.range.log", "Hist. sp. losses", "Foodweb lineage", "Temp.seasonality", "Prod.seasonality", "Hab.diversity")

###(1) Neo
lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_smooth(method = method, color = "grey80", fill = "grey90", alpha = 0.8,...) +
    geom_point(pch=21, fill = "#adc9bc", color = "grey60", alpha = 0.8, stroke = 0.1)
}

ggpairs(food.web.metircs.env.evo.hist[food.web.metircs.env.evo.hist$Region == "Neotropic", predictor], 
        lower = list(continuous = wrap(lowerFn, method = "lm")),
        diag = list(continuous = wrap("barDiag", colour = "NA", fill = "#adc9bc")),
        upper = list(continuous = wrap("cor", size = 5)),
        columnLabels = pred.label) +
  theme_bw()+ theme(panel.grid = element_blank(), panel.background  = element_rect(fill = "white"))

###(2) Afro
lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_smooth(method = method, color = "grey80", fill = "grey90", alpha = 0.8,...) +
    geom_point(pch=21, fill = "#ffb453", color = "grey60", alpha = 0.8, stroke = 0.1) 
  p
}

ggpairs(food.web.metircs.env.evo.hist[food.web.metircs.env.evo.hist$Region == "Afrotropic", predictor], 
        lower = list(continuous = wrap(lowerFn, method = "lm")),
        diag = list(continuous = wrap("barDiag", colour = "NA", fill = "#ffb453")),
        upper = list(continuous = wrap("cor", size = 5)),
        columnLabels = pred.label) +
  theme_bw()+ theme(panel.grid = element_blank(), panel.background  = element_rect(fill = "white"))

###(3) Indo
lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_smooth(method = method, color = "grey80", fill = "grey90", alpha = 0.8,...) +
    geom_point(pch=21, fill = "#e6aaf0", color = "grey60", alpha = 0.8, stroke = 0.1) 
  p
}

ggpairs(food.web.metircs.env.evo.hist[food.web.metircs.env.evo.hist$Region == "Indomalayan", predictor], 
        lower = list(continuous = wrap(lowerFn, method = "lm")),
        diag = list(continuous = wrap("barDiag", colour = "NA", fill = "#e6aaf0")),
        upper = list(continuous = wrap("cor", size = 5)),
        columnLabels = pred.label) +
  theme_bw()+ theme(panel.grid = element_blank(), panel.background  = element_rect(fill = "white"))


#Figure S15 Correlation among response variables of predator dietary niche breadth and overlap ----------
lowerFn <- function(data, mapping, method = "lm",  ...) {
  p <- ggplot(data = data, mapping = mapping, aes(fill = Region, color = Region)) +
    geom_smooth(method = method, alpha = 0.2) +
    geom_point(pch=21, color = "grey60", alpha = 0.8, stroke = 0.1) 
  p
}

food.web.metircs.env.evo.hist$Region <- factor(food.web.metircs.env.evo.hist$Region, levels = c("Neotropic", "Afrotropic", "Indomalayan"))
Region.color <- c("Neotropic" = "#adc9bc","Afrotropic" = "#ffb453" , "Indomalayan" = "#e6aaf0")

ggpairs(food.web.metircs.env.evo.hist[,c("Region",  "N.Sp.c.inf95.10p", "resource.hv.alpha.m.inf.log", "resource.hv.beta.m.inf") ], 
        columns = 2:4, columnLabels = c("Predator richness", "log10-Mean dietary niche breadth", "Mean dietary niche overlap"),
        lower = list(continuous = wrap(lowerFn, method = "lm")),
        diag = list(continuous = wrap("barDiag", colour = "grey80")),
        upper = list(continuous = wrap("cor", size = 5)),
        ggplot2::aes(colour=Region, fill = Region, alpha = 0.6)) + scale_fill_manual(values = Region.color) + scale_color_manual(values = Region.color) + 
  theme_bw()+ theme(panel.grid = element_blank(), panel.background  = element_rect(fill = "white"),
                    strip.text = element_text(size = 10, face = "bold"))

