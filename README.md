-----------------------------------------------

# [Historical legacies shape continental variation in contemporary mammal food webs](https://doi.org/10.1073/pnas.2519938123)

### [Chia Hsieh](https://chiahsieheeb.wixsite.com/chiahsieh), [Evan C. Fricke](https://www.evanfricke.com), [Wei-Hao Lee](https://sites.google.com/view/wei-haolee/home), [Daniel Gorczynski](https://danielgorczynski.weebly.com), and [Lydia Beaudrot](https://lydiabeaudrot.weebly.com)


### Latest version of Data and Code DOI: [https://doi.org/10.5281/zenodo.17821933](https://doi.org/10.5281/zenodo.17821933)

### Please contact the corresponding author, Chia Hsieh (hsiehch@msu.edu or chiahsieheeb@gmail.com), for questions about the code or data.


-----------------------------------------------

## Abstract

Contemporary ecological networks reflect the influence of evolutionary and ecological processes on species interactions. Additionally, late Quaternary megafaunal extinctions and subsequent range contractions have interrupted mammalian predator-prey interactions to differing degrees among continents. However, the extent to which species losses have shaped geographic variation in the vertical and horizontal structure of contemporary mammalian food webs remains poorly understood. Furthermore, the relative influences of species loss, evolutionary history, and contemporary environmental drivers on the structure of mammalian food webs have not been tested. Here, we assembled mammal food webs for 389 sites throughout the Neotropics, Afrotropics, and Indomalaya. We first examined variation among continents in food web properties, including predator-prey richness ratios, predator generality, prey vulnerability, and predator dietary specialization. We found that Neotropical food webs bear imprints of severe species loss, featuring smaller predators, lower predator generality, higher prey vulnerability, and reduced prey species availability compared to Afrotropical food webs, in which most large mammals persisted. Furthermore, we reveal that multiple historical legacies, including community lineage history, paleoclimatic variability, mountain uplift, and species losses, as well as contemporary environmental variability, collectively predicted variation in predator dietary niche breadth within continents. Our findings offer new perspectives on how trophic interactions within food webs respond to species losses and how macroevolutionary and macroecological processes have shaped the biogeography of contemporary mammal food webs. 

---

## Repository Directory 

1. Code: Contains code files used to derive non-proprietary data for identifying potential prey and quantifying hypervolume of predator dietary niche space based on prey trait space, as well as regression model covariates. 
2. Data: Contains derived non-proprietary data supporting this project

---

## Code and Workflow

### I. Extract environmental covariates and estimate species loss for each community

  1. [GGE_10kmBuff_EVI_and_SRTMDEM.txt](Code/GGE_10kmBuff_EVI_and_SRTMDEM.txt): Google Earth Engine script to obtain EVI and elevation values
       a. Input: [Rowan_2020_Community_10kmBuff.zip](Data/Rowan_2020_Community_10kmBuff.zip)
       b. Output: [Rowan_Comm_10kmBuff_EVI_2000-2016.csv](Data/Rowan_Comm_10kmBuff_EVI_2000-2016.csv) and [Rowan_Comm_10kmBuff_elevation.csv](Data/Rowan_Comm_10kmBuff_elevation.csv)
  
  2. [1_1_Community_Env_10km.R](Code/1_1_Community_Env_10km.R): R script used to extract environmental covariates within 10-km radius buffer of each community location. 410 communities across Neotropics, Afrotropics, and Indomalaya. 
       a. Input: Published community checklist from [Rowan et al. 2020 PNAS](<https://www.pnas.org/doi/10.1073/pnas.1910489116#supplementary-materials>) and public environmental layers with repository link described in script
       a. Output [Rowan_Comm_10kmBuff_Env_cov.csv](Data/Rowan_Comm_10kmBuff_Env_cov.csv)
       
  3. [1_2_Community_LQspeciesloss_10km.R](Code/1_2_Community_LQspeciesloss_10km.R): R script used to estimate historical species loss per community and per region. 
       a. Input: Public species range map from [PHYLACINE 1.2 dataset](https://megapast2future.github.io/PHYLACINE_1.2/)
       b. Output [PHYLACINE_LQ_sp_gain_loss_no_ocean_no_homo_410community.RData](Data/PHYLACINE_LQ_sp_gain_loss_no_ocean_no_homo_410community.RData)
  

### II. Regional meta-web based on documented and potential predator-prey pairs
       
  1. [2_Regional_metaweb_potential_prey.R](Code/2_Regional_metaweb_potential_prey.R): R script used to measure food web metrics: predator-basal prey ratio, mean body mass of predator and prey, normalized generality and normalized vulnerability, and mean predator dietary niche breadth and overlap. 
       a. Input [Region_sp_intx_500g.RData](Data/Region_sp_intx_500g.RData) 
       b. Output [Regional_carni_links_500g_cosine_similarity_in95_10p_all.csv](Data/Regional_carni_links_500g_cosine_similarity_in95_10p_all.csv)

### III. Local food web and measur food web metrics
  1. [3_1_Community_trophic_links.R](Code/3_1_Community_trophic_links.R): R script to augment local food webs with identified potential prey based on cosine similarity approach. 
       a. Input [Region_sp_intx_500g.RData](Data/Region_sp_intx_500g.RData), [Rowan_Comm_10kmBuff_Env_cov.csv](Data/Rowan_Comm_10kmBuff_Env_cov.csv), and [Regional_carni_links_500g_cosine_similarity_in95_10p_all.csv](Data/Regional_carni_links_500g_cosine_similarity_in95_10p_all.csv).
       b. Output [Community408_carni_web_doc_inf95_500g_SR3_commcollect.Rdata](Data/Community408_carni_web_doc_inf95_500g_SR3_commcollect.Rdata)
       
  2. [3_2_Foodweb_metrics.R](Code/3_2_Foodweb_metrics.R): R script to measure food web metrics: predator-basal prey ratio, mean body mass of predator and prey, normalized generality and normalized vulnerability, and mean predator dietary niche breadth and overlap. R script includes code to make *Supplementary Figure S6* to demonstrate density distribution of proportion of the total predator species and basal prey species across the 389 mammal food webs by realm and *Supplementary Figure S8* to demonstrate Association between the number of prey and dietary niche breadth for each predator within regional meta-webs. Code for *Supplementary Figure S14*.
      a. Input [Community408_carni_web_doc_inf95_500g_SR3_commcollect.Rdata](Data/Community408_carni_web_doc_inf95_500g_SR3_commcollect.Rdata), [Region_sp_intx_500g.RData](Data/Region_sp_intx_500g.RData), and [Regional_carni_links_500g_cosine_similarity_in95_10p_all.csv](Data/Regional_carni_links_500g_cosine_similarity_in95_10p_all.csv). 
      b. Output for [3_3_Predator_dietary_niche_hypervolume.py](Code/3_3_Predator_dietary_niche_hypervolume.py) to estimate predator dietary niche breadth: [Predator_prey_pairs_for_dietary_niche_space.csv](Data/Predator_prey_pairs_for_dietary_niche_space.csv) and [Gower_11traits_PCoA_Axis13_inf95_10p_500g.csv](Data/Gower_11traits_PCoA_Axis13_inf95_10p_500g.csv). 
      c. Final output [Foodweb_metrics.csv](Data/Foodweb_metrics.csv). 
      
  3. [3_3_Predator_dietary_niche_hypervolume.py](Code/3_3_Predator_dietary_niche_hypervolume.py): Python script used to quantify dietary niche hypervolume for each predator and the overlapping hypervolume between predator pairs using prey trait space. 
      a. Input [Predator_prey_pairs_for_dietary_niche_space.csv](Data/Predator_prey_pairs_for_dietary_niche_space.csv) and [Gower_11traits_PCoA_Axis13_inf95_10p_500g.csv](Data/Gower_11traits_PCoA_Axis13_inf95_10p_500g.csv). 
      b. Output [Predator_dietary_niche_hypervolume.csv](Data/Predator_dietary_niche_hypervolume.csv) and [Predator_dietary_niche_overlap_hypervolume.csv](Data/Predator_dietary_niche_overlap_hypervolume.csv). Replace predator-specific prey to total or regional prey species to measure hypervolume as [Regional_predator_dietary_hypervolume.csv](Data/Regional_predator_dietary_hypervolume.csv).

### IV. Food web lineage history 
 1. [4_Foodweb_lineage_history.R](Code/4_Foodweb_lineage_history.R): R script to estimate natural present community lineage history and food web lineage structure using species-specific DR metric. R script includes code for *Supplementary Figure S4* for phylogeny and *Supplementary Figure S16* to demonstrate distribution of Relationship between species-level body mass and species-specific lineage diversification history.
      a. Input [Community408_carni_web_doc_inf95_500g_SR3_commcollect.Rdata](Data/Community408_carni_web_doc_inf95_500g_SR3_commcollect.Rdata), [PHYLACINE_LQ_sp_gain_loss_no_ocean_no_homo_410community.RData](Data/PHYLACINE_LQ_sp_gain_loss_no_ocean_no_homo_410community.RData), and [Rowan_Comm_10kmBuff_Env_cov.csv](Data/Rowan_Comm_10kmBuff_Env_cov.csv).
      a. Output [Foodweb_lineage_history.csv](Data/Foodweb_lineage_history.csv). 

### V. Spatial regression analyses
 1. [5_Spatial_regression.R](Code/5_Spatial_regression.R): R script to perform spatial regression models for results in *Figures 1-3*. R script includes code to demonstrate regional variations in model covariates in *Supplementary Figure S10*, correlation plots of model covariates as *Supplementary Figures S11-13*, correlation plots among response variable as *Supplementary Figures S15*, and summed AIC model weights as *Supplementary Figure S17*. Code for spatial maps of *Figure 1* and *Supplementary Figure S9*.
      a. Input [Rowan_Comm_10kmBuff_Env_cov.csv](Data/Rowan_Comm_10kmBuff_Env_cov.csv), [PHYLACINE_LQ_sp_gain_loss_no_ocean_no_homo_410community.RData](Data/PHYLACINE_LQ_sp_gain_loss_no_ocean_no_homo_410community.RData), [Foodweb_lineage_history.csv](Data/Foodweb_lineage_history.csv), and [Foodweb_metrics.csv](Data/Foodweb_metrics.csv). 
      
 2. [Sarlm_MuMIn_addon.R](Code/Sarlm_MuMIn_addon.R): R script to manually override model average functions for spatial regression model output
 3. [MuMIn_sarlm_model.avg.R](Code/MuMIn_sarlm_model.avg.R): R script to support model average with spatial regression model output
 
### VI. Supplemtary Tables and Figures
 1. [Supp_Regional_SpLoss.R](Code/Supp_Regional_SpLoss.R): R script to summarize regional species loss for *Supplementary Tables S3 and S4* and *Supplementary Figure S5*
     a. Input [Rowan_Comm_10kmBuff_Env_cov.csv](Data/Rowan_Comm_10kmBuff_Env_cov.csv), [Foodweb_metrics.csv](Data/Foodweb_metrics.csv), and [PHYLACINE_LQ_sp_gain_loss_no_ocean_no_homo_410community.RData](Data/PHYLACINE_LQ_sp_gain_loss_no_ocean_no_homo_410community.RData). 
     
 2. [Supp_Regional_metaweb_trait_Vis.R](Code/Supp_Regional_metaweb_trait_Vis.R): R script to visualize regional meta-web and trait values used to quantify prey hypervolume for *Supplementary Figures S1 to S3*.
     a. Input [Community408_carni_web_doc_inf95_500g_SR3_commcollect.Rdata](Data/Community408_carni_web_doc_inf95_500g_SR3_commcollect.Rdata)
     
 3. [Supp_Predator_dietary_niche_Vis.R](Code/Supp_Predator_dietary_niche_Vis.R): R script to visualize predator dietary niche overlap in *Figure 2* 
      a. Input [Predator_prey_pairs_for_dietary_niche_space.csv](Data/Predator_prey_pairs_for_dietary_niche_space.csv) and [Gower_11traits_PCoA_Axis13_inf95_10p_500g.csv](Data/Gower_11traits_PCoA_Axis13_inf95_10p_500g.csv). 
 
 
 
----

## Data 
### Shapefile of community buffer for Google Earth Engine
1. [Rowan_2020_Community_10kmBuff.zip](Data/Rowan_2020_Community_10kmBuff.zip)
    a. Polygon of community-specific 10-km-radius buffer of community checklist location to extract EVI and elevation values

### Model covariates 
1. [Rowan_Comm_10kmBuff_EVI_2000-2016.csv](Data/Rowan_Comm_10kmBuff_EVI_2000-2016.csv)
     a. Community: Community name to match community checklists and their coordinates at [Rowan et al. 2020 PNAS](<https://www.pnas.org/doi/10.1073/pnas.1910489116#supplementary-materials>)
     b. monthly EVI values by columns

2.  [Rowan_Comm_10kmBuff_elevation.csv](Data/Rowan_Comm_10kmBuff_elevation.csv)
     a. Community: Community name to match community checklists and their coordinates at [Rowan et al. 2020 PNAS](<https://www.pnas.org/doi/10.1073/pnas.1910489116#supplementary-materials>)
     b. mean, standard deviation, minimum, and maximum elevation values

3. [Rowan_Comm_10kmBuff_Env_cov.csv](Data/Rowan_Comm_10kmBuff_Env_cov.csv) 
    a. Community: Community name to match community checklists and their coordinates at [Rowan et al. 2020 PNAS](<https://www.pnas.org/doi/10.1073/pnas.1910489116#supplementary-materials>)
    b. Community_unify: unified community name with special character to match communities among data
    c. Long: longitude of community coordinates from [Rowan et al. 2020 PNAS](<https://www.pnas.org/doi/10.1073/pnas.1910489116#supplementary-materials>)
    d. Lat: latitude of community coordinates from [Rowan et al. 2020 PNAS](<https://www.pnas.org/doi/10.1073/pnas.1910489116#supplementary-materials>)
    e. Region: Biogeographic realms 
    f. Temp_deep_Mean: average annual mean temperature (BIO1) over twelve periods since the late Pliocene (~3.3 Mya)
    g. Temp_deep_SD: standard deviation of annual mean temperature (BIO1) over twelve periods since the late Pliocene (~3.3 Mya) 
    h. present_temp: annual mean temperature (BIO1) for 1979–2013  
    i. present_temp_sea: temperature seasonality (BIO4) for 1979–2013
    j. EVI_month.mean: mean monthly productivity during 2000 to 2016
    k. EVI_seasonality: productivity seasonality during 2000 to 2016
    l. Elev_mean: mean elevantion
    m. Elev_range: elevation range
    n. Shannon: habitat diversity

4. [PHYLACINE_LQ_sp_gain_loss_no_ocean_no_homo_410community.RData](Data/PHYLACINE_LQ_sp_gain_loss_no_ocean_no_homo_410community.RData)
    a. nature.present.comm.SpList: species list of present natural communities (N = 410) 
    b. current.comm.SpList: species list of current communities (N = 410)
    c. LQ_SR_anomaly: species richness difference (SR_LQ_anomaly) between present natural (N.sp.nature) and current (N.sp.current) communities
    d. LQ.loss.gain.comm.SpList: identified lost and gained (LossGain) species in current community relative to present natural community 

### Regional meta webs 
1. [Region_sp_intx_500g.RData](Data/Region_sp_intx_500g.RData): Input data in [2_2_Regional_metaweb_potential_prey.R](Code/2_2_Regional_metaweb_potential_prey.R)
     a. community.sp.list.408comm.500g: species list (>= 500g) of 408 communities containing trophic links 
     b. comm.web: community collection of documented predator-prey pairs of 408 communities, in format prepared for food web metric measurements with R package *cheddar*
     c. df_region.match: unique species list of 408 community checklists in Neotropic (neo), Afrotropic (afr), and Indomalaya (ind) at [Rowan et al. 2020 PNAS](<https://www.pnas.org/doi/10.1073/pnas.1910489116#supplementary-materials>) with 686 species with body mass >= 500g and sufficient ecological traits. This file is used to identify potential prey for each predator across communities that are not recorded.
     d. df_region_intx.match: unique documented pairs (intx) of predator (consumer) prey (resource) among three realms (Region). Data source of these predator-prey pairs in regional meta-webs: Global meta-web from [Fricke et al. 2022 Science](https://doi.org/10.1126/science.abk3510) (Global.web)
          
2. [Regional_carni_links_500g_cosine_similarity_in95_10p_all.csv](Data/Regional_carni_links_500g_cosine_similarity_in95_10p_all.csv): Output data from [2_1_Regional_metaweb_links.R](Code/2_1_Regional_metaweb_links.R) 
     a. Cosine similarity value (similarity) between a documented prey (resource; taxonomic order of prey: order_r) with potential prey (potential_resources; taxonomic order of potential prey: order_p_r) for a focal predator (consumer; taxonomic order of predator: order_p_r) based on regional meta-web (region). Potential prey is identified by criteria of cosine similarity >= 0.95 (in95) and 10% body mass differences (10p) relative to documented prey (0.9 <=body.ratio.prey.pairs <= 1.1), and marked as undocumented (documented = 0).    
 
### Local food webs
1. [Community408_carni_web_doc_inf95_500g_SR3_commcollect.Rdata](Data/Community408_carni_web_doc_inf95_500g_SR3_commcollect.Rdata)
     a. Tropic.carni.webs.doc: community collection of mammal food webs with documented trophic links with only predators in Order of Carnivora and food webs contain at least three predator species, in format prepared for food web metric measurements with R package *cheddar*
     b. Tropic.carni.webs.inf95.10p: community collection of mammal food webs with documented trophic links augmented with potential prey identified by cosine similarity approaches with only predators in Order of Carnivora and food webs contain at least three predator species, in format prepared for food web metric measurements with R package *cheddar*
 
2. [Predator_prey_pairs_for_dietary_niche_space.csv](Data/Predator_prey_pairs_for_dietary_niche_space.csv) 
     a. Regional (Region) co-occurring pair of predator (trophic.links.consumer) and prey (trophic.links.resource) across communities, with potential prey identified by cosine similarity approach marked as undocumented (trophic.links.documeted = 0). This data is used to quantify hypervolume of predator-specific dietary niche breadth as input data for [Predator_dietary_niche_hypervolume.py](Code/Predator_dietary_niche_hypervolume.py)   

3. [Gower_11traits_PCoA_Axis13_inf95_10p_500g.csv](Data/Gower_11traits_PCoA_Axis13_inf95_10p_500g.csv) 
     a. Thirteen PCoA axes (Axis.1 - Axis.13) based on the Gower dissimilarity values of 11 traits to construct the global multidimensional trait space the total 423 prey species (Species) from both documented prey and potential prey identified by cosine similarity >= 0.95 (in95) and 10% body mass differences (10p), as input data for [Predator_dietary_niche_hypervolume.py](Code/Predator_dietary_niche_hypervolume.py)   

4. [Regional_predator_dietary_hypervolume.csv](Data/Regional_predator_dietary_hypervolume.csv)
    a. Global (region = all) and regional (region = Neotropic, Afrotropic, or Indomalayan) prey trait space hypervolume (hypervolume)

5. [Predator_dietary_niche_hypervolume.csv](Data/Predator_dietary_niche_hypervolume.csv) 
     a. Predator(consumer)-specific dietary niche breadth (hypervolume) based on prey trait space (number of prey species in regional meta-web; n.species) by region 

6. [Predator_dietary_niche_overlap_hypervolume.csv](Data/Predator_dietary_niche_overlap_hypervolume.csv) 
     a. Pairwise dietary niche overlap (intersect_hypervolume/union_hypervolume = overlap_hypervolume) between predator pairs (consumer1, consumer2) based on prey trait space by region

7. [Foodweb_metrics.csv](Data/Foodweb_metrics.csv) 
     i. Food web metric values with documented trophic links only: number of food web species richness (N.Sp.doc: total richness; N.Sp.c.doc: predator richness; N.Sp.r.doc: prey richness; N.Sp.basal.doc: basal prey richness) and food web structure (Per.basal.doc: percentage of basal prey species; Per.inter.doc: percentage of intermediate predator species; Per.top.doc: percentage of top predator species; N.links.doc: total links), and those augmented with potential prey (same variable labeled with inf95.10p). 
     ii. Food web metrics used for formal analyses are all based on documented links augmented with potential prey, include predator-basal prey ratio (cr.ratio), log10-transformed mean predator mass in g (mass.log.c), log10-transformed mean prey mass in g	(mass.log.r). 
     iii. Mean normalized predator generality (generality.sp.normalized) is characterized by mean predator generality (generality.sp) and linkage density (link) based on binary link, as well as those based on probabilistic links with body size model [Rohr et al 2010 - Am Nat](https://doi.org/10.1086/653667) (labeled with *body.size.model.p*)and allometric niche model [Williams et al 2010. PLOS One](https://doi.org/10.1371/journal.pone.0012092) (labeled with *allometric.model.p*) . 
     iv. Mean normalized prey vulnerability (vulnerability.sp.normalized) is characterized by mean prey vulnerability (vulnerability.sp) and linkage density (link) based on binary link, as well as those based on probabilistic links with body size model [Rohr et al 2010 - Am Nat](https://doi.org/10.1086/653667) (labeled with *body.size.model.p*) and allometric niche model [Williams et al 2010. PLOS One](https://doi.org/10.1371/journal.pone.0012092) (labeled with *allometric.model.p*). 
     vi. Mean predator dietary niche breadth (resource.hv.alpha.m.inf) and mean predator dietary niche overlap (resource.hv.beta.m.inf)
     
### Food web lineage history
1. [Foodweb_lineage_history.csv](Data/Foodweb_lineage_history.csv) 
    a. Lineage history of present natural community (naturalpresent_comm_meanDR) and lineage structure of contemporary food web (foodweb_meanDR)
