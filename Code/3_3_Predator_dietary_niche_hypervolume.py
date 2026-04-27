
import pandas as pd
from pygmo import hypervolume
import numpy as np
from itertools import combinations

#trophic interactions by realm
df_metaweb = pd.read_csv('Predator_prey_pairs_for_dietary_niche_space.csv')
#PCoA axes to quantity hypervolume of predator dietary niche breadth based on prey trait space 
df_pcoa = pd.read_csv('Gower_11traits_PCoA_Axis13_inf95_10p_500g.csv', index_col=0)

df_pcoa.drop(columns=["Species"], inplace=True)
df_pcoa = df_pcoa - df_pcoa.max(axis=0)

regions = df_metaweb["Region"].unique()

res = []
for region in regions:

    df_region_metaweb = df_metaweb.loc[df_metaweb["Region"] == region]
    consumers = df_region_metaweb["trophic.links.consumer"].unique()

    for consumer in consumers:
        df_consumer_metaweb = df_region_metaweb.loc[df_region_metaweb["trophic.links.consumer"] == consumer]
        resources = df_consumer_metaweb["trophic.links.resource"].unique()
        resources_pcoa = df_pcoa.loc[resources, :].to_numpy()

        # hypervolume
        ref_point = [0 for _ in range(resources_pcoa.shape[1])]

        hv = hypervolume(resources_pcoa)
        hyperv = hv.compute(ref_point)

        n_prey = resources_pcoa.shape[0]

        # results
        res.append((region, consumer, hyperv, resources_pcoa.shape[0]))

df_hypervolume = pd.DataFrame(res, columns=['region', 'consumer', 'hypervolume', 'n.species'])
df_hypervolume.to_csv("Predator_dietary_niche_hypervolume.csv")


# for overlap hypervolume
res = []
for region in regions:

    df_region_metaweb = df_metaweb.loc[df_metaweb["Region"] == region]
    consumers = df_region_metaweb["trophic.links.consumer"].unique()

    for (consumer1, consumer2) in combinations(consumers, 2):
        consumer1_metaweb = df_region_metaweb.loc[df_region_metaweb["trophic.links.consumer"] == consumer1]
        resources1 = consumer1_metaweb["trophic.links.resource"].unique()
        

        consumer2_metaweb = df_region_metaweb.loc[df_region_metaweb["trophic.links.consumer"] == consumer2]
        resources2 = consumer2_metaweb["trophic.links.resource"].unique()

        # union hypervolume
        union_resources = list(set(resources1).union(set(resources2)))
        union_resources_pcoa = df_pcoa.loc[union_resources, :].to_numpy()
        ref_point = [0 for _ in range(union_resources_pcoa.shape[1])]
        hv_union = hypervolume(union_resources_pcoa)
        hyperv_union = hv_union.compute(ref_point)

        # intersect hypervolume
        intersect_resources = list(set(resources1).intersection(set(resources2)))
        if len(intersect_resources) == 0:
            hyperv_intersect = 0
        else:
            intersect_resources_pcoa = df_pcoa.loc[intersect_resources, :].to_numpy()
            ref_point = [0 for _ in range(intersect_resources_pcoa.shape[1])]
            hv_intersect = hypervolume(intersect_resources_pcoa)
            hyperv_intersect = hv_intersect.compute(ref_point)
                    
        # results
        res.append((region, consumer1, consumer2, hyperv_intersect, hyperv_union, hyperv_intersect / hyperv_union))

df_hypervolume = pd.DataFrame(res, columns=['region', 'consumer1', 'consumer2', 'intersect_hypervolume', 'union_hypervolume', 'overlap_hypervolume'])
df_hypervolume.to_csv("Predator_dietary_niche_overlap_hypervolume.csv")

#Integrate these measures with community checklist from Rowan et al. 2020 PNAS https://www.pnas.org/doi/10.1073/pnas.1910489116#supplementary-materials
#to measure mean dietary niche breadth per predator and mean dietary niche overlap among predators by communities
