# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
import numpy as np 
import pandas as pd
import scipy 
import matplotlib.pyplot as plt
from tqdm import tqdm
import warnings

with warnings.catch_warnings():
    warnings.simplefilter('ignore')

# %matplotlib inline
# -

# # Goals of additional replication analyses 
#
# * One of the main goals is to find genes with minimal overlap in CS between GTEx and KGPEx
# * Key Q: how do we combine the results over mulitiple tissues for DAPG to create meta-CS?
#
#
# NOTES:
# * Q: for each CS in our data, does it overlap with any CS in DAPG?
#       * take all CS per-tissue in DAPG, and ask what the SNP-level overlap is with KGPEX Susie Credible Set?
#
# 1. KGPEX CS that don't replicate in GTEx 
# Pseudocode:
#
# ```
# for gene in gene_gtex:
#     all_variants = aggregate variants in all DAPG CS (unique variants)
#     for cs in kgpex_CS[gene]:
#         intersection_variants = set(cs, all_variants).intersect()
#         if len(intersection_variants) == 0:
#             # add to list of non-replicated cs...
#         # also store the replicated ones ... 
# ```
#
# * for non-replicating CS we can then extract the top eVariants
#
# 2. GTEx credible sets that don't replicate in our data 
#
#     a. Combine CS in DAPG across tissues? (take set overlap between independent credible sets?)
#     b. You will still end up with sets of variants across Tissues in DAPG
#     c. For metaCS that don't intersect with any CS in KGPEx, move them to special list
#     d. Take top-eVariant across the metaCS (defined as maximal PIP?)
#     e. For max-PIP variant, keep track of cell-type in which it is maximal PIP as well ...  
#

# zgrep -v "chrX" /data/rmccoy22/GTEx_Analysis_v8_eQTL/GTEx_v8_finemapping_DAPG/GTEx_v8_finemapping_DAPG.txt.gz | awk '$4 >= 0.95' > GTEx_v8_finemapping_DAPG.autosome_pip95.txt
# awk 'FNR == NR {a[$3] = $1; next} {OFS="\t";print $0,a[$5]}' /scratch16/rmccoy22/1KGP_expression/sharedData/variant_lists/kgpex_variants.all.with_rsID_GTEx_ID.txt GTEx_v8_finemapping_DAPG.autosome_pip95.txt > GTEx_v8_finemapping_DAPG.autosome_pip95.kgpex_ids.txt 
#
 
#1. Reading in the DAPG Results 
dapg_fp = 'GTEx_v8_finemapping_DAPG.autosome_pip95.kgpex_ids.txt'
dapg_finemap_df = pd.read_csv(dapg_fp, low_memory=True, sep="\t")
dapg_finemap_df["canonicalGene"] = [g.split(".")[0]  for g in dapg_finemap_df.gene_id]
dapg_finemap_df.rename(columns={"gene_id": "geneID", "Unnamed: 6": "variantID"}, inplace=True)
dapg_finemap_df.head()

# zcat /scratch16/rmccoy22/1KGP_expression/QTL_analysis/eQTL_discovery/eQTL_mapping_results/susie_results/eQTL_finemapping.significantAssociations.txt.gz > eQTL_finemapping.significantAssociations.txt
# awk 'FNR == NR {a[$1] = $3; next} {OFS="\t";print $0,a[$2]}' /scratch16/rmccoy22/1KGP_expression/sharedData/variant_lists/kgpex_variants.all.with_rsID_GTEx_ID.txt eQTL_finemapping.significantAssociations.txt >  eQTL_finemapping.significantAssociations.anno.txt

# +
# 2. Reading in the kgpex results and merging appropriately
susie_kgpex_fp = '/scratch16/rmccoy22/1KGP_expression/QTL_analysis/eQTL_discovery/eQTL_mapping_results/susie_results/eQTL_finemapping.credibleSets.txt.gz'
finemapping_res_df = pd.read_csv(susie_kgpex_fp, sep="\t")
# 3. Small renaming necessary prior to merging ... 
finemapping_res_df.rename(columns={"credibleSet": "variantCredibleSet"}, inplace=True)
kgpex_finemapped_df = pd.read_csv('eQTL_finemapping.significantAssociations.anno.txt', sep="\t").rename(columns={'Unnamed: 4': 'GTEx_ID'})

kgpex_finemapped_df["canonicalGene"] = [g.split(".")[0]  for g in kgpex_finemapped_df.geneID]
finemapping_res_df['canonicalGene'] = [g.split(".")[0]  for g in finemapping_res_df.geneID]

# This file defines the "best Hits" used by Dylan for consistency - mainly for breaking ties ... 
core_eVariants_df = pd.read_csv('/scratch16/rmccoy22/1KGP_expression/QTL_analysis/eQTL_discovery/eQTL_mapping_results/susie_results/eQTL_finemapping.bestHits.txt.gz', sep="\t")
core_eVariants_df["canonicalGene"] = [g.split(".")[0]  for g in core_eVariants_df.geneID]
core_eVariants_df['bestHit'] = True

kgpex_replication_gtex_info_df = kgpex_finemapped_df.merge(finemapping_res_df, how='left', on=['canonicalGene', 'geneID', 'variantCredibleSet'])
kgpex_replication_gtex_info_df = kgpex_replication_gtex_info_df.merge(core_eVariants_df, how='left').fillna(value={'bestHit': False})

def isolate_dapg_variants(dapg_df, pip_thresh=0.0):
    """Isolate DAPG variants for a specific gene."""
    assert 'canonicalGene' in dapg_df
    assert 'variant_pip' in dapg_df
    assert 'variant_id' in dapg_df
    return dapg_df[dapg_df.variant_pip > pip_thresh].groupby('canonicalGene')['variant_id'].agg(lambda x: list(x)).to_dict()


def kgpex_set_overlap(kgpex_df, gene_id="ENSG00000000457", compare_variants=[], pip_thresh=0.0, p_thresh=0.0):
    assert 'canonicalGene' in kgpex_df.columns
    assert 'GTEx_ID' in kgpex_df.columns
    assert 'variantCredibleSet' in kgpex_df.columns
    non_replicated_vars = []
    replicated_vars = []
    cur_df = kgpex_df[(kgpex_df.canonicalGene == gene_id) & (kgpex_df.variantPIP > pip_thresh)].groupby('variantCredibleSet')['GTEx_ID'].agg(lambda x: set(x)).reset_index()
    for cs_label, cs_ids in zip(cur_df.variantCredibleSet, cur_df.GTEx_ID):
        n_intersect = len(cs_ids.intersection(compare_variants))
        prop_intersect = n_intersect / len(cs_ids)
        if prop_intersect <= p_thresh:
            non_replicated_vars.append([gene_id, cs_label, n_intersect, list(cs_ids)])
        else:
            replicated_vars.append([gene_id, cs_label, n_intersect, list(cs_ids)])
    rep_df = pd.DataFrame(replicated_vars, columns=['canonicalGene', 'variantCredibleSet', 'n_DAPG_intersect', 'variants'])
    non_rep_df = pd.DataFrame(non_replicated_vars, columns=['canonicalGene', 'variantCredibleSet', 'n_DAPG_intersect', 'variants'])
    return rep_df, non_rep_df 

# The first function implements a dictionary for the total variants found via DAPG across all tissues...
dapg_variant_dict = isolate_dapg_variants(dapg_finemap_df)

uniq_genes = kgpex_replication_gtex_info_df.canonicalGene.unique()
rep_df_list = []
nonrep_df_list = []
for gene in tqdm(uniq_genes):
    try:
        rep_df, nonrep_df = kgpex_set_overlap(kgpex_replication_gtex_info_df, gene_id=gene, compare_variants=dapg_variant_dict[gene])
    except KeyError:
        rep_df, nonrep_df = kgpex_set_overlap(kgpex_replication_gtex_info_df, gene_id=gene, compare_variants=[])
    rep_df_list.append(rep_df)
    nonrep_df_list.append(nonrep_df)

total_rep_df = pd.concat(rep_df_list)
total_nonrep_df = pd.concat(nonrep_df_list)

nonrep_kgpex_dapg_df = total_nonrep_df.explode('variants').rename(columns={'variants':'GTEx_ID'}).merge(kgpex_replication_gtex_info_df)
rep_kgpex_dapg_df = total_rep_df.explode('variants').rename(columns={'variants':'GTEx_ID'}).merge(kgpex_replication_gtex_info_df)

# Filter this to the top eVariants in each credible set ...  
nonrep_kgpex_dapg_top_eVariants_df = nonrep_kgpex_dapg_df.sort_values(['coverage', 'variantPIP', 'bestHit'], ascending=False).groupby(['canonicalGene', 'variantCredibleSet']).head(1)
rep_kgpex_dapg_top_eVariants_df = rep_kgpex_dapg_df.sort_values(['coverage', 'variantPIP', 'bestHit'], ascending=False).groupby(['canonicalGene', 'variantCredibleSet']).head(1)

# Some sanity checks that the number of credible sets are preserved 
print(nonrep_kgpex_dapg_df.groupby(['canonicalGene', 'variantCredibleSet']).agg('sum').shape[0], rep_kgpex_dapg_df.groupby(['canonicalGene', 'variantCredibleSet']).agg('sum').shape[0])
print(kgpex_replication_gtex_info_df.groupby(['canonicalGene', 'variantCredibleSet']).agg('sum').shape[0])
print(nonrep_kgpex_dapg_df.groupby(['canonicalGene', 'variantCredibleSet']).agg('sum').shape[0] + rep_kgpex_dapg_df.groupby(['canonicalGene', 'variantCredibleSet']).agg('sum').shape[0])

# check that credible sets are preserved for each of the "top eVariant" files as well
print(rep_kgpex_dapg_top_eVariants_df.groupby(['canonicalGene', 'variantCredibleSet']).agg('sum').shape[0])
print(nonrep_kgpex_dapg_top_eVariants_df.groupby(['canonicalGene', 'variantCredibleSet']).agg('sum').shape[0])

# Verification that the unique number of genes should be the same here?
print(uniq_genes.size)
not_in_dapg = []
for g in uniq_genes:
    try:
        dapg_variant_dict[g]
    except KeyError:
        not_in_dapg.append(g)

nonrep_kgpex_dapg_gene_list = set(total_nonrep_df.canonicalGene.unique())
rep_kgpex_dapg_gene_list = set(total_rep_df.canonicalGene.unique())
intersect_genes = rep_kgpex_dapg_gene_list.intersection(nonrep_kgpex_dapg_gene_list)

rep_kgpex_dapg_df.to_csv('kgpex_dapg_metaCS_replicates_pip95.092623.tsv.gz', sep="\t", index=None, na_rep="NA")
nonrep_kgpex_dapg_df.to_csv('kgpex_dapg_metaCS_nonreplicates_pip95.092623.tsv.gz', sep="\t", index=None, na_rep="NA")

rep_kgpex_dapg_top_eVariants_df.to_csv('kgpex_dapg_metaCS_replicates_pip95.top_eVariants.092623.tsv.gz', sep="\t", index=None, na_rep="NA")
nonrep_kgpex_dapg_top_eVariants_df.to_csv('kgpex_dapg_metaCS_nonreplicates_pip95.top_eVariants.092623.tsv.gz', sep="\t", index=None, na_rep="NA")

# ------  Credible Sets in DAP-G that are not found in KGPEx ------ #

from copy import deepcopy

def merge_sets(list_of_sets):
	pooled = deepcopy(list_of_sets)
	merging = True
	while merging:
		merging=False
		for i,group in enumerate(pooled):
			matching_group = next((g for g in pooled[i+1:] if g.intersection(group)),None)
			if not matching_group: continue
			group.update(matching_group)
			pooled.remove(matching_group)
			merging = True
	return pooled

dapg_finemap_metaCS_draft_df = dapg_finemap_df.groupby(['canonicalGene', 'tissue_id', 'cluster_id'])['variant_id'].agg(lambda x: set(x)).reset_index().rename(columns={'variant_id': 'variant_set'})
dapg_finemap_metaCS_draft_dict = dapg_finemap_metaCS_draft_df.groupby('canonicalGene').agg(lambda x: list(x))[['tissue_id', 'variant_set']].to_dict()
dapg_finemap_metaCS_draft_df.head()


def get_variant_tissue_dict(tissues, variant_sets):    
    tissue_dict = {}
    for t, v in zip(tissues, variant_sets):
        for s in v:
            if s in tissue_dict:
                tissue_dict[s] = tissue_dict[s] + [t]
            else:
                tissue_dict[s] = [t]
    return tissue_dict

def tissue_snp_count(tissue_dict, snp_set = []):
    """Count the number of tissues for which this SNP appears significant in a CS."""
    tissue_counts = {}
    for s in snp_set:
        if s in tissue_dict:
            for a in tissue_dict[s]:
                if a not in tissue_counts:
                    tissue_counts[a] = 1
                else:
                    tissue_counts[a] += 1
    return tissue_counts

# Create the metaCS for DAP-G GTEx Genes (this aggregates CS across tissues to be similar to one another)
dapg_genes = dapg_finemap_metaCS_draft_df.canonicalGene.unique()

dapg_metaCS_dict = {}
dapg_metaCS_tissue_dict = {}
for gene in tqdm(dapg_genes):
    tissue_list = dapg_finemap_metaCS_draft_dict['tissue_id'][gene]
    test_set_list = dapg_finemap_metaCS_draft_dict['variant_set'][gene]
    cur_tissue_dict = get_variant_tissue_dict(tissue_list, test_set_list)
    dapg_metaCS_dict[gene] = merge_sets(test_set_list)
    dapg_metaCS_tissue_dict[gene] = {}
    for i,s in enumerate(dapg_metaCS_dict[gene]):
        dapg_metaCS_tissue_dict[gene][i] = tissue_snp_count(cur_tissue_dict, s)

kgpex_variant_dict = kgpex_replication_gtex_info_df[kgpex_replication_gtex_info_df.variantPIP > 0.0].groupby('canonicalGene')['GTEx_ID'].agg(lambda x: list(x)).to_dict()

rep_dapg_df_list = []
nonrep_dapg_df_list = []

for gene in tqdm(dapg_genes):
    try:
        kgpex_variants = kgpex_variant_dict[gene]
    except:
        kgpex_variants = []
    non_replicated_vars = []
    replicated_vars = []
    for i,s in enumerate(dapg_metaCS_dict[gene]):
        n_intersect = len(s.intersection(kgpex_variants))
        prop_intersect = n_intersect / len(s)
        if prop_intersect <= 0.0:
            non_replicated_vars.append([gene, n_intersect, list(s), i, len(s), dapg_metaCS_tissue_dict[gene][i]])
        else:
            replicated_vars.append([gene, n_intersect, list(s), i, len(s), dapg_metaCS_tissue_dict[gene][i]])        
    rep_df = pd.DataFrame(replicated_vars, columns=['canonicalGene', 'n_KGPEx_intersect', 'variants', 'metaCS', 'nsnps_metaCS', 'tissues'])
    nonrep_df = pd.DataFrame(non_replicated_vars, columns=['canonicalGene', 'n_KGPEx_intersect', 'variants', 'metaCS', 'nsnps_metaCS', 'tissues'])
    rep_dapg_df_list.append(rep_df)
    nonrep_dapg_df_list.append(nonrep_df)

rep_dapg_df = pd.concat(rep_dapg_df_list)
nonrep_dapg_df = pd.concat(nonrep_dapg_df_list)

rep_dapg_kgpex_df = rep_dapg_df.explode('variants').rename(columns={"variants": "GTEx_ID"}).merge(kgpex_replication_gtex_info_df)
nonrep_dapg_kgpex_df = nonrep_dapg_df.explode('variants').rename(columns={"variants": "variant_id"}).merge(dapg_finemap_df)
nonrep_dapg_kgpex_df['n_tissues'] = [len(x) for x in nonrep_dapg_kgpex_df.tissues]
rep_dapg_kgpex_df['n_tissues'] = [len(x) for x in rep_dapg_kgpex_df.tissues]

# NOTE: we filter by DAPG cluster pip and variant_pip for the non-replicating variants ... 
nonrep_dapg_kgpex_top_eVariants_df = nonrep_dapg_kgpex_df.sort_values(['cluster_pip', 'variant_pip'], ascending=False).groupby(['canonicalGene', 'metaCS']).head(1)
rep_dapg_kgpex_top_eVariants_df = rep_dapg_kgpex_df.sort_values(['coverage', 'variantPIP'], ascending=False).groupby(['canonicalGene', 'metaCS']).head(1)

# Do some sanity checks ... 
# tot_num_metaCS = np.sum([len(dapg_metaCS_dict[x]) for x in dapg_metaCS_dict])
# tot_rep_num_dapg_metacs = rep_dapg_kgpex_df.groupby(['canonicalGene', 'metaCS']).agg('sum').shape[0]
# tot_nonrep_num_dapg_metacs = nonrep_dapg_kgpex_df.groupby(['canonicalGene', 'metaCS']).agg('sum').shape[0]
# print(tot_num_metaCS, tot_rep_num_dapg_metacs, tot_nonrep_num_dapg_metacs, tot_rep_num_dapg_metacs+tot_nonrep_num_dapg_metacs)

# Check that the number of metaCS is preserved across top eVariant designations as well ...
# print(nonrep_dapg_kgpex_top_eVariants_df.groupby(['canonicalGene', 'metaCS']).agg('sum').shape[0], tot_nonrep_num_dapg_metacs)
# print(rep_dapg_kgpex_top_eVariants_df.groupby(['canonicalGene', 'metaCS']).agg('sum').shape[0], tot_rep_num_dapg_metacs)

rep_dapg_kgpex_df.to_csv('dapg_kgpex_metaCS_replicates_pip95.092623.tsv.gz', sep="\t", index=None, na_rep="NA")
nonrep_dapg_kgpex_df.to_csv('dapg_kgpex_metaCS_nonreplicates_pip95.092623.tsv.gz', sep="\t", index=None, na_rep="NA")

rep_dapg_kgpex_top_eVariants_df.to_csv('dapg_kgpex_metaCS_replicates_pip95.top_eVariants.092623.tsv.gz', sep="\t", index=None, na_rep="NA")
nonrep_dapg_kgpex_top_eVariants_df.to_csv('dapg_kgpex_metaCS_nonreplicates_pip95.top_eVariants.092623.tsv.gz', sep="\t", index=None, na_rep="NA")