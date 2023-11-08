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
import matplotlib.pyplot as plt
from geovar import *
from tqdm import tqdm
import gzip
from scipy.stats import mannwhitneyu
from arjun_plot.utils import *

plt.rcParams.update({'font.family': 'sans-serif'})
plt.rcParams.update({'font.sans-serif': 'Helvetica'})

# %matplotlib inline 
# -

# Reading in the replicates vs. non-replicators 
kgpex_dapg_replicates_top_eVariants_df = pd.read_csv('kgpex_dapg_metaCS_replicates_pip95.top_eVariants.092623.tsv.gz', sep="\t")
kgpex_dapg_nonreplicates_top_eVariants_df = pd.read_csv('kgpex_dapg_metaCS_nonreplicates_pip95.top_eVariants.092623.tsv.gz', sep="\t")
nrep = kgpex_dapg_replicates_top_eVariants_df.shape[0]
nnonrep = kgpex_dapg_nonreplicates_top_eVariants_df.shape[0]


finemap_afs_df = pd.read_csv('results/af_tables/finemapped_kgpex.eqtl.tsv', sep="\t")
nominal_afs_df = pd.read_csv('results/af_tables/nominal_kgpex.eqtl.tsv', sep="\t")
finemap_sqtl_afs_df = pd.read_csv('results/af_tables/finemapped_single_kgpex.sqtl.tsv', sep="\t")
afs_df = pd.concat([finemap_afs_df, nominal_afs_df, finemap_sqtl_afs_df]).drop_duplicates()

# Isolating the AFS 
afs_replicates_df = afs_df[afs_df.SNP.isin(kgpex_dapg_replicates_top_eVariants_df.variantID)]
afs_nonreplicates_df = afs_df[afs_df.SNP.isin(kgpex_dapg_nonreplicates_top_eVariants_df.variantID)]
afs_replicates_df.to_csv('/tmp/kgpex_top_evariants_rep.092623.tsv', sep="\t", index=None)
afs_nonreplicates_df.to_csv('/tmp/kgpex_top_evariants_nonrep.092623.tsv', sep="\t", index=None)

geovar_kgpex_eVariants_rep = GeoVar()
geovar_kgpex_eVariants_rep.add_freq_mat(freq_mat_file='/tmp/kgpex_top_evariants_rep.092623.tsv')
geovar_kgpex_eVariants_rep.generate_bins(bins=[(0, 0), (0.0, 0.05), (0.05, 1.0)])
geovar_kgpex_eVariants_rep.geovar_binning()

geovar_kgpex_eVariants_nonrep = GeoVar()
geovar_kgpex_eVariants_nonrep.add_freq_mat(freq_mat_file='/tmp/kgpex_top_evariants_nonrep.092623.tsv')
geovar_kgpex_eVariants_nonrep.generate_bins(bins=[(0, 0), (0.0, 0.05), (0.05, 1.0)])
geovar_kgpex_eVariants_nonrep.geovar_binning()

# The actual plotting code here ...
geovar_plot_eVariants_rep = GeoVarPlot()
geovar_plot_eVariants_rep.add_data_geovar(geovar_kgpex_eVariants_rep)
geovar_plot_eVariants_rep.filter_data(max_freq=0.0005)
geovar_plot_eVariants_rep.add_cmap(str_labels=['U','R', 'C'], lbl_colors=['black', 'black','white'])
geovar_plot_eVariants_rep.reorder_pops(np.array(['AFR','EUR','SAS','EAS','AMR']))

geovar_plot_eVariants_nonrep = GeoVarPlot()
geovar_plot_eVariants_nonrep.add_data_geovar(geovar_kgpex_eVariants_nonrep)
geovar_plot_eVariants_nonrep.filter_data(max_freq=0.0005)
geovar_plot_eVariants_nonrep.add_cmap(str_labels=['U','R','C'], lbl_colors=['black', 'black','white'])
geovar_plot_eVariants_nonrep.reorder_pops(np.array(['AFR','EUR','SAS','EAS','AMR']))

fig, ax = plt.subplots(1, 2, figsize=(6,6), sharey=True, layout='constrained')
# The full plotting routine here ... 
geovar_plot_eVariants_rep.plot_geovar(ax[0]);
geovar_plot_eVariants_nonrep.plot_geovar(ax[1]);

ax[0].set_xticklabels(geovar_plot_eVariants_rep.poplist);
ax[1].set_xticklabels(geovar_plot_eVariants_nonrep.poplist);
ax[1].set_ylabel('')
ax[0].set_title(f'eVariants replicated GTeX\n(n={nrep})', fontsize=11)
ax[1].set_title(f'eVariants non-replicated GTeX\n(n={nnonrep})', fontsize=11)
ax[0].tick_params(axis='y', which='minor', labelsize=9)
ax[1].tick_params(axis='y', which='minor', labelsize=9)
ax[0].tick_params(axis='x', which='minor', labelsize=11)
ax[1].tick_params(axis='x', which='minor', labelsize=11)
ax[0].set_ylabel('Cumulative fraction of variants', fontsize=11)

plt.savefig('kgpex_geovar_rep_nonrep_GTExDAPG.102623.pdf', dpi=300, bbox_inches='tight')
plt.savefig('kgpex_geovar_rep_nonrep_GTExDAPG.102623.png', dpi=300, bbox_inches='tight')
# -

# Tissue Plot for DAPG replicating vs non-replicating 
# Reading in the replicates vs. non-replicators 
rep_dapg_kgpex_top_eVariants_df = pd.read_csv('dapg_kgpex_metaCS_replicates_pip95.top_eVariants.092623.tsv.gz', sep="\t")
nonrep_dapg_kgpex_top_eVariants_df = pd.read_csv('dapg_kgpex_metaCS_nonreplicates_pip95.top_eVariants.092623.tsv.gz', sep="\t")
print("rep / non-rep:", rep_dapg_kgpex_top_eVariants_df.shape[0], nonrep_dapg_kgpex_top_eVariants_df.shape[0])
fig, ax = plt.subplots(1,1,figsize=(4,4), layout='constrained')
ax.hist(nonrep_dapg_kgpex_top_eVariants_df.n_tissues, bins=25, density=True, alpha=0.5, label='Non-replicating top eVariants\n(GTEx to MAGE)')
ax.hist(rep_dapg_kgpex_top_eVariants_df.n_tissues, bins=25, density=True, alpha=0.5, label='Replicating top eVariants\n(GTEx to MAGE)')
ax.legend(fontsize=8, frameon=False)
ax.set_xlabel('Tissues in cross-tissue \nmerged credible sets (GTEx)', fontsize=11)
ax.set_ylabel(r'Density', fontsize=11)
ax.tick_params(axis='both', which='major', labelsize=9)
debox(ax)
plt.savefig('dapg_non_replicate_num_tissues.103023.pdf', dpi=300, bbox_inches='tight')
plt.savefig('dapg_non_replicate_num_tissues.103023.png', dpi=300, bbox_inches='tight')
print(mannwhitneyu(rep_dapg_kgpex_top_eVariants_df.n_tissues, nonrep_dapg_kgpex_top_eVariants_df.n_tissues))

# +
# Supplementary plots of restricted GeoVar of variants
afs_replicates_df = pd.concat([afs_replicates_df[afs_replicates_df.SNP.isin([x])] for x in tqdm(kgpex_dapg_replicates_top_eVariants_df.variantID)])
afs_nonreplicates_df = pd.concat([afs_nonreplicates_df[afs_nonreplicates_df.SNP.isin([x])] for x in tqdm(kgpex_dapg_nonreplicates_top_eVariants_df.variantID)])

tot_afs_df = pd.concat([afs_replicates_df, afs_nonreplicates_df])
non_globally_common_df = tot_afs_df.iloc[~np.all(tot_afs_df[['AFR', 'EUR', 'SAS', 'EAS', 'AMR']].values >= 0.05, axis=1), :]
eur_rare_df = tot_afs_df.iloc[tot_afs_df['EUR'].values <= 0.00, :]
eur_afr_rare_df = tot_afs_df.iloc[np.all(tot_afs_df[['EUR', 'AFR']].values <= 0.00, axis=1), :]

tot_afs_df.to_csv('/tmp/kgpex_tot_variants.092623.tsv', sep="\t", index=None)
non_globally_common_df.to_csv('/tmp/kgpex_nonglobal.092623.tsv', sep="\t", index=None)
eur_rare_df.to_csv('/tmp/kgpex_eur_rare.092623.tsv', sep="\t", index=None)
eur_afr_rare_df.to_csv('/tmp/kgpex_eur_afr_rare.092623.tsv', sep="\t", index=None)

# +
geovar_kgpex_eVariants_tot = GeoVar()
geovar_kgpex_eVariants_tot.add_freq_mat(freq_mat_file='/tmp/kgpex_tot_variants.092623.tsv')

geovar_kgpex_eVariants_noglobal = GeoVar()
geovar_kgpex_eVariants_noglobal.add_freq_mat(freq_mat_file='/tmp/kgpex_nonglobal.092623.tsv')

geovar_kgpex_eVariants_eur_rare = GeoVar()
geovar_kgpex_eVariants_eur_rare.add_freq_mat(freq_mat_file='/tmp/kgpex_eur_rare.092623.tsv')

geovar_kgpex_eVariants_eur_afr_rare = GeoVar()
geovar_kgpex_eVariants_eur_afr_rare.add_freq_mat(freq_mat_file='/tmp/kgpex_eur_afr_rare.092623.tsv')

for g in [geovar_kgpex_eVariants_tot,geovar_kgpex_eVariants_noglobal, geovar_kgpex_eVariants_eur_rare, geovar_kgpex_eVariants_eur_afr_rare]:
    g.generate_bins(bins=[(0, 0), (0.0, 0.05), (0.05, 1.0)])
    g.geovar_binning()

# The actual plotting code here ...
geovar_plot_eVariants_tot = GeoVarPlot()
geovar_plot_eVariants_tot.add_data_geovar(geovar_kgpex_eVariants_rep)

geovar_plot_eVariants_noglobal = GeoVarPlot()
geovar_plot_eVariants_noglobal.add_data_geovar(geovar_kgpex_eVariants_noglobal)

geovar_plot_eVariants_eur_rare = GeoVarPlot()
geovar_plot_eVariants_eur_rare.add_data_geovar(geovar_kgpex_eVariants_eur_rare)

geovar_plot_eVariants_eur_afr_rare = GeoVarPlot()
geovar_plot_eVariants_eur_afr_rare.add_data_geovar(geovar_kgpex_eVariants_eur_afr_rare)

for g in [geovar_plot_eVariants_tot,geovar_plot_eVariants_noglobal, geovar_plot_eVariants_eur_rare, geovar_plot_eVariants_eur_afr_rare]:
# These can be repeated in loop
    g.filter_data(max_freq=0.0005)
    g.add_cmap(str_labels=['U','R', 'C'], lbl_colors=['black', 'black','white'])
    g.reorder_pops(np.array(['AFR','EUR','SAS','EAS','AMR']))


fig, ax = plt.subplots(1, 4, figsize=(8.5,5), layout='constrained', sharey=True)
# The full plotting routine here ... 
geovar_plot_eVariants_tot.plot_geovar(ax[0]);
geovar_plot_eVariants_noglobal.plot_geovar(ax[1]);
geovar_plot_eVariants_eur_rare.plot_geovar(ax[2]);
geovar_plot_eVariants_eur_afr_rare.plot_geovar(ax[3]);


ax[0].set_xticklabels(geovar_plot_eVariants_tot.poplist);
ax[1].set_xticklabels(geovar_plot_eVariants_noglobal.poplist);

ax[2].set_xticklabels(geovar_plot_eVariants_eur_rare.poplist);
ax[3].set_xticklabels(geovar_plot_eVariants_eur_afr_rare.poplist);

ax[1].set_ylabel('')
ax[2].set_ylabel('')
ax[3].set_ylabel('')
ax[0].set_title(f'lead eVariants\n(n={15664})', fontsize=11)
ax[1].set_title(f'lead eVariants\n(not globally common)\n(n={non_globally_common_df.shape[0]})', fontsize=11)
ax[2].set_title(f'lead eVariants\n(unobserved in EUR)\n(n={eur_rare_df.shape[0]})', fontsize=11)
ax[3].set_title(f'lead eVariants\n(unobserved in EUR + AFR)\n(n={eur_afr_rare_df.shape[0]})', fontsize=11)

ax[0].tick_params(axis='y', which='minor', labelsize=9)
ax[1].tick_params(axis='y', which='minor', labelsize=9)
ax[0].tick_params(axis='x', which='minor', labelsize=11)
ax[1].tick_params(axis='x', which='minor', labelsize=11)
ax[0].set_ylabel('Cumulative fraction of variants', fontsize=11)
label_multipanel(ax, ['A', 'B', 'C', 'D'], fontweight='bold')
plt.savefig('supp_fig_draft_geovar.unobserved.102623.pdf', bbox_inches='tight')

kgpex_dapg_replicates_cs_df = pd.read_csv('kgpex_dapg_metaCS_replicates_pip95.092623.tsv.gz', sep="\t")
kgpex_dapg_nonreplicates_cs_df = pd.read_csv('kgpex_dapg_metaCS_nonreplicates_pip95.092623.tsv.gz', sep="\t")
kgpex_dapg_replicates_cs_df['replicate'] = 1.0
kgpex_dapg_nonreplicates_cs_df['replicate'] = 0.0
full_kgpex_dapg_cs_df = pd.concat([kgpex_dapg_replicates_cs_df, kgpex_dapg_nonreplicates_cs_df])
full_kgpex_dapg_cs_df

dapg_kgpex_replicates_cs_df = pd.read_csv('dapg_kgpex_metaCS_replicates_pip95.092623.tsv.gz', sep="\t")
dapg_kgpex_nonreplicates_cs_df = pd.read_csv('dapg_kgpex_metaCS_nonreplicates_pip95.092623.tsv.gz', sep="\t")
full_dapg_kgpex_cs_df = pd.concat([dapg_kgpex_replicates_cs_df, dapg_kgpex_nonreplicates_cs_df])

# +
# number of MAGE credible sets that replicate vs those that do not
print(full_kgpex_dapg_cs_df[full_kgpex_dapg_cs_df['replicate'] == 1.0][['canonicalGene', 'variantCredibleSet']].drop_duplicates().shape[0])
print(full_kgpex_dapg_cs_df[full_kgpex_dapg_cs_df['replicate'] == 0.0][['canonicalGene', 'variantCredibleSet']].drop_duplicates().shape[0])

# Number of genes involved in replicating vs non-replicating 
print(full_kgpex_dapg_cs_df[full_kgpex_dapg_cs_df['replicate'] == 1.0][['canonicalGene']].drop_duplicates().shape[0])
print(full_kgpex_dapg_cs_df[full_kgpex_dapg_cs_df['replicate'] == 0.0][['canonicalGene']].drop_duplicates().shape[0])
# -

# number of genes with > 0 CS and none found in GTEx
full_kgpex_dapg_cs_df[~full_kgpex_dapg_cs_df.canonicalGene.isin(full_dapg_kgpex_cs_df.canonicalGene)]['canonicalGene'].drop_duplicates()

# Number of CS only found in GTEx & number of genes affected
print(full_dapg_kgpex_cs_df[full_dapg_kgpex_cs_df.n_KGPEx_intersect == 0][['canonicalGene', 'metaCS']].drop_duplicates().shape[0])
print(full_dapg_kgpex_cs_df[full_dapg_kgpex_cs_df.n_KGPEx_intersect == 0]['canonicalGene'].drop_duplicates().shape[0])

full_dapg_kgpex_cs_df['canonicalGene'].drop_duplicates()

nonrep_dapg_kgpex_top_eVariants_df.shape

rep_dapg_kgpex_top_eVariants_df.shape

# Plots for splicing QTLs
kgpex_splice_qtl_fp = '/scratch16/rmccoy22/1KGP_expression/QTL_analysis/sQTL_discovery/sQTL_mapping_results/susie_results/sQTL_finemapping.gene.mergedCredibleSets.bestHits.txt.gz'
kgpex_splice_qtl_top_sVariants_df = pd.read_csv(kgpex_splice_qtl_fp, sep="\t")
kgpex_splice_qtl_top_sVariants_df

afs_sQTLs_raw_df = afs_df[afs_df.SNP.isin(kgpex_splice_qtl_top_sVariants_df.variantID)]
afs_sQTLs_df = pd.concat([afs_sQTLs_raw_df[afs_sQTLs_raw_df.SNP.isin([x])] for x in tqdm(kgpex_splice_qtl_top_sVariants_df.variantID)])

# +
# Supplementary plots of restricted GeoVar of variants 
non_globally_common_sqtl_df = afs_sQTLs_df.iloc[~np.all(afs_sQTLs_df[['AFR', 'EUR', 'SAS', 'EAS', 'AMR']].values >= 0.05, axis=1), :]
eur_rare_sqtl_df = afs_sQTLs_df.iloc[afs_sQTLs_df['EUR'].values <= 0.00, :]
eur_afr_rare_sqtl_df = afs_sQTLs_df.iloc[np.all(afs_sQTLs_df[['EUR', 'AFR']].values <= 0.00, axis=1), :]

afs_sQTLs_df.to_csv('/tmp/kgpex_tot_sqtl_variants.092623.tsv', sep="\t", index=None)
non_globally_common_sqtl_df.to_csv('/tmp/kgpex_nonglobal_sqtl.092623.tsv', sep="\t", index=None)
eur_rare_sqtl_df.to_csv('/tmp/kgpex_eur_rare_sqtl.092623.tsv', sep="\t", index=None)
eur_afr_rare_sqtl_df.to_csv('/tmp/kgpex_eur_afr_rare_sqtl.092623.tsv', sep="\t", index=None)

geovar_kgpex_sVariants_tot = GeoVar()
geovar_kgpex_sVariants_tot.add_freq_mat(freq_mat_file='/tmp/kgpex_tot_sqtl_variants.092623.tsv')

geovar_kgpex_sVariants_noglobal = GeoVar()
geovar_kgpex_sVariants_noglobal.add_freq_mat(freq_mat_file='/tmp/kgpex_nonglobal_sqtl.092623.tsv')

geovar_kgpex_sVariants_eur_rare = GeoVar()
geovar_kgpex_sVariants_eur_rare.add_freq_mat(freq_mat_file='/tmp/kgpex_eur_rare_sqtl.092623.tsv')

geovar_kgpex_sVariants_eur_afr_rare = GeoVar()
geovar_kgpex_sVariants_eur_afr_rare.add_freq_mat(freq_mat_file='/tmp/kgpex_eur_afr_rare_sqtl.092623.tsv')

for g in [geovar_kgpex_sVariants_tot,geovar_kgpex_sVariants_noglobal, geovar_kgpex_sVariants_eur_rare, geovar_kgpex_sVariants_eur_afr_rare]:
    g.generate_bins(bins=[(0, 0), (0.0, 0.05), (0.05, 1.0)])
    g.geovar_binning()

# The actual plotting code here ...
geovar_plot_sVariants_tot = GeoVarPlot()
geovar_plot_sVariants_tot.add_data_geovar(geovar_kgpex_sVariants_tot)

geovar_plot_sVariants_noglobal = GeoVarPlot()
geovar_plot_sVariants_noglobal.add_data_geovar(geovar_kgpex_sVariants_noglobal)

geovar_plot_sVariants_eur_rare = GeoVarPlot()
geovar_plot_sVariants_eur_rare.add_data_geovar(geovar_kgpex_sVariants_eur_rare)

geovar_plot_sVariants_eur_afr_rare = GeoVarPlot()
geovar_plot_sVariants_eur_afr_rare.add_data_geovar(geovar_kgpex_sVariants_eur_afr_rare)

for g in [geovar_plot_sVariants_tot,geovar_plot_sVariants_noglobal, geovar_plot_sVariants_eur_rare, geovar_plot_sVariants_eur_afr_rare]:
# These can be repeated in loop
    g.filter_data(max_freq=0.0005)
    g.add_cmap(str_labels=['U','R', 'C'], lbl_colors=['black', 'black','white'])
    g.reorder_pops(np.array(['AFR','EUR','SAS','EAS','AMR']))


fig, ax = plt.subplots(1, 4, figsize=(8.5,5), layout='constrained', sharey=True)
# The full plotting routine here ... 
geovar_plot_sVariants_tot.plot_geovar(ax[0]);
geovar_plot_sVariants_noglobal.plot_geovar(ax[1]);
geovar_plot_sVariants_eur_rare.plot_geovar(ax[2]);
geovar_plot_sVariants_eur_afr_rare.plot_geovar(ax[3]);


ax[0].set_xticklabels(geovar_plot_sVariants_tot.poplist);
ax[1].set_xticklabels(geovar_plot_sVariants_noglobal.poplist);

ax[2].set_xticklabels(geovar_plot_sVariants_eur_rare.poplist);
ax[3].set_xticklabels(geovar_plot_sVariants_eur_afr_rare.poplist);

ax[1].set_ylabel('')
ax[2].set_ylabel('')
ax[3].set_ylabel('')
ax[0].set_title(f'lead sVariants\n(n={afs_sQTLs_df.shape[0]})', fontsize=11)
ax[1].set_title(f'lead sVariants\n(not globally common)\n(n={non_globally_common_sqtl_df.shape[0]})', fontsize=11)
ax[2].set_title(f'lead sVariants\n(unobserved in EUR)\n(n={eur_rare_sqtl_df.shape[0]})', fontsize=11)
ax[3].set_title(f'lead sVariants\n(unobserved in EUR + AFR)\n(n={eur_afr_rare_sqtl_df.shape[0]})', fontsize=11)

ax[0].tick_params(axis='y', which='minor', labelsize=9)
ax[1].tick_params(axis='y', which='minor', labelsize=9)
ax[0].tick_params(axis='x', which='minor', labelsize=11)
ax[1].tick_params(axis='x', which='minor', labelsize=11)
ax[0].set_ylabel('Cumulative fraction of variants', fontsize=11)
label_multipanel(ax, ['A', 'B', 'C', 'D'], fontweight='bold')
plt.savefig('supp_fig_draft_geovar.sqtls.unobserved.102623.pdf', bbox_inches='tight')