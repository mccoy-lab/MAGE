# GeoVar Figure Generation

The scripts in this directory generate the following figures in the main MAGE manuscript:

* Figure 5A
* Extended data Figure 4
* Extended data Figure 5
* Extended data Figure 6

## Running Figure Generation

To create the appropriate figures, we recommend the following steps:

1. Create the conda environment:

```
conda env create -f env.yaml
conda activate mage-geovar
```

2. Run the snakemake pipeline for data processing and generation of intermediate tables
```
snakemake -s afs.smk
```

3. Run the primary scripts for table generation: 

```
python3 finemapping_metaCS_dapg_mage_replicate.py
```

NOTE: this step specifically creates our meta credible sets to allow comparison between GTEx multi-tissue results and MAGE results.  

4. Create figures 

```
python3 finemapping_metaCS_geovar_plots.py
```

## Description of Intermediate files

### Allele Frequency Tables

The files detailing allele frequencies have the following columns: 
- `CHR`: chromosome in hg38 
- `SNP`: variant ID in CHROM:POS:REF:ALT
- `A1`: major allele
- `A2`: minor allele
- `MAC`: counts of the minor allele in MAGE 
- `MAF`: minor allele frequency (global)
- `AFR`: minor allele frequency in the AFR ancestry group 
- `EUR`: minor allele frequency in the EUR ancestry group 
- `SAS`: minor allele frequency in the SAS ancestry group 
- `EAS`: minor allele frequency in the EAS ancestry group
- `AMR`: minor allele frequency in the AMR ancestry group

### DAPG + MAGE meta credible set tables 

These are the tables created by `finemapping_metaCS_dapg_mage_replicate.py` and highlight the overlap (or non-overlap) of DAPG credible sets from GTEx and the SuSiE credible sets in MAGE for eQTLs / sQTLs. They have the following columns: 

- `canonicalGene`: Ensembl gene ID 
- `n_MAGE_intersect`: Number of 
- `GTEx_ID`: ID of the variant in GTEx
- `metaCS`: index of the GTEx metaCS  
- `nsnps_metaCS`: number of SNPs in the metaCS
- `tissues`: dictionary of tissues that the eQTL in GTEx is found in
- `geneID`: Full ensembl ID for the gene
- `variantID`: variant ID within MAGE (in chrom:pos:ref:alt format on hg38)
- `variantPIP`: posterior inclusion probability from SuSiE for MAGE eQTL variant
- `variantCredibleSet`: indicator of MAGE eQTL variant corresponding credible set from SuSiE
- `numVariants`: number of variants in the SuSiE 
- `coverage`: coverage of the SuSiE credible set for this eQTL in MAGE, or the probability of a credible set containing a true causal signal 
- `bestHit`: whether the variant is the "lead variant" 
- `n_tissues`: number of tissues in GTEx the focal variant is called as an eQTL
