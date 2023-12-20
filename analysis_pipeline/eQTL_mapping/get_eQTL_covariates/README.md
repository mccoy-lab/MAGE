# Preparing covariates for eQTL mapping

This directory contains code to replicate the eQTL mapping covariate preparation steps done for the initial MAGE publication. This procedure comprises three steps:
1. Calculation of genotype PCs with [Plink](https://www.cog-genomics.org/plink/)
2. Calculation of PEER factors with [peertool](https://github.com/PMBio/peer)
3. Collection and formatting of all eQTL mapping covariates

*Note: these steps should be run in order*<br><br>

## Calculation of top genotype PCs

The `01_get_genotype_PCs.sh` script will calculate the top genotype PCs with [plink](https://www.cog-genomics.org/plink/). We only use variants with MAF > 0.01 for calculation of genotype PCs. We only used 5 PCs as covariates for eQTL mapping, but the script calculates the top 20 for analysis.

The script takes three arguments:
1. `inVCF`: A VCF file containing variants from only the autosomes. This is the `autosomal.MAGE.v0.1.vcf.gz` file created in the [`subset_VCFs`](../../data_preparation/subset_VCFs/) step
2. `outDir`: The directory to write plink output files to
3. `outPrefix`: Prefix of output plink files in output directory
<br><br>

## Calculation of PEER factors

The `02_get_PEER_covariates.sh` script will calculate 60 PEER factors from the filtered expression quantifications used for eQTL mapping. The number of PEER factors to calculate was based on optimizations previously done by GTEx for this sample size.

The script takes two arguments:
1. `invNormTMMBed`: A BED file containing inverse normal transformed TMM values used as input for eQTL mapping. This is the `<prefix>.filtered.inverse-normal-tmm.bed.gz` file created in the [`expression_quantification`](../../data_preparation/expression_quantification/) step
2. `outDir`: The directory to write peertool output files to

### Software requirements

This script assumes you have the `peertool` executable (installation instructions [here](https://github.com/PMBio/peer/wiki/Installation-instructions)) installed and on your PATH.<br><br>

## Collection and formatting of covariates for eQTL mapping

The `03_collect_covariates.sh` script will collect  the genotype PCs and PEER factors calculated in the previous steps and any other covariates and format them for eQTL mapping with FastQTL.

The script takes four arguments:
1. `pcaEigenvec`: The `<prefix>.eigenvec` file created by the `01_get_genotype_PCs.sh` script
2. `peerX`: The `X.csv` file created by the `02_get_PEER_covariates.sh` script
3. `addCovFile`: A tab-separated file with additional covariates to be included in the final covariate file. There should be one line per sample, one column per covariate. The first column should be the sample IDs. For the initial MAGE publication, we included only sample sex as an additional covariate (i.e. in addition to genotype PCs and PEER factors)
4. `outFile`: The formatted output file, to be used for eQTL mapping with FastQTL
<br><br>
