# Preparing covariates for sQTL mapping

This directory contains code to replicate the sQTL mapping covariate preparation steps done for the initial MAGE publication. This procedure comprises three steps:
1. Calculation of PEER factors with [peertool](https://github.com/PMBio/peer)
2. Collection and formatting of all sQTL mapping covariates

*Note: these steps should be run in order, and require that genotype PCs were previously calculated during eQTL mapping preparation*<br><br>

## Calculation of PEER factors

The `01_get_PEER_covariates.sh` script will calculate 15 PEER factors from the filtered intron excision ratios used for sQTL mapping. The number of PEER factors to calculate was based on optimizations previously done by GTEx for this sample size.

The script takes two arguments:
1. `normedRatiosBed`: A BED file containing normalized intron excision ratios used as input for sQTL mapping. This is the `<prefix>_perind.counts.normed.genes.filtered.bed.gz` file created in the [`splicing_quantification`](../../data_preparation/splicing_quantification/) step
2. `outDir`: The directory to write peertool output files to

### Software requirements

This script assumes you have the `peertool` executable (installation instructions [here](https://github.com/PMBio/peer/wiki/Installation-instructions)) installed and on your PATH.<br><br>

## Collection and formatting of covariates for sQTL mapping

The `02_collect_covariates.sh` script will collect  the genotype PCs and PEER factors calculated in the previous steps and any other covariates and format them for sQTL mapping with FastQTL.

The script takes four arguments:
1. `pcaEigenvec`: The `<prefix>.eigenvec` file created in the [`get_eQTL_covariates`](../../eQTL_mapping/get_eQTL_covariates/) step
2. `peerX`: The `X.csv` file created by the `02_get_PEER_covariates.sh` script
3. `addCovFile`: A tab-separated file with additional covariates to be included in the final covariate file. There should be one line per sample, one column per covariate. The first column should be the sample IDs. For the initial MAGE publication, we included only sample sex as an additional covariate (i.e. in addition to genotype PCs and PEER factors)
4. `outPrefix`: The prefix of the formatted output file, to be used for sQTL mapping with FastQTL
<br><br>
