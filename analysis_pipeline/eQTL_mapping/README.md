# MAGE *cis*-eQTL Mapping Pipeline

This directory contains all the componenets of the *cis*-eQTL mapping pipeline used as part of the initial MAGE publication. This pipeline is broadly split into X sections:
1. Prepare covariates (sex, genotype PCs, PEER factors) for QTL mapping
2. Run nominal and permutation passes with [FastQTL](https://github.com/francois-a/fastqtl)
3. Perform fine-mapping with [SuSiE](https://github.com/stephenslab/susieR)
4. Calculate lead variant effect sizes with [aFC-n](https://github.com/PejLab/aFCn) (though note that we use a slightly modified version, found [here](https://github.com/dtaylo95/aFCn/))
5. Run genotype-continental group interaction test
<br><br>

## Prepare eQTL-mapping covariates

Code to calculate and prepare covariates for eQTL mapping is in the [`get_eQTL_covariates`](get_eQTL_covariates/) directory.

For *cis*-eQTL mapping, we included sample sex, the top five genotype PCs, and 60 [PEER factors](https://doi.org/10.1038%2Fnprot.2011.457) as covariates. The number of genotype PCs was chosen based on the proportion of genetic variance explained and their correlation with population labels. The number of PEER factors was chosen based on optimizations previously performed by the GTEx consortium.<br><br>

## Nominal and Permutation pass with FastQTL

Code to run nominal and permutation *cis*-eQTL mapping passes with [FastQTL](https://github.com/francois-a/fastqtl) is in the [`run_FastQTL`](run_FastQTL/) directory.

*cis*-eQTL mapping was performed for all autosomal genes that passed filtering thresholds, using inverse-normal transformed TMM values. Significant eGene-eVariant assocations were identified at a 5% FDR.<br><br>

## Fine-mapping with SuSiE

