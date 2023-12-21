# MAGE *cis*-sQTL Mapping Pipeline

This directory contains all the componenets of the *cis*-sQTL mapping pipeline used as part of the initial MAGE publication. This pipeline is broadly split into three sections:
1. Prepare covariates (sex, genotype PCs, PEER factors) for QTL mapping
2. Run nominal and permutation passes with [FastQTL](https://github.com/francois-a/fastqtl)
3. Perform fine-mapping with [SuSiE](https://github.com/stephenslab/susieR)
<br><br>

## Prepare sQTL mapping covariates

Code to calculate and prepare covariates for sQTL mapping is in the [`get_sQTL_covariates`](get_sQTL_covariates/) directory.

For *cis*-eQTL mapping, we included sample sex, the top five genotype PCs, and 15 [PEER factors](https://doi.org/10.1038%2Fnprot.2011.457) as covariates. The number of genotype PCs was chosen based on the proportion of genetic variance explained and their correlation with population labels. The number of PEER factors was chosen based on optimizations previously performed by the GTEx consortium.<br><br>

## Nominal and Permutation pass with FastQTL

Code to run nominal and permutation *cis*-sQTL mapping passes with [FastQTL](https://github.com/francois-a/fastqtl) is in the [`run_FastQTL`](run_FastQTL/) directory.

*cis*-sQTL mapping was performed for all autosomal introns that passed filtering thresholds, using normalized intron exision ratios. Significant eGene-eVariant assocations were identified at a 5% FDR.<br><br>

## Fine-mapping with SuSiE

Code to perform *cis*-sQTL fine-mapping with [SuSiE](https://github.com/stephenslab/susieR) is in the ['run_susie'](run_susie/) directory.

Fine-mapping was performed for all sGenes identified in the FastQTL permutation pass using summary statistics from the FastQTL nominal pass. Up to 10 credible sets were identified per intron, at a minimum coverage of 95%.

Intron-level credible sets were then iteratively merged into gene-level credible sets.<br><br>
