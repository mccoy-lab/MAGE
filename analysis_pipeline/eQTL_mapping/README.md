# MAGE *cis*-eQTL Mapping Pipeline

This directory contains all the componenets of the *cis*-eQTL mapping pipeline used as part of the initial MAGE publication. This pipeline is broadly split into 5 sections:
1. Prepare covariates (sex, genotype PCs, PEER factors) for QTL mapping
2. Run nominal and permutation passes with [FastQTL](https://github.com/francois-a/fastqtl)
3. Perform fine-mapping with [SuSiE](https://github.com/stephenslab/susieR)
4. Calculate credible set lead variant effect sizes with [aFC-n](https://github.com/PejLab/aFCn) (though note that we use a slightly modified version, found [here](https://github.com/dtaylo95/aFCn/))
5. Run genotype-continental group interaction test
<br><br>

## Prepare eQTL-mapping covariates

Code to calculate and prepare covariates for eQTL mapping is in the [`get_eQTL_covariates`](get_eQTL_covariates/) directory.

For *cis*-eQTL mapping, we included sample sex, the top five genotype PCs, and 60 [PEER factors](https://doi.org/10.1038%2Fnprot.2011.457) as covariates. The number of genotype PCs was chosen based on the proportion of genetic variance explained and their correlation with population labels. The number of PEER factors was chosen based on optimizations previously performed by the GTEx consortium.<br><br>

## Nominal and Permutation pass with FastQTL

Code to run nominal and permutation *cis*-eQTL mapping passes with [FastQTL](https://github.com/francois-a/fastqtl) is in the [`run_FastQTL`](run_FastQTL/) directory.

*cis*-eQTL mapping was performed for all autosomal genes that passed filtering thresholds, using inverse-normal transformed TMM values. Significant eGene-eVariant assocations were identified at a 5% FDR.<br><br>

## Fine-mapping with SuSiE

Code to perform *cis*-eQTL fine-mapping with [SuSiE](https://github.com/stephenslab/susieR) is in the ['run_susie'](run_susie/) directory.

Fine-mapping was performed for all eGenes identified in the FastQTL permutation pass using summary statistics from the FastQTL nominal pass. Up to 10 credible sets were identified per gene, at a minimum coverage of 95%.<br><br>

## Fine-mapped eQTL effect size estimation

Code to estimate eQTL effect sizes with [aFC-n](https://github.com/PejLab/aFCn) is in the [`run_aFCn`](run_aFCn/) directory.

We calculated the effect size of fine-mapped eQTL signals as log2 allelic fold change (log2(aFC); read more [here](https://doi.org/10.1101/gr.216747.116)). We used a slightly modified version of the [aFC-n](https://github.com/PejLab/aFCn) tool (our fork is [here](https://github.com/dtaylo95/aFCn/)) to calculate effect sizes for the lead variant of each SuSiE credible set. This tool allows for multiple causal signals per gene so that the effect size of any given lead variant is calculated conditional on all other SuSiE credible sets for the same gene.<br><br>

## Genotype-continental group interaction test

Code to perform the genotype-continental group interaction tests (single and multiple causal variant models) is in the [`run_interaction_test`](run_interaction_test/) directory.

In the first model, we test for an interaction between genotype and continental group for each lead variant from SuSiE (i.e. the same set of variants for which log2(aFC) was calculated). We require that each variant has MAF > 0.05 in at least two continental groups to test. The second model mirrors the first, except that we condition on all other SuSiE credible sets for the same gene, controlling for the additive effects of multiple causal variants.<br><br>
