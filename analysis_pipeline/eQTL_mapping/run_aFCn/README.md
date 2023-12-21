# eQTL effect size estimation

This directory contains code to replicate eQTL effect size estimation done for the initial MAGE publication. We used [aFC-n](https://github.com/PejLab/aFCn) (though note that we use a slightly modified version, found [here](https://github.com/dtaylo95/aFCn/)) to calculate the effect of each fine-mapped lead eQTL from SuSiE as log2 allelic fold change (log2(aFC)). You can read more about the interpretation of log2(aFC) [here](https://doi.org/10.1101/gr.216747.116). The aFC-n tool allows for multiple causal signals per gene such that the effect size of any given lead variant is calculated conditional on all other SuSiE credible sets for the same gene. We note that our fork of the aFC-n tool does not modify the underlying model at all, it simply relaxes some of the file formatting requirements of the aFC-n tool.<br><br>

# log2(aFC) calculation with aFC-n

The `01_run_aFCn.sh` script formats inputs for and runs the aFC-n tool to calculate log2(aFC) effect sizes for each lead variant from SuSiE.

The script takes seven arguments:
1. `leadHitSummaryFile`: This should be the lead hit summary output file from the ['run_susie`](../run_susie/) step
2. `variantInfoFile`: Tab-separated file with position info for SuSiE lead eQTLs. The file must have header row. First column must be the variantID as it appears in the SuSiE output. It must also contain at least a "CHROM" column and a (1-based) "POS" column
3. `inVCF`: The same VCF file used for eQTL mapping and fine-mapping. This should be the `autosomal.MAGE.v0.1.vcf.gz` file created in the [`subset_VCFs`](../../data_preparation/subset_VCFs/) step
4. `expressionBed`: A BED file containing log2 DESeq2-normalized pseudocounts. This is the `<prefix>.filtered.DESeq2-normalized-logged-pseudocounts.bed.gz` file created in the [`expression_quantification`](../../data_preparation/expression_quantification/) step
5. `covFile`: The same covariate file used eQTL mapping and fine-mapping. This is the `.tab.gz` file created in the [`get_eQTL_covariates`](../get_eQTL_covariates/) step
6. `out_eQTLs`: This is a four-column tab-delimited file describing the eQTLs to test with aFC-n. It is used as input for both aFC-n and the [interaction test](../run_interaction_test/)
7. `out_aFCnFile`: The file to write output log2(aFC) results to

### Software requirements

This script assumes you have the `afcn.py` script (installed as part of [aFC-n](https://github.com/dtaylo95/aFCn)) on your PATH<br><br>
