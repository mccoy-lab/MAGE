# eQTL genotype-by-continental group interaction test

This directory contains code to replicate the eQTL genotype-by-continental group interaction tests done as part of the initial MAGE publication to discover heterogenous effect eQTLs (he-eQTLs). We test for significant genotype-by-continental group interactions for each lead eQTL from SuSiE (as long as it passes MAG thresholds). Broadly, we consider two models:
1. A model that, for each lead variant, tests for significant genotype-by-continental group interactions, conditional only on the eQTL mapping covariates
2. A second model that mirrors the first, but ALSO conditions on the additive effects of all other SuSiE credible sets for the same gene, thereby controlling for the additive effect of multiple causal variants.

For both models, we only test interactions for variants with MAF > 0.05 in at least two continental groups. For the second model, all other lead variants (for the same gene) are included as covariates, regardless of MAF.<br><br>

## Running the interaction tests

The `01_run_contGroup_interaction.sh` script runs both of the models described above.

The script takes seven arguments:
1. `eQTLFile`: Four-column tab-delimted file of eQTLs to test. This is the `out_eQTLs` file created in the [`run_aFCn`](../run_aFCn/) step
2. `inVCF`: The same VCF file used for eQTL mapping and fine-mapping. This should be the `autosomal.MAGE.v0.1.vcf.gz` file created in the [`subset_VCFs`](../../data_preparation/subset_VCFs/) step
3. `expressionBed`: A BED file containing inverse normal transformed TMM values. This is the `<prefix>.filtered.inverse-normal-tmm.bed.gz` file created in the [`expression_quantification`](../../data_preparation/expression_quantification/) step, and used for eQTL mapping and fine-mapping
4. `populationTable`: A two-column tab-separated files with sampleIDs and the corresponding population labels. A header line is expected
5. `covFile`: The same covariate file used eQTL mapping and fine-mapping. This is the `.tab.gz` file created in the [`get_eQTL_covariates`](../get_eQTL_covariates/) step
6. `outFileSingle`: File to write single-variant output results to
7. `outFileMulti`: File to write multi-variant output results to
<br><br>
