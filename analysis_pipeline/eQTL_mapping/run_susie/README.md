# eQTL fine-mapping with SuSiE

This directory contains code to replicate eQTL fine-mapping analyses done for the initial MAGE publication. We used [SuSiE](https://github.com/stephenslab/susieR)
to perform fine-mapping, specifically using the `susie_rss` function to do fine-mapping using nominal results from FastQTL. The procedue compries two steps:
1. Run SuSiE to determine independent sets of putatively causal variants (credible sets)
2. Summarize and collate raw results from SuSiE

*Note: These steps should be run in order*<br><br>

## Run SuSiE using nominal results from FastQTL

The `01_run_susie.sh` script will perform all of the steps necessary to run SuSiE for a single gene. You can run this script in parallel across genes to speed up this step. For the initial MAGE publication, we ran SuSiE for all eGenes from FastQTL (at a 5% FDR).

The script takes six arguments:
1. `geneID`: The ID of the gene to run SuSiE with. This should be the geneID as it appears in the nominal results file and expression bed file
2. `inVCF`: The VCF file with variants to run eQTL mapping on. This should be the `autosomal.MAGE.v0.1.vcf.gz` file created in the [`subset_VCFs`](../../data_preparation/subset_VCFs/) step
3. `invNormTMMBed`: A BED file containing inverse normal transformed TMM values. This is the `<prefix>.filtered.inverse-normal-tmm.bed.gz` file created in the [`expression_quantification`](../../data_preparation/expression_quantification/) step, and previously used as input for FastQTL
4. `covFile`: This is the `.tab.gz` file created in the [`get_eQTL_covariates`](../get_eQTL_covariates/) step, and previously used as input for FastQTL
5. `nomResultsFile`: The output of running the FastQTL nominal pass
6. `outDir`: The directory to write results to. Two files will be written per gene: `<geneID>_full_pip_results.txt` and `<geneID>_cs_results.txt`
<br><br>

## Collect and format single-gene SuSiE results

The `02_collate_results.sh` script will collect the single-gene results created by the `01_run_susie.sh` script and determine lead eQTLs for each SuSiE credible set.

The script takes four arguments:
1. `geneOutFileList`: A three-column tab-separated file with geneID, path to "cs_results.txt" file, path to "_full_pip_results.txt". The file should have header fields: `geneID`, `csResultsFile`, and `fullResultsFile`
2. `outCSSummary`: File to write credible set-level summary of SuSiE results
3. `outFullSNPSummary`: File to write full SNP-level results to
4. `outSigSNPSummary`: A subset of `outFullSNPSummary` containing only the eQTLs that are in a SuSiE credible set
5. `outLeadHitSummary`: A subset of `outLeadGitSummary` containing only a single eQTL with the highest PIP from each credible set.

**Note: The `outLeadHitSummary` file  is the file used for most downstream analyses**

After running the `02_collate_results.sh` script, we reccomend either compressing or deleting the `outDir` directory from `01_run_susie.sh`, as it is likely to be quite large, and better represented by the summary files created by the `02_collate_results.sh` script.<br><br>
