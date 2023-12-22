# sQTL fine-mapping with SuSiE

This directory contains code to replicate sQTL fine-mapping analyses done for the initial MAGE publication. We used [SuSiE](https://github.com/stephenslab/susieR)
to perform fine-mapping, specifically using the `susie_rss` function to do fine-mapping using nominal results from FastQTL. The procedure compries two steps:
1. Run SuSiE to determine independent sets of putatively causal variants (credible sets)
2. Summarize and collate raw results from SuSiE

*Note: These steps should be run in order*<br><br>

## Index nominal results for fast look-up

The nominal results output from FastQTL can be quite large. As we are using this file for fine-mapping, it is beneficial to index the nominal results file with tabix to speed up look up of records in this file.

The `00_index_nominal_results.sh` file will index the nominal results file from FastQTL for use with SuSiE.

The script takes three arguments:
1. `nominalResultsFile`: The output of running the FastQTL nominal pass
2. `normedRatiosBed`: A BED file containing normalized intron excision ratios used as input for sQTL mapping. This is the `<prefix>_perind.counts.normed.genes.filtered.bed.gz` file created in the [`splicing_quantification`](../../data_preparation/splicing_quantification/) step, and previously used as input for FastQTL
3. `outSortedFile`: File to write tabix-indexed nominal results to
<br><br>

## Run SuSiE using nominal results from FastQTL

The `01_run_susie.sh` script will perform all of the steps necessary to run SuSiE for all introns of a single gene. You can run this script in parallel across genes to speed up this step. For the initial MAGE publication, we ran SuSiE for all sGenes from FastQTL (at a 5% FDR).

The script takes seven arguments:
1. `geneID`: The ID of the gene to run SuSiE with. This should be the geneID as it appears in the nominal results file and expression bed file
2. `inVCF`: The VCF file with variants to run sQTL mapping on. This should be the `autosomal.MAGE.v0.1.vcf.gz` file created in the [`subset_VCFs`](../../data_preparation/subset_VCFs/) step
3. `normedRatiosBed`: A BED file containing normalized intron excision ratios used as input for sQTL mapping. This is the `<prefix>_perind.counts.normed.genes.filtered.bed.gz` file created in the [`splicing_quantification`](../../data_preparation/splicing_quantification/) step, and previously used as input for FastQTL
4. `covFile`: This is the `.tab.gz` file created in the [`get_sQTL_covariates`](../get_esTL_covariates/) step, and previously used as input for FastQTL
5. `phenGroupsFile`: A file mapping introns to genes. This is the `<prefix>_perind.counts.normed.genes.filtered.phenotype_groups.txt` file created in the [`splicing_quantification`](../../data_preparation/splicing_quantification/) step, and previously used as input for FastQTL
6. `nomResultsFile`: The output of running the FastQTL nominal pass
7. `outDir`: The directory to write results to. Two files will be written per intron: `<intronID>_full_pip_results.txt` and `<intronID>_cs_results.txt`
<br><br>

## Collect and format single-gene SuSiE results

The `02_collate_results.sh` script will collect the single-intron results created by the `01_run_susie.sh` script, merge intron-level credible sets into gene-level credible sets, and determine lead sQTLs for each gene-level credible set.

The script takes six arguments:
1. `intronOutFileList`: A four-column tab-separated file with intronID, geneID, path to "cs_results.txt" file, path to "_full_pip_results.txt". The file should have header fields: `intronID`, `geneID`, `csResultsFile`, and `fullResultsFile`
2. `outCSSummary`: File to write credible set-level summary of SuSiE results
3. `outFullSNPSummary`: File to write full SNP-level results to
4. `outSigSNPSummary`: A subset of `outFullSNPSummary` containing only the sQTLs that are in a SuSiE credible set
5. `outGeneSigSNPSummary`: File to write gene-level credible set fine-mapped sQTLs to
6. `outLeadHitSummary`: A subset of `outGeneSigSNPSummary` containing only a single sQTL from each gene-level credible set. This is the sQTL with the highest PIP in the intron-level credible set (of those making up the gene-level credible set) with the highest coverage

**Note: The `outLeadHitSummary` file  is the file used for most downstream analyses**

After running the `02_collate_results.sh` script, we reccomend either compressing or deleting the `outDir` directory from `01_run_susie.sh`, as it is likely to be quite large, and better represented by the summary files created by the `02_collate_results.sh` script.<br><br>
