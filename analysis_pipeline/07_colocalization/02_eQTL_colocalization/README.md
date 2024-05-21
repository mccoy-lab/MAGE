# MAGE eQTL colocalization

This directory contains code to perform colocalization between MAGE cis-eQTLs and complex trait GWAS. We used [Coloc](https://chr1swallace.github.io/coloc/) tool with SuSiE (as implemented in the `coloc` package) to perform colocalization between cis-eQTLs in MAGE and 25 complex trait GWAS from [the PAGE study](https://www.nature.com/articles/s41586-019-1310-4). The code in this directory is organized into two steps:
1. Colocalization (one sentinel variant at a time)
2. Collation of sentinel variant-level results
<br><br>

## Colocalization between eQTLs and complex trait GWAS

The `01_run_colocalization.sh` script performs colocalization between MAGE eQTLs and GWAS results in a window around a chosen GWAS sentinel variant. At a high level, the script first generates a window around the input sentinel variant, identifies eGenes whose TSS falls in that window, subsets eQTL nominal pass summary statistics for those genes to that window, and performs colocalization between the GWAS results and eQTL results for each of those eGenes. Notably colocalization is performed using credible sets from SuSiE, and colocalization is performed between each pair of credible sets between the eQTL and GWAS results.

We report "moderate" colocalizations as those with a posterior colocalization probability (`PP.H4`) of at least 0.5, and strong colocalizations as those with a `PP.h4` of at least 0.8.

The script akes 10 arguments:
1. `sentinelVariant`: The sentinel variant to construct an LD matrix around. Should be formatted chrom:pos:ref:alt (e.g. chr19:16092494:C:A). Position is expected to be 1-based
2. `sumstatsFile`: Formatted summary statistics for the corresponding trait
3. `ldPrefix`: Prefix (including path) of GWAS LD matrix (and snp list) for this sentinel variant (produced by the [`03_get_sentVar_LD_matrix.sh`](../01_format_and_sentinel_variants/03_get_sentVar_LD_matrix.sh) script
4. `numSamples`: Number of samples in the GWAS
5. `traitType`: GWAS trait type. "quant" for quantitative traits, "cc" for case/control traits
6. `eQTL_VCF`: VCF file used for eQTL mapping
7. `eQTL_permutationResultsFile`: Tabixed and formatted eQTL results from FastQTL permutation pass. This script assumes these files are formatted as they are in the publicly available[MAGEv0.1 Zenodo](https://zenodo.org/doi/10.5281/zenodo.10535719)
8. `eQTL_nominalResultsFile`: Tabixed and formatted eQTL results from FastQTL nominal pass. As before, this script assumes these files are formatted as they are in the  [MAGEv0.1 Zenodo](https://zenodo.org/doi/10.5281/zenodo.10535719)
9. `eQTL_covFile`: Covariates file used for eQTL mapping (again, formatted as in the  [MAGEv0.1 Zenodo](https://zenodo.org/doi/10.5281/zenodo.10535719))
10. `outPrefix`: Prefix (including path) of output files

### Output

This script will produce four files:
1. a `.coloc.summary.txt` file describing the number of credible sets and moderate/strong colocalizations for each gene
2. a `.coloc.fullResults.txt` file describing the results of each tested colocalizations, including posterior probabilities
3. a `.coloc.GWAS_CS.txt` file describing all SuSiE credible sets from the GWAS, AS LONG AS there was at least one moderate colocalization
4. a `.coloc.eQTL_CS.txt` file describing all SuSiE credible sets for each eGene with at least one moderate colocalization

We recommend parallelizing across sentinel variants if possible when running this script.<br><br>

## Collating colocalization results

The `02_collate_coloc_results.sh` script collates the sentinel variant-level results file into a single set of files across all tested sentinel variants.

The script takes two arguments:
1. `mapFile`: A four-column tab-separated file (no header) with trait accession, trait description, sentinel variant id, and coloc results prefix (including path) for each sentinel variant
2. `outPrefix`: Prefix (including path) of output files

### Output

The script will output four files:
1. Collated results from all `.coloc.fullResults.txt` files
2. The subset of the above file with moderate (`PP.H4 >= 0.5`) colocalizations
3. Collated results from all `.coloc.GWAS_CS.txt` files
4. Collated results from all `.coloc.eQTL_CS.txt` files
<br><br>
