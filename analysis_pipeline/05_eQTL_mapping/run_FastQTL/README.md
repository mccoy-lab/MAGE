# *cis*-eQTL mapping nominal/permutation pass with FastQTL

This directory contains code to replicate the nominal eQTL mapping analyses done for the initial MAGE publication. We used [FastQTL](https://github.com/francois-a/fastqtl) to discover significant nominal eQTLs, and the procedure comprises three broad steps:
1. A nominal pass to map associations between each gene and surrounding variants
2. A permutation pass to identify significant eGenes (genes with at least one eQTL) and determine significance thresholds for nominal associations
3. A step to subset nominal associations to only those that are statistically significant

*Note: steps 1 and 2 can be run in any order (or in parallel), step 3 should be run after both previous steps are completed*<br><br>

## FastQTL nominal pass

The `01_run_fastqtl_nominal.sh` script will run [FastQTL](https://github.com/francois-a/fastqtl) in its default "nominal" mode. For each gene, FastQTL will regress expression of that gene onto the genotype of each variant near that gene. For the intial MAGE publication, we filter to variants with MAF > 0.01 within 1Mbp up- or down-stream of the TSS of the gene.

The script takes five arguments:
1. `inVCF`: The VCF file with variants to run eQTL mapping on. This can be any of the `.MAGE.v0.1.vcf.gz` files created in the [`subset_VCFs`](../../data_preparation/subset_VCFs/) step. For the initial MAGE publication, we ran each chromosome separately. You can run all chromosomes together, but you will need to change the `chunks` argument when running FastQTL
2. `invNormTMMBed`: A BED file containing inverse normal transformed TMM values. This is the `<prefix>.filtered.inverse-normal-tmm.bed.gz` file created in the [`expression_quantification`](../../data_preparation/expression_quantification/) step. For the initial MAGE publication, we ran each chromosome separately, subsetting this bed file to a single chromosome, which we then used as input here. 
3. `covFile`: This is the `.tab.gz` file created in the [`get_eQTL_covariates`](../get_eQTL_covariates/) step.
4. `outDir`: Directory to write output nominal results files to
5. `outPrefix`: Prefix of output nominal results files

### Software requirements

This script assumes you have the `run_FastQTL_threaded.py` script (installed as part of [FastQTL](https://github.com/francois-a/fastqtl/blob/master/INSTALL)) on your PATH.<br><br>

## FastQTL permutation pass

The `02_run_fastqtl_permuation.sh` script will run [FastQTL](https://github.com/francois-a/fastqtl) in its "permutation" mode which allows the user to 1) determine signficant eGenes at some FDR threshold, and 2) set a significance threshold for nominal assocations at the same FDR.

The script takes the same five arguments as the `01_run_fastqtl_nominal.sh` script above. As before, we ran the permutation pass on each chromosome separately, but users can choose to run all chromosomes together, making sure to modify the `chunks` argument in the script. 

### Software requirements

As above, this script assumes you have the `run_FastQTL_threaded.py` script (installed as part of [FastQTL](https://github.com/francois-a/fastqtl/blob/master/INSTALL)) on your PATH.<br><br>

## Finding significant eGenes/eQTLs

After running the previous two steps, the `03_get_sig_results.sh` script will determine significant eGenes, as well as significant eQTL associations. For the initial MAGE publication, we used a 5% FDR significance threshold to determine eGenes/eQTLs.

The script takes four arguments:
1. `nominalResultsFile`: The output of running the FastQTL nominal pass. If you ran each chromosome separately (as we did), you should merge the results into a single file to be used as input here
2. `permutationResultsFile`: The output of running the FastQTL permutation pass. If you ran each chromosome separately (as we did), you should merge the results into a single file to be used as input here
3. `out_eGeneListFile`: File to write the list of eGenes (at a 5% FDR) to. This file will be used when running SuSiE
4. `out_sigNominalResultsFile`: File to write significant eQTL associations (at a 5% FDR) to. This will be a subset of the `nominalResultsFile` file
<br><br>
