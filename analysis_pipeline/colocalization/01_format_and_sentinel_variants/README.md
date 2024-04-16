# Formatting of GWAS Catalog sumstats and selection of sentinel variants

This directory contains code to re-format harmonized GWAS summary statistics downloaded from the [GWAS Catalog](https://www.ebi.ac.uk/gwas), and to select sentinel variants based on these summary statistics that will serve as windows of colocalization in the next steps. Additionally, using these sentinel variants, we construct LD matrices from genotypes for samples in the GWAS. As such, there are three broad steps:
1. Re-format GWAS catalog summary statistics
2. Select sentinel variants from re-formatted summary statistics
3. Construct LD matrices from GWAS sample genotypes around each sentinel variant
<br><br>

## Formatting GWAS Catalog sumstats

The `01_format_sumstats.sh` script will re-format harmonized GWAS catalog summary statistics into a structure that is more useful for downstream analyses.

The script takes two arguments:
1. `inFile`: The input sumstats file downloaded from the GWAS Catalog. This filename should end with ".h.tsv.gz". We used summary statistics reported on the GRCh38 reference
2. `outFile`: The file to write re-formatted summary statistics to

This script assumes the input file contains summary statistics for just a single trait. You should run it once for each trait sumstats file.<br><br>

## Selecting sentinel variants

The `02_get_sentinel_variants.sh` script will select sentinel variants from the re-formatted sumstats file produced in the previous step. Sentinel variants are selected by iteratively choosing the variant with the most signficant p-value, and then removing from selection all variants within some distance of the selected variant.

The script takes two arguments:
1.  `formattedSumstatsFile`: The formatted summary statistics from the previous step
2.  `outSentinelVariantsFile`: An output file listing each of the selected sentinel variants

The script asssumes a colocalization window of +/-500kbp around each variant. This means sentinel variants will be selected that are at least twice that distance away from eachother (i.e. 1Mbp). Additionally, the script only selects sentinel variants with `log(p-value) >= 5e-8`. Feel free to change these parameters if you use this script for your own work.<br><br>

## Calculating LD matrices around sentinel variants

The `03_get_sentVar_LD_matrix.sh` script will construct LD matrices around each sentinel variant from sample genotypes. This is done so that you can use a matched LD matrix when running fine-mapping on the GWAS results. For our analysis, we used GWAS summary statistics from [PAGE](https://www.nature.com/articles/s41586-019-1310-4), and we constructed LD matrices from overlapping samples in the [TOPMed imputation panel](https://www.nature.com/articles/s41586-021-03205-y). For each trait, we identified the PAGE samples used in GWAS of that trait, and subset the TOPMed genotypes to samples that overlap between the two datasets. We note that these genotypes are protected access.

This script is run for a single sentinel variant takes five arguments:
1. `sentinelVariant`: The sentinel variant to construct an LD matrix around. Should be formatted chrom:pos:ref:alt (e.g. chr19:16092494:C:A). Position is expected to be 1-based
2. `sumstatsFile`: The formatted summary statistics from the first step
3. `inVCF`: The VCF to pull genotypes from
4. `selectedSampleList`: A file with the samples to subset from the input VCF, one per line
5. `outPrefix`: The prefix (including path) of the output files

### Output

This script will produce two files: a `.ld` file with the LD matrix, and a `.snps` file listing the corresponding SNPs in the LD matrix.

We recommend parallelizing across sentinel variants if possible when running this script.<br><br>
