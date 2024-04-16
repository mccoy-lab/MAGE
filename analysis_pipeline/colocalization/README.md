# MAGE colocalization pipeline

This directory contains code to replicate the colocalization analysis done for the initial MAGE publication. We used [Coloc](https://chr1swallace.github.io/coloc/) tool with SuSiE (as implemented in the `coloc` package) to perform colocalization between cis-e/sQTLs in MAGE and 25 complex trait GWAS from [the PAGE study](https://www.nature.com/articles/s41586-019-1310-4). This pipeline is broadly split into three sections:
1. Selection of "sentinel variant" from GWAS summary statistics
2. Colocalization between MAGE cis-eQTLs and PAGE GWAS
3. Colocalization between MAGE cis-sQTLs and PAGE GWAS

The first step must be run first, but the second and third steps do not need to be run sequentially.<br><br>

## Selecting sentinel variants

Harmonized GWAS summary statistics (on GRCh38) for PAGE traits were downloaded from the [GWAS Catalog](https://www.ebi.ac.uk/gwas/publications/31217584).

Code to format downloaded summary statistics and select sentinel variants is in the [`01_format_and_sentinel_variants/`](01_format_and_sentinel_variants/) directory.

Colocalization is inherently built on fine-mapping. Unforunately, for a GWAS, attempting to fine-map the entire genome is prohibitive. As such, we attempt to first select indepent GWAS signals (termed "sentinel variants"), and then run fine-mapping and colocalization in a window around each of those independent GWAS signals.

For the initial publication, we selected sentinel variants that are at least 1Mbp away from each other, and run colocalization in a +/-500kbp window around each sentinel variant.<br><br>

## eQTL colocalization

Code to perform colocalization between MAGE cis-eQTLs and PAGE GWAS is in the [`02_eQTL_colocalization/`](02_eQTL_colocalization/) directory.

For each GWAS sentinel variant, we perform colocalization with each nominal eGene whose TSS is within a +/- 500kbp window around the sentinel variant.<br><br>

## sQTL colocalization

Code to perform colocalization between MAGE cis-sQTLs and PAGE GWAS is in the [`03_sQTL_colocalization/`](03_sQTL_colocalization/) directory.

For each GWAS sentinel variant, we perform colocalization with all introns of each nominal sGene whose TSS is within a +/- 500kbp window around the sentinel variant.<br><br>
