# ADMIXTURE and PCA comparison of datasets

This directory contains code to replicate the PCA and ADMIXTURE analysis performed in the initial MAGE publication. Specifically, the code loops over all autosomes and extracts a set of LD-pruned variants that exceed a specified minor allele frequency threshold.

These variants are then extracted from each of the GTEx and AFGR (MKK population, as these samples were not part of the 1000 Genomes Project) sample VCFs.

Finally, the three datasets are combined and PCA and ADMIXTURE are applied, varying k over a range from 6 to 10. 

The script takes four arguments:

1. `kgpVCFdir`: The directory containing 1000 Genomes VCF files (one per chromosome) from the New York Genome Center 8/5/2020 release with names such as `CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz`
2. `gtexVCF`: The path to the GTEx VCF file
3. `afgrVCF`: The path to the AFGR MKK VCF file
4. `outDir`: The directory where output files will be written

Note that the code assumes that PLINK and ADMIXTURE are in the user's PATH.
