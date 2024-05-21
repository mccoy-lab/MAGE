# MAGE VCF Preparation

This directory contains code to prepare the MAGE VCF used for QTL mapping. Because the samples included in MAGE are a subset of samples in the [1000 Genomes Project](https://doi.org/10.1038/nature15393), we can generate a MAGE VCF by subsetting the 1KGP VCF produced by the [New York Genome Center (NYGC)](https://doi.org/10.1016/j.cell.2022.08.004).

## Subset NYGC 1KGP VCF to MAGE samples

The `01_subset_VCFs.sh` script will subset a 1KGP VCF to just the samples included in MAGE. We used the phased VCF from the NYGC (download [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/)).

This script takes three arguments:
1. `in_vcfDir`: Directory containing chromosome-split VCFs from 1KGP
2. `sampleListFile`: File with MAGE samples to subset 1KGP VCFs to (one per line)
3. `out_vcfDir`: Directory to write subset VCFs to
<br><br>

## Subset NYGC 1KGP chrX VCF and transform genotypes

The `02_prep_chrX_VCF.sh` script will subset the 1KGP chrX VCF to just the samples included in MAGE. We used the phased VCF from the NYGC (download [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/)).

In addition to subsetting the VCF to a pre-selected set of samples, the script will 1) assign variant IDs, because these are missing from the version of the 1KGP NYGC VCF we used and 2) convert haploid genotypes in XY samples to diploid genotypes. The second step is necessary to include XX and XY samples together when performing QTL mapping. There are some assumptions that this approach makes:
1. The chrX pseudo-autosomal regions (PARs) represent both the chrX PARs as well as the chrY PARs. As such, all samples (XX and XY) are diploid for all variants within the PARs
2. All genes within the PARs escape X-inactivation
3. All genes outside the PARs are subject to X-inactivation in XX samples

This script takes three arguments:
1. `inVCF`: The chrX VCF from 1KGP
2. `sampleListFile`: File with MAGE samples to subset 1KGP VCFs to (one per line)
3. `outVCFgz`: The output subset, genotype-transformed VCF (should end with .gz)
<br><br>
