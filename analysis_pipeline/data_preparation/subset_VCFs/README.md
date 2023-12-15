# MAGE VCF Preparation

This directory contains code to prepare the MAGE VCF used for QTL mapping. Because the samples included in MAGE are a subset of samples in the [1000 Genomes Project](https://doi.org/10.1038/nature15393), we can generate a MAGE VCF by subsetting the 1KGP VCF produced by the [New York Genome Center (NYGC)](https://doi.org/10.1016/j.cell.2022.08.004).

## Subset NYGC 1KGP VCF to MAGE samples

The `01_subset_VCFs.sh` script will subset a 1KGP VCF to just the samples included in MAGE. We used the phased VCF from the NYGC (download [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/)).

This script takes three arguments:
1. `in_vcfDir`: Directory containing chromosome-split VCFs from 1KGP
2. `sampleListFile`: File with MAGE samples to subset 1KGP VCFs to (one per line)
3. `out_vcfDir`: Directory to write subset VCFs to
<br><br>
