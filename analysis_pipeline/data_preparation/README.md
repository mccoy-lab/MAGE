# MAGE Data Preparation

This directory contains all the components of the data preparation pipeline used as part of the initial MAGE publication. This pipeline is broadly split into three sections:
1. Preparing MAGE VCFs from 1KGP VCFs
2. Quantifying expression from raw reads using Salmon
3. Alignment with STAR and splicing quantification with Leafcutter

*Note: all analyses were performed on the GRCh38 reference*<br><br>

## Preparing MAGE VCFs

Code to produce MAGE VCFs is in the [`subset_VCFs`](subset_VCFs/) directory.

The samples included in MAGE are a subset of the [1000 Genomes Project (1KGP)](https://doi.org/10.1038/nature15393). The [New York Genome Center (NYGC)](https://doi.org/10.1016/j.cell.2022.08.004) previously generated a set of high-quality variant calls for the full 1KGP data set using high-coverage whole genome sequencing, and we generate a MAGE VCF as a subset of the NYGC VCF for these samples (downloaded from the [1KGP FTP site](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/).<br><br>

## Expression Quantification

Code to perform expression quantification is in the [`expression_quantification`](expression_quantification/) directory.

Expression quantification was performed with raw reads using [Salmon](https://combine-lab.github.io/salmon/) with [GENCODE v38](https://www.gencodegenes.org/human/release_38.html) transcript annotations.

Transcript-level pseudocounts from Salmon were merged to gene-level estimates, filtered, and normalized (normalization depends on downstream application).<br><br>

## Splicing Quantification

Code to perform splicing quantification is in the [`splicing_quantification`](splicing_quantifcation/) directory.

Reads were first aligned to the GRCh38 reference using the [STAR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/) aligner. These alignments were then used as input with [Leafcutter](https://www.nature.com/articles/s41588-017-0004-9) to quantify splicing in an annotation agnostic manner. Splicing quantifications were then filtered. For sQTL discovery, filtered quants were mapped to genes and normalized.<br><br>
