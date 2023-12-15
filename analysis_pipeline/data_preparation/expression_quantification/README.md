# MAGE Expression Quantification

This directory contains code to replicate the expression quantification procedure does as part of the initial MAGE publication. This procedure is broadly split into three sections:
1. Initial transcript-level quantification with [Salmon](https://combine-lab.github.io/salmon/)
2. Filtering and normalization for eQTL mapping analysis
3. Filtering and normalization for Differential Expression analyses

*Note: these analyses steps should be run in order*<br><br>

## Transcript-level quantification with Salmon

The `01_run_salmon.sh` script will perform expression quantification with [Salmon](https://combine-lab.github.io/salmon/), an alignment free expression quantification tool. For the intial MAGE publication, we used the [GENCODE v38]() transcript sequences file for expression quantification (download [here](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz))

The script takes four arguments:
1. `fastqListFile`: A tab-separated file with paths to MAGE FASTQ files. This file should have three columns: 1) libraryID, 2) path to pair1 FASTQ file, 3) path to pair2 FASTQ file
2. `gencodeTranscriptsFasta`: The gencode transcript sequences FASTA file (we used the Gencode v38 file [here](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz)
3. `outIndexDir`: Directory to write salmon index files to
4. `outSalmonDir`: Directory to write salmon quantifications to

### Software requirements
No additional software needed.<br><br>

## Filtering and normalization for eQTL analyses

The `02_filter_normalize.sh` script will merge transcript-level quantifications to gene-level quantifications and filter out lowly expressed genes (only genes with counts >= 6 and TPM >= 0.1 i at least 20% of samples are kept). It will also output four different quantifications:
1. Raw counts (no normalization)
2. TMM values
3. Inverse-normalized TMM values (used for eQTL mapping)
4. DESeq2-normalized logged counts (used for calculation of allelic fold change [aFC])

The script takes four arguments:
1. `transcriptGTF`: The Gencode transcript GTF file, used for mapping transcripts to genes (we used the Gencide v38 files [here](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz))
2. `salmonQuantListFile`: A tab-separated file with paths to Salmon `quant.sf` files from `01_run_salmon.sh`. This files should have two columns: 1) sampleID (to be used in output files) and 2) path to corresponding `quant.sf` file. For eQTL analyses, for each of the 24 samples sequenced in triplicate, we selected only the replicate with the most reads for further analysis.
3. `outQuantsDir`: Directory to write gene-level filtered expression quantification files to
4. `outPrefix`: Prefix of output files in `outQuantsDir`

### Software requirements
No additional software needed.<br><br>
