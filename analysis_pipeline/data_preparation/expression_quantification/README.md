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
The script assumes you have `salmon` installed and on your PATH.<br><br>

## Filtering and normalization for eQTL analyses

