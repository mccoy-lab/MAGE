# MAGE Expression Quantification

This directory contains code to replicate the expression quantification procedure done as part of the MAGE v1.1 release. This procedure is broadly split into three sections:
1. Initial transcript-level quantification with [Salmon](https://combine-lab.github.io/salmon/)
2. Filtering and normalization for eQTL mapping analysis
3. Differential expression analysis and generating expression matrices

*Note: these steps should be run in order*<br><br>

## Transcript-level quantification with Salmon

The `01_run_salmon.sh` script will perform expression quantification with [Salmon](https://combine-lab.github.io/salmon/), an alignment free expression quantification tool. For the intial MAGE publication, we used the [GENCODE v38](https://www.gencodegenes.org/human/release_38.html) transcript sequences file for expression quantification (download [here](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz))

The script takes four arguments:
1. `fastqListFile`: A tab-separated file with paths to MAGE FASTQ files. This file should have three columns: 1) sampleID, 2) path to pair1 FASTQ file, 3) path to pair2 FASTQ file
2. `gencodeTranscriptsFasta`: The gencode transcript sequences FASTA file (we used the Gencode v38 file [here](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz))
3. `outIndexDir`: Directory to write salmon index files to
4. `outSalmonDir`: Directory to write salmon quantifications to
<br><br>

## Filtering and normalization for eQTL analyses

The `02_filter_normalize.sh` script will merge transcript-level quantifications to gene-level quantifications and filter out lowly expressed genes (only genes with counts >= 6 and TPM >= 0.1 i at least 20% of samples are kept). It will also output four different quantifications:
1. Raw counts (no normalization)
2. TMM values
3. Inverse-normalized TMM values (used for eQTL mapping)
4. DESeq2-normalized logged counts (used for calculation of allelic fold change [aFC])

The script takes four arguments:
1. `transcriptGTF`: The Gencode transcript GTF file, used for mapping transcripts to genes (we used the Gencode v38 files [here](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz))
2. `salmonQuantListFile`: A tab-separated file with paths to Salmon `quant.sf` files from `01_run_salmon.sh`. This files should have two columns: 1) sampleID (to be used in output files) and 2) path to corresponding `quant.sf` file. For eQTL analyses, for each of the 24 samples sequenced in triplicate, we selected only the replicate with the most reads for further analysis.
3. `outQuantsDir`: Directory to write gene-level filtered expression quantification files to
4. `outPrefix`: Prefix of output files in `outQuantsDir`
<br><br>

## Differential expression analysis and generating expression matrices

The `03_deseq.ipynb` Jupyter notebook details the steps for generating expression matrices and performing differential expression analyses with the following contrasts:
1. Focal continental group vs. all other groups (e.g. AFR vs. non-AFR)
2. Focal population vs. all other populations (within the same continental group, e.g. YRI vs. non-YRI-in-AFR)

The environment used to perform this analysis are detailed in the [`03_deseq.yaml`](resources/03_deseq.yaml) anaconda environment recipe, and the input parameters are available in the [`03_deseq_config.yaml`](resources/03_deseq_config.yaml) config file. The inputs for this notebook are the following:
1. Directory of salmon quantifications generated from [`01_run_salmon.sh`](01_run_salmon.sh)
2. BED file containing genes retained after filtering for low-expression or low-complexity (generated from [`02_filter_normalize.sh`](02_filter_normalize.sh))
3. Gencode annotation GFF3 file (available from [here](https://www.gencodegenes.org/human/release_38.html))
4. Metadata file containing sample information

