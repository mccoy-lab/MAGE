# MAGE Splicing Quantification

This directory contains code to replicate the splicing quantification procedure done as part of the initial MAGE publication. This procedure is broadly split into three sections:
1. Variant-aware read mapping with STAR
2. Annotation-agnostic estimation of intron excision with [Leafcutter](https://davidaknowles.github.io/leafcutter/)
3. Filtering and normalization for sQTL analsyses

*Note: these steps should be run in order*<br><br>

## Variant-aware read mapping with STAR

The `01_align_dev.sh` script uses VCF files from 1000 Genomes samples to perform variant-aware alignment with STAR. The script assumes that `bcftools`, `htslib`, `vcftools`, and `STAR` are in the user's path. The script requires seven arguments:

1. `workDir`: Working directory where outputs will be generated
2. `metaData`: The path to the `updated_metadata.tsv` file, provided in this directory
3. `sampleID`: The identifier of a sample to align, such as `GM18861` or `HG00148`; replicates are denoted based on batch and replicate number, such as `GM20815-13-1` and `GM20815-13-2`
4. `threads`: The number of threads to use for computation (set to 24 for our analysis)
5. `kgpVCFdir`: Directory containing per-chromosome 1000 Genomes VCF files with names such as `CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.eagle2-phased.v2.vcf.gz`
6. `starIndex`: Path to the STAR index directory (generated from [`GCA_000001405.15_GRCh38_no_alt_analysis_set.fna`](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz))
7. `fastqDir`: Path to directory containing FASTQ files, organized into subdirectories based on sequencing batch

Choices of parameters were inspired by code from [TOPMed](https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md). The script outputs BAM alignment files.

## Splicing quantification with Leafcutter

The `02_intron_usage.sh` script will use the `STAR` alignments generated in the previous section to estimate intron cluster usage per sample. Following, intron usage estimates are aggregated across all samples and consolidated into a single file. The focal outputs from this script are `*._perind_numers.counts.gz` and `*perind.counts.gz` which contain ***counts*** of intron cluster usage and ***ratios***, respectively.<br>

The script takes a single argument that is the path to the STAR output directory containing BAM files and their respective indices.<br>
*Before running this script, install `leafcutter` following the instructions [here](https://davidaknowles.github.io/leafcutter/articles/Installation.html).*<br><br>

## Filtering and normalization for sQTL analyses

The `03_map_and_filter_clusters.sh` script will filter out lowly expressed and low-complexity splicing clusters. It will output three files:
1. Raw filtered intron excision ratios (no normalization)
2. Normalized filtered intron excision ratios
3. A file mapping introns to genes

The script takes seven arguments:
1. `gencodeGTF`: The Gencode transcript GTF file, used for mapping introns to genes (we used the Gencode v38 files [here](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz))
2. `leafcutterOutDir`: The directory containing leafcutter files from the previous step
3. `leafcutterOutPrefix`: Prefix of Leafcutter output files (i.e. everything before `_perind_numers.counts.gz` and `_perind.counts.gz`)
4. `sampleLookupFile`: A file mapping columns of the leafcutter files to sampleIDs (to be included in the output file). Should be a two-column tab-separated file
5. `selectedSampleList`: A file listing the samples to include in the output file. If one column, will subset to the corresponding samples (after renaming from the previous argument). If two columns, will subset and re-name. We used this to subset to only one replicate for samples that were sequenced in triplicate.
6. `leafcutterPrepHelperScript`: Path to the Leafcutter `prepare_phenotype_table.py` helper script, which is used for normalization (and some filtering)
7. `outDir`: Path to directory to write output files to
<br><br>
