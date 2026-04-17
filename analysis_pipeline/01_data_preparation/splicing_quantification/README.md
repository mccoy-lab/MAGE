# MAGE Splicing Quantification

This directory contains code to replicate the splicing quantification procedure done as part of the MAGE v.1.1 release. This procedure is broadly split into two sections:
1. Variant-aware read-mapping with [STAR](https://github.com/alexdobin/STAR)
2. Annotation-agnostic estimation of intron excision with [Leafcutter](https://davidaknowles.github.io/leafcutter/) and subsequent filtering and normalization for sQTL mapping analyses

*Note: these steps should be run in order*<br><br>

## Variant-aware read mapping with STAR

The `01_align_dev.sh` script uses VCF files from 1000 Genomes samples to perform variant-aware alignment with STAR. The script assumes that `bcftools`, `htslib`, `vcftools`, and `STAR` are in the user's path. The script requires seven arguments:

1. `workDir`: Working directory where outputs will be generated
2. `metaData`: The path to the `updated_metadata.tsv` file, provided in the `scripts/` directory
3. `sampleID`: The identifier of a sample to align, such as `GM18861` or `HG00148`; replicates are denoted based on batch and replicate number, such as `GM20815-13-1` and `GM20815-13-2`
4. `threads`: The number of threads to use for computation (set to 24 for our analysis)
5. `kgpVCFdir`: Directory containing per-chromosome 1000 Genomes VCF files with names such as `CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.eagle2-phased.v2.vcf.gz`
6. `starIndex`: Path to the STAR index directory (generated from [`GCA_000001405.15_GRCh38_no_alt_analysis_set.fna`](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz))
7. `fastqDir`: Path to directory containing FASTQ files, organized into subdirectories based on sequencing batch

Choices of parameters were inspired by code from [TOPMed](https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md). The script outputs BAM alignment files.

## Splicing quantification with Leafcutter

The [`02_intron_usage/`](./02_intron_usage/) directory contains a snakemake pipeline for quantifying intron excision ratios with [Leafcutter](https://davidaknowles.github.io/leafcutter/) using the `STAR` alignments generated in the previous section. Intron usage estimates are aggregated across all samples and then filtered for clusters bearing either low-expression or low-complexity.

The primary outputs of this pipeline are the following:
1. Raw filtered intron excision ratios (no normalization)
2. Normalized filtered intron excision ratios
3. A file mapping introns to genes

Technical details on the pipeline are documented in [`02_intron_usage/README.md`](./02_intron_usage/README.md).