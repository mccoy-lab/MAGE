# MAGE Splicing Quantification

This directory contains code to replicate the splicing quantification procedure done as part of the initial MAGE publication. This procedure is broadly split into three sections:
1. Variant-aware read mapping with STAR
2. Annotation-agnostic estimation of intron excision with [Leafcutter](https://davidaknowles.github.io/leafcutter/)
3. Filtering and normalization for sQTL analsyses

*Note: these steps should be run in order*<br><br>

## Variant-aware read mapping with STAR

Work-in-progress
<br><br>

## Splicing quantification with Leafcutter

Work-in-progress
<br><br>

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
