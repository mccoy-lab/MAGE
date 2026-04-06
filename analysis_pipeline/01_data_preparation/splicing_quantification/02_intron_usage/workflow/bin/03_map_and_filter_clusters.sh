#!/usr/bin/bash

#==========#
# Get args #
#==========#

gencodeGTF=$1 # Gencode GTF (not gzipped)
leafcutterOutDir=$2 # Path to Leafcutter output directory
leafcutterOutPrefix=$3 # Prefix of Leafcutter output files (i.e. everything before `_perind_numers.counts.gz` and `_perind.counts.gz`)
sampleLookupFile=$4
selectedSampleList=$5
leafcutterPrepHelperScript=$6 # Path to leafcutter `prepare_phenotype_table.py` helper script
outDir=$7


#=============#
# Set scripts #
#=============#

# Get the directory of this script
thisDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

scriptDir=$thisDir/scripts

extractExonScript=$scriptDir/extract_exons_from_gtf.py
mappingScript=$scriptDir/map_clusters_to_genes.R
filterScript=$scriptDir/filter_quants.py


#=========================#
# Create output directory #
#=========================#

if [[ ! -d $outDir ]]; then
	mkdir -p $outDir
fi


#========================#
# Extract exons from GTF #
#========================#

exonFile=$(mktemp)

$extractExonScript $gencodeGTF --outputFile $exonFile


#================================#
# Map splicing clusters to genes #
#================================#

$mappingScript $leafcutterOutDir/${leafcutterOutPrefix}_perind.counts.gz $exonFile $outDir/clusters_to_genes.txt


#======================================#
# Filter and normalize splicing quants #
#======================================#

$filterScript $leafcutterOutDir \
			  $leafcutterOutPrefix \
			  $sampleLookupFile \
			  $gencodeGTF \
			  $outDir/clusters_to_genes.txt \
			  --keep_samples $selectedSampleList \
			  --pct_cutoff 10 \
			  --leafcutter_prep_script $leafcutterPrepHelperScript \
			  --output_dir $outDir \
			  --autosomes_chrX


#==========#
# Clean up #
#==========#

rm $exonFile
