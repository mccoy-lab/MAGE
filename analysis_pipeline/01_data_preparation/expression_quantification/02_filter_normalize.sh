#!/usr/bin/bash

#==========#
# Get args #
#==========#

transcriptGTF=$1 # GTF file with transcript annotations (e.g. gencode.v38.annotation.gtf.gz)
salmonQuantListFile=$2 # Two-column tab-separated file with sampleID and path to salmon quant.sf file
outQuantsDir=$3 # Directory to write gene-level quantifications to
outPrefix=$4 # Prefix of output files


#=============#
# Set scripts #
#=============#

# Get the directory of this script
thisDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

scriptDir=$thisDir/scripts

mappingScript=$scriptDir/map_transcript_to_gene.py
getGeneTssBedScript=$scriptDir/getGeneTSSBed.py
normFilterScript=$scriptDir/normalize_and_filter.R
invNormTransformScript=$scriptDir/invNormTransform.py
tableToBedScript=$scriptDir/table_to_bed.py


#=========================#
# Create output directory #
#=========================#

if [[ ! -d $outQuantsDir ]]; then
	mkdir -p $outQuantsDir
fi


#=====================#
# Create mapping file #
#=====================#

tmpMapFile=$(mktemp)
$mappingScript --gtf $transcriptGTF --out $tmpMapFile


#============================#
# Get autosomal + chrX genes #
#============================#

tmpGeneTSSBed=$(mktemp)
$getGeneTssBedScript --gtf $transcriptGTF --out $tmpGeneTSSBed

tmpGeneListFile=$(mktemp)
for i in {1..22} X; do
	chrom=chr${i}
	cat $tmpGeneTSSBed | awk -v chrom=$chrom '$1==chrom {print $4}' >> $tmpGeneListFile
done


#======================================#
# Merge transcripts, normalize, filter #
#======================================#

$normFilterScript --salmonFiles $salmonQuantListFile \
				  --tx2gene $tmpMapFile \
				  --keepGenes $tmpGeneListFile \
				  --outPrefix $outQuantsDir/$outPrefix


#==========================#
# Inverse normal transform #
#==========================#

$invNormTransformScript --tmm $outQuantsDir/$outPrefix.filtered.tmm.tab \
						--out $outQuantsDir/$outPrefix.filtered.inverse-normal-tmm.tab


#===================#
# Generate TSS Beds #
#===================#

for suffix in "DESeq2-normalized-logged-pseudocounts" "inverse-normal-tmm" "pseudocounts" "tmm"; do
	expressionTable=$outQuantsDir/$outPrefix.filtered.$suffix.tab
	outBed=$outQuantsDir/$outPrefix.filtered.$suffix.bed

	$tableToBedScript --inputTable $expressionTable --gtf $transcriptGTF --out $outBed

	bgzip $outBed
	tabix -p bed $outBed.gz

done


#==========#
# Clean up #
#==========#

rm $tmpMapFile $tmpGeneTSSBed $tmpGeneListFile
