#!/usr/bin/bash

#==========#
# Get args #
#==========#

normedRatiosBed=$1 # Bed file with filtered, normalized splicing ratios (autosomal + chrX)
outDir=$2 # Directory to write PEER output to

numFactors=15


#=========================#
# Create output directory #
#=========================#

if [[ ! -d $outDir ]]; then
	mkdir -p $outDir
fi


#===================#
# Format expression #
#===================#

tmpSplicing=$(mktemp --suffix .splicing.tab) # .tab suffix is NECESSARY for peertool to function properly

printf "\t" > $tmpSplicing
zcat $normedRatiosBed | head -n1 | cut -f 5- >> $tmpSplicing
zcat $normedRatiosBed | tail -n +2 | cut -f 4- >> $tmpSplicing


#=====================#
# Get PEER covariates #
#=====================#

peertool --transpose \
		 --has_rownames \
		 --has_header \
		 -i 100 \
		 -f $tmpSplicing \
		 -n $numFactors \
		 -o $outDir


#==========#
# Clean up #
#==========#

rm $tmpSplicing
