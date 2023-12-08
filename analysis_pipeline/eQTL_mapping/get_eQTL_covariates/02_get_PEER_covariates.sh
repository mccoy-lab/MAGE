#!/usr/bin/bash

#==========#
# Get args #
#==========#

invNormTMMBed=$1 # Inverse normal transformed TMM values BED file (autosomal + chrX)
outDir=$2 # Directory to write PEER output to

numFactors=60


#=========================#
# Create output directory #
#=========================#

if [[ ! -d $outDir ]]; then
	mkdir -p $outDir
fi


#===================#
# Format expression #
#===================#

tmpExpression=$(mktemp --suffix .expression.tab) # .tab suffix is NECESSARY for peertool to function properly

printf "\t" > $tmpExpression
zcat $invNormTMMBed | head -n1 | cut -f 5- >> $tmpExpression
zcat $invNormTMMBed | tail -n +2 | cut -f 4- >> $tmpExpression


#=====================#
# Get PEER covariates #
#=====================#

peertool --transpose \
		 --has_rownames \
		 --has_header \
		 -i 100 \
		 -f $tmpExpression \
		 -n $numFactors \
		 -o $outDir


#==========#
# Clean up #
#==========#

rm $tmpExpression
