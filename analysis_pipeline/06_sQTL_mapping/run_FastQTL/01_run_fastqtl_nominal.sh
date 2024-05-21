#!/usr/bin/bash

#==========#
# Get args #
#==========#

inVCF=$1 # VCF file with sample genotypes
normedRatiosBed=$2 # Bed file with filtered, normalized splicing ratios
covFile=$3 # Formatted covariates file
phenGroupsFile=$4
outDir=$5 # Directory to write output files to
outPrefix=$6 # Prefix of output files

window=1e6 # window up- and downstream of TSS to test for associations 
mafThreshold=0.01 # MAF threshold to test for associations
chunks=1 # This assumes each chromosome being run separately. If together, this will need to be set to the number of chromosomes (as a minimum)
threads=1


#=============#
# Run FastQTL #
#=============#

run_FastQTL_threaded.py $inVCF $normedRatiosBed $outPrefix \
						--covariates $covFile \
						--phenotype_groups $phenGroupsFile \
						-o $outDir \
						--window $window \
						--maf_threshold $mafThreshold \
						--chunks $chunks \
						--threads $chunks
