#!/usr/bin/bash

#==========#
# Get args #
#==========#

inVCF=$1 # VCF file with sample genotypes
invNormTMMBed=$2 # TSS Bed file with inverse normal transformed TMM values for samples
covFile=$3 # Formatted covariates file
outDir=$4 # Directory to write output files to
outPrefix=$5 # Prefix of output files

window=1e6 # window up- and downstream of TSS to test for associations 
mafThreshold=0.01 # MAF threshold to test for associations
chunks=1 # This assumes each chromosome being run separately. If together, this will need to be set to the number of chromosomes (as a minimum)
threads=1


#=============#
# Run FastQTL #
#=============#

run_FastQTL_threaded.py $inVCF $invNormTMMBed $outPrefix \
						--covariates $covFile \
						-o $outDir \
						--window $window \
						--maf_threshold $mafThreshold \
						--chunks $chunks \
						--threads $chunks
