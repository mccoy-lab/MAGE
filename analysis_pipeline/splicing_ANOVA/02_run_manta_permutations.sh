#!/usr/bin/bash

#==========#
# Get args #
#==========#

ratiosFile=$1 # Filtered splicing ratios file from Leafcutter
metadataFile=$2 # Tab-separated file with sample metadata: (at least batch, sex, continentalGroup, and population)
nPerms=$3 # Number of permutations to run for each cluster
outDir=$4 # Directory to write MANTA permutation results to
outPrefix=$5 # Prefix of output manta files

perms=1000
threads=16


#=============#
# Set scripts #
#=============#

# Get the directory of this script
thisDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

scriptDir=$thisDir/scripts

mantaPermScript=$scriptDir/run_manta_permutations.R


#=========================#
# Create output directory #
#=========================#

if [[ ! -d $outDir ]]; then
	mkdir -p $outDir
fi

#===========#
# Run MANTA #
#===========#

$mantaPermScript --splicing_frac $ratiosFile --covariates $metadataFile --perms $perms --threads $threads --outPrefix $outDir/$outPrefix
