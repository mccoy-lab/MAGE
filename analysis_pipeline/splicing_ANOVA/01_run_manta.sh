#!/usr/bin/bash

#==========#
# Get args #
#==========#

ratiosFile=$1 # Filtered splicing ratios file from Leafcutter
metadataFile=$2 # Tab-separated file with sample metadata: (at least batch, sex, continentalGroup, and population)
outFile=$3 # File to write MANTA results to

threads=16


#=============#
# Set scripts #
#=============#

# Get the directory of this script
thisDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

scriptDir=$thisDir/scripts

mantaScript=$scriptDir/run_manta.R


#===========#
# Run MANTA #
#===========#

$mantaScript --splicing_frac $ratiosFile --covariates $metadataFile --threads $threads --output $outFile
