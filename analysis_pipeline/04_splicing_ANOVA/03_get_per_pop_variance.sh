#!/usr/bin/bash

#==========#
# Get args #
#==========#

ratiosFile=$1 # Filtered splicing ratios file from Leafcutter
metadataFile=$2 # Tab-separated file with sample metadata: (at least batch, sex, continentalGroup, and population)
outContGroupFile=$3 # File to write within-continental group variance results to
outPopFile=$4 # File to write within-population variance results to

threads=16


#=============#
# Set scripts #
#=============#

# Get the directory of this script
thisDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

scriptDir=$thisDir/scripts

mantaScript=$scriptDir/per_pop_variance.R


#===========#
# Run MANTA #
#===========#

$mantaScript --splicing_frac $ratiosFile --covariates $metadataFile --threads $threads --contGroupOut $outContGroupFile --popOut $outPopFile
