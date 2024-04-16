#!/usr/bin/bash

#==========#
# Get args #
#==========#

mapFile=$1 # Four-column tab-separated file with trait accession, trait description, sentinel variant, coloc results prefix
outPrefix=$2 # Prefix of output files


#=============#
# Set scripts #
#=============#

# Get the directory of this script
thisDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
scriptDir=$thisDir/scripts

collateScript=$scriptDir/collate_coloc_results.py


#=======================#
# Collate coloc results #
#=======================#

$collateScript $mapFile $outPrefix
