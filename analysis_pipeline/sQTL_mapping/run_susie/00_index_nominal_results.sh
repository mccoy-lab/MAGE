#!/usr/bin/bash

#==========#
# Get args #
#==========#

nominalResultsFile=$1 # Nominal results file from FastQTL
normedRatiosBed=$2 # The same filtered, normalized splicing ratios bed file used with FastQTL (or, at least the first four columns)
outSortedFile=$3 # File to write position sorted results to. Will be gzipped, so include ".gz"


#=============#
# Set scripts #
#=============#

# Get the directory of this script
thisDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

scriptDir=$thisDir/scripts

bedScript=$scriptDir/make_bed.py


#===========================================#
# Add intron location info, sort, and index #
#===========================================#

tmpOutFile=$(mktemp --suffix ".txt")

$bedScript $nominalResultsFile $normedRatiosBed $tmpOutFile

bgzip $tmpOutFile
mv $tmpOutFile.gz $outSortedFile

tabix -p bed -S1 $outSortedFile
