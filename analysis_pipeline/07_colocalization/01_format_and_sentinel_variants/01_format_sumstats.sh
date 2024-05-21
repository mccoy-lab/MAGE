#!/usr/bin/bash

#==========#
# Get args #
#==========#

inFile=$1 # Path to harmonized summary statistics from the GWAS catalog. The filename should end with ".h.tsv.gz"\
outFile=$2 # Path to output file

outFile=${outFile%".gz"} # Remove ".gz" suffix if necessary


#=============#
# Set scripts #
#=============#

# Get the directory of this script
thisDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
scriptDir=$thisDir/scripts

formatScript=$scriptDir/format_sumstats.py


#=================#
# Format sumstats #
#=================#

$formatScript $inFile $outFile

bgzip -f $outFile
tabix -f -S1 -b2 -e2 $outFile.gz
