#!/usr/bin/bash

#==========#
# Get args #
#==========#

formattedSumstatsFile=$1 # Formatted sumstats file from the previous script
outSentinelVariantsFile=$2 # Output file with sentinel variants

outSentinelVariantsFile=${outSentinelVariantsFile%".gz"} # Remove ".gz" suffix if necessary

# Global variables
colocWindow=500000 # +/- 500kb around each sentinel variant will be tested for colocalization
sentinelVariantSigCutoff=5e-8


#=============#
# Set scripts #
#=============#

# Get the directory of this script
thisDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
scriptDir=$thisDir/scripts

sentinelVariantScript=$scriptDir/get_sentinel_variants.py


#=======================#
# Get sentinel variants #
#=======================#

$sentinelVariantScript $formattedSumstatsFile $((2*colocWindow)) $sentinelVariantSigCutoff $outSentinelVariantsFile

bgzip -f $outFile
