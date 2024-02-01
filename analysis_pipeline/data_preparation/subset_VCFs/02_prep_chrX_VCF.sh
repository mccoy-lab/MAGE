#!/usr/bin/bash

#==========#
# Get args #
#==========#

inVCF=$1 # chrX VCF from 1KGP
sampleListFile=$2 # Text file with selected samples
outVCFgz=$3 # Output VCF, should end in .gz


#=============#
# Set scripts #
#=============#

# Get the directory of this script
thisDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

scriptDir=$thisDir/scripts

hap2dipScript=$scriptDir/haploid2diploid.py


#==============#
# Get chrX VCF #
#==============#

tmpVCF=$(mktemp --suffix=".vcf.gz")

bcftools view -S $sampleListFile --min-ac=1 -O z -o $tmpVCF $inVCF


#================================================#
# Convert haploid genotypes to diploid genotypes #
#================================================#

tmpVCF2=$(mktemp --suffix=".vcf")

$hap2dipScript $tmpVCF > $tmpVCF2

bgzip $tmpVCF2
mv $tmpVCF2.gz $outVCFgz
tabix -p vcf $outVCFgz


#==========#
# Clean-up #
#==========#

rm $tmpVCF
