#!/usr/bin/bash

#==========#
# Get args #
#==========#

inVCF=$1 # chrX VCF from 1KGP
sampleListFile=$2 # Text file with selected samples
out_chrX_VCFgz=$3 # Output chrX VCF (no genotype conversion), should end in .gz
out_chrX_hap2dip_VCFgz=$4 # Output chrX VCF with haploid genotypes converted to diploid genotypes, should end in .gz


#=============#
# Set scripts #
#=============#

# Get the directory of this script
thisDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

scriptDir=$thisDir/scripts

addIDsScript=$scriptDir/add_ids.py
hap2dipScript=$scriptDir/haploid2diploid.py


#==============#
# Get chrX VCF #
#==============#

tmpVCF=$(mktemp --suffix=".vcf.gz")

bcftools view -S $sampleListFile --min-ac=1 -O z -o $tmpVCF $inVCF


#=================#
# Add variant IDs #
#=================#

tmpVCF2=$(mktemp --suffix=".vcf")

$addIDsScript $tmpVCF > $tmpVCF2

bgzip $tmpVCF2
mv $tmpVCF2.gz $out_chrX_VCFgz
tabix -p vcf $out_chrX_VCFgz

rm $tmpVCF


#================================================#
# Convert haploid genotypes to diploid genotypes #
#================================================#

tmpVCF3=$(mktemp --suffix=".vcf")

$hap2dipScript $out_chrX_VCFgz > $tmpVCF3

bgzip $tmpVCF3
mv $tmpVCF3.gz $out_chrX_hap2dip_VCFgz
tabix -p vcf $out_chrX_hap2dip_VCFgz
