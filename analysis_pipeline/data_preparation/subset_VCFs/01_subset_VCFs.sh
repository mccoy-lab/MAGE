#!/usr/bin/bash

#==========#
# Get args #
#==========#

in_vcfDir=$1 #Directory containing chromosome-split VCFs from 1KGP
sampleListFile=$2 # Text file with selected samples
out_vcfDir=$3 #Output VCF directory


#=========================#
# Create output directory #
#=========================#

if [[ ! -d $out_vcfDir ]]; then
	mkdir -p $out_vcfDir
fi


#====================#
# Get autosomal VCFs #
#====================#

for i in {1..22}; do
	chrom=chr${i}
	inVCF=$in_vcfDir/CCDG_14151_B01_GRM_WGS_2020-08-05_${chrom}.filtered.shapeit2-duohmm-phased.vcf.gz
	outVCF=$out_vcfDir/$chrom.MAGE.v0.1.vcf.gz

	bcftools view -S $sampleListFile --min-ac=1 -O z -o $outVCF $inVCF
	tabix -p vcf $outFile
done


#======================#
# Merge autosomal VCFs #
#======================#

outMergedVCF=$out_vcfDir/autosomal.MAGE.v0.1.vcf.gz

inputString=""
for i in {1..22}; do
	inputString="$inputString $out_vcfDir/$chrom.MAGE.v0.1.vcf.gz"
done

bcftools concat -o $outMergedVCF -O z $inputString
tabix -p vcf $outMergedVCF
