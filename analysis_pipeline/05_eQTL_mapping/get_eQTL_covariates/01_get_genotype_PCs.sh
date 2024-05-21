#!/usr/bin/bash

#==========#
# Get args #
#==========#

inVCF=$1 # autosomal VCF (for samples to be used in eQTL mapping)
outDir=$2 # Directory to write plink PCA files to
outPrefix=$3 # Prefix for plink output files

MAF_threshold=0.01
num_PCs=20


#=========================#
# Create output directory #
#=========================#

if [[ ! -d $outDir ]]; then
	mkdir -p $outDir
fi


#=============================================#
# Get top genotype PCs and relatedness matrix #
#=============================================#

plink --vcf $sampleVCF \
--pca $num_PCs \
--make-rel \
--maf $MAF_threshold \
--out $outDir/$outPrefix
