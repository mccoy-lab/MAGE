#!/usr/bin/env bash

# Inputs
metadata=$1
vcf=$2
finemapped_variants=$3
chrom_str=$4
output_dir=$5

# Use chrom_str to create a temporary txt file with finemapped loci for chrom_str.
# Note that the finemapped_variants file is gzipped and tabix indexed.
echo "Extracting finemapped loci for ${chrom_str}"
finemapped_loci=$(mktemp)
tabix -h $finemapped_variants ${chrom_str} |\
    cut -f5 > $finemapped_loci #This output will be removed by Snakemake by marking as temp()

# Subset vcf to finemapped loci and index
# Note that using `ID==@$finemapped_loci` will subset to match loci and variant states.
echo "Subsetting VCF to finemapped loci for ${chrom_str}"
subset_vcf=${output_dir}finemapped_loci.${chrom_str}.vcf.gz

bcftools view --include ID==@$finemapped_loci -Oz -o $subset_vcf $vcf
bcftools index --tbi $subset_vcf
echo "Finemapped loci written to ${subset_vcf}"

# Get sample list from vcf
sample_list=$(bcftools query -l "$subset_vcf")

# Loop through continental group labels
for group in AFR AMR EAS EUR SAS; do
    echo "Processing group: $group"
    
    # Extract sample IDs for the current group
    focal_samples=$(awk -v grp="$group" '$6 == grp {print $4}' "$metadata")
    background_samples=$(awk -v grp="$group" '$6 != grp {print $4}' "$metadata")

    # Get the index of each element in focal_samples and 
    # background_samples in the sample_list
    focal_sample_indices=$(for id in $focal_samples; do
        echo "${sample_list}" | grep -n $id | cut -d: -f1
    done | sort -n | paste -sd "," -)

    background_sample_indices=$(for id in $background_samples; do
        echo "${sample_list}" | grep -n $id | cut -d: -f1
    done | sort -n | paste -sd "," -)
    
    echo "Focal sample indices: $focal_sample_indices"
    echo "Background sample indices: $background_sample_indices"

    # Perform Fst calculation using wcFst from vcflib
    wcFst --file $subset_vcf \
        --target $focal_sample_indices \
        --background $background_sample_indices \
        --type GT > ${output_dir}wcFst.${group}.${chrom_str}.tsv
done