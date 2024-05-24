#!/usr/bin/env bash

# Bash Script for VEP Cache Building and sQTL Data VEP Annotation

# Description:
# This script builds VEP cache for GRCh38 and performs sQTL variant annotation using the VEP tool.

# Example Usage:
# ./script.sh < <vcf_file> <output_file> <cache_directory> 

# Check for required arguments
if [ "$#" -ne 3 ]; then
    echo "Incorrect usage!"
    echo "Usage: $0 <vcf_file> <output_file> <cache_directory>"
    exit 1
fi


# User-specified inputs
vcf_file=$1
output_file=$2
cache_dir=$3


# Initialize and activate conda environment
initialize_conda() {
    echo "Initializing Conda environment..."
    eval "$(conda shell.bash hook)"
    conda activate vep_env
    echo "Conda environment activated."
}


# Function to build VEP cache for GRCh38
build_vep_cache() {
    echo -e "\nBuilding VEP cache for GRCh38 ....\n"

    INSTALL.pl --CACHEDIR $cache_dir \
               -a cfp \
               -s homo_sapiens \
               -y GRCh38 \
               --CONVERT \
               --PLUGINS TSSDistance,LOEUF,LoF,LoFtool,pLI,Conservation,SpliceAI,SpliceRegion,LOFTEE \
               -PLUGINSDIR $cache_dir/Plugins/
}


# Function to annotate sQTLs with VEP
annotate_sqtls_with_vep() {
    echo -e "\nAnnotation of sQTLs with VEP ....\n"

    vep --input_file $vcf_file \
        --format vcf \
        --output_file $output_file \
        --vcf \
        --cache \
        --regulatory \
        --fork 4 \
        --dir_cache $cache_dir \
        --plugin LoF,loftee_path:$cache_dir/Plugins/loftee \
        --dir_plugins $cache_dir/Plugins/loftee
}


# Call functions
initialize_conda
build_vep_cache
annotate_sqtls_with_vep

