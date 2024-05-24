#!/usr/bin/env bash

# Bash Script Pipeline for eQTL Data Processing and GREGOR Analysis

# Description:
# This script processes the eQTL-input files, preparation of input files for GREGOR, 
# updates to configuration files, and the execution of GREGOR for enrichment analysis.

# Example Usage:
# ./script.sh <input_directory> <output_directory> <config_base_directory> <work_directory>

# Exit if not all required arguments are provided
if [ "$#" -ne 4 ]; then
    echo "Incorrect usage!"
    echo "Usage: $0 <input_directory> <output_directory> <config_base_directory> <work_directory>"
    exit 1
fi

# User-specified directories
input_dir=$1
out_dir=$2
config_base_dir=$3
work_dir=$4


# Initialize and activate conda environment
initialize_conda() {
    echo "Initializing Conda environment..."
    eval "$(conda shell.bash hook)"
    conda activate perl_env
    echo "Conda environment activated."
}


# Process input files and convert them to gregor format
process_input_files() {
    echo "Processing input files..."
    local input_file="${input_dir}/pip95.finalmage_eqtl.txt"
    local output_filename="${input_file%.txt}.gregformat.str.txt"

    awk -F'\t' 'NR > 1 {
        split($2, arr, ":");
        print "chr"arr[1], arr[2] - 1, arr[2], arr[3], arr[4], $2, $1, $3, $4
    }' OFS='\t' $input_file | sort -k1,1 -k2,2n > $output_filename

    echo "Output prepared at $output_filename"
}


# Generate a unique SNP index and input the gregor format file
generate_snp_index() {
    echo "Generating unique SNP index..."
    mkdir -p $out_dir
    local input_file="${input_dir}/snps.index.finalmage_eqtl.gregformat.str.txt"
    local base_file=$(basename "$input_file")
    local new_file=$(echo "$base_file" | sed -E 's/finalmage_eqtl/finalmage_eqtluniq/')

    sort $input_file | uniq > "${out_dir}/${new_file}"
    echo "SNP index generated at ${out_dir}/${new_file}"
}


# Update configuration files with the new file paths
update_config_files() {
    echo "Updating configuration files with new paths..."
    local source_string="INDEX_SNP_FILE = ${out_dir}/snps.index.finalmage_eqtl"
    local target_string="INDEX_SNP_FILE = ${out_dir}/snps.index.finalmage_eqtluniq"

    for file in "${config_base_dir}"/*.conf; do
        if [ -f "$file" ]; then
            sed -i "s|$source_string|$target_string|g" "$file"
        fi
    done
    echo "Configuration files updated."
}


# Submit GREGOR analysis jobs using slurm
submit_gregor_jobs() {
    echo "Submitting GREGOR analysis jobs..."
    local gregor_script="${work_dir}/script/GREGOR.pl"
    local datavar="mage_eqtl"

    # Ensure log directory exists
    mkdir -p "${work_dir}/logfiles"

    # Loop to submit jobs for different configurations
    for i in {1..8}; do
        local conf_file="${work_dir}/conf/${datavar}/config_set${i}.conf"
        sbatch -p defq -n 12 --mem 40G -t 72:00:00 \
               -o "${work_dir}/logfiles/${datavar}log${i}.out" \
               -e "${work_dir}/logfiles/${datavar}log${i}.err" \
               -J "${datavar} ${i}_ENRICH ANALYSIS" \
               --wrap="perl $gregor_script --conf $conf_file"
    done
    echo "GREGOR analysis jobs have been submitted."
}


# Call functions
initialize_conda
process_input_files
generate_snp_index
update_config_files
submit_gregor_jobs

