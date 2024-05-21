#!/usr/bin/bash

#==========#
# Get args #
#==========#

nominalResultsFile=$1 
permutationResultsFile=$2
out_eGeneListFile=$3 # Output file containing significant eGenes
out_sigNominalResultsFile=$4 # Output file contains significant nominal associations

FDR=0.05


#=============#
# Set scripts #
#=============#

# Get the directory of this script
thisDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

scriptDir=$thisDir/scripts

getSigResultsScript=$scriptDir/get_sig_results.py


#=========================#
# Get significant results #
#=========================#

gzip_check=False

if [[ $out_sigNominalResultsFile == *.gz ]]; then
	out_sigNominalResultsFile=${out_sigNominalResultsFile%".gz"}
	gzip_check=True
fi

$getSigResultsScript --nominalResults $nominalResultsFile \
					 --permResults $permutationResultsFile \
					 --out_eGene $out_eGeneListFile \
					 --out_sigResults $out_sigNominalResultsFile \
					 --FDR $FDR

if [[ $gzip_check == "True" ]]; then
	gzip $out_sigNominalResultsFile
fi
