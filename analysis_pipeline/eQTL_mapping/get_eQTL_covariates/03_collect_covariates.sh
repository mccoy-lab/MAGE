#!/usr/bin/bash

#==========#
# Get args #
#==========#

pcaEigenvec=$1 # output .eigenvec file from plink
peerX=$2 # X.csv file from PEER
addCovFile=$3 # File with additional covariates to include (e.g. sex)
outFile=$4 # File to write formatted covariates file to

numPCs=5

# Note that the PEER X.csv file does not have sample IDs.
# This script assumes that the samples in X.csv are in the same order as
# in the plink .eigenvec file (i.e. that the sample order in the expression
# bed matches the sample order in the VCF used with plink)


#=============#
# Set scripts #
#=============#

# Get the directory of this script
thisDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

scriptDir=$thisDir/scripts

collectionScript=$scriptDir/collect_covariates.py


#=====================#
# Subset Genotype PCs #
#=====================#

tmpPCs=$(mktemp)

numCols=$(expr $numPCs + 2)

cat $pcaEigenvec | cut -d' ' -f -$numCols > $tmpPCs


#====================#
# Collect covariates #
#====================#

$collectionScript --eigenvecs $tmpPCs \
				  --peerX $peerX \
				  --addCov $addCovFile \
				  --out $outFile


#==========#
# Clean up #
#==========#

rm $tmpPCs
