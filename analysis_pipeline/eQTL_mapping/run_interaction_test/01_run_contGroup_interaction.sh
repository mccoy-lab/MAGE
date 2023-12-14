#!/usr/bin/bash

#==========#
# Get args #
#==========#

eQTLFile=$1 # Four-column tab-delimted file of eQTLs to test (the same file used as input to aFCn)
inVCF=$2 # Same VCF file used for eQTL mapping
expressionBed=$3 # TSS BED with inv-norm transformed TMM values (same file used for eQTL mapping)
populationTable=$4 # Two column table with sampleIDs and corresponding population label
covFile=$5 # Same file of covariates used for eQTL mapping
outFileSingle=$6 # File to write single-variant output results to
outFileMulti=$7 # file to write multi-variant output results to

mafThreshold=0.05 # Variants must have MAF above this threshold in at least two populations to be tested


#=============#
# Set scripts #
#=============#

# Get the directory of this script
thisDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

scriptDir=$thisDir/scripts

interactionScript=$scriptDir/QTL_interaction.py
multivariantInteractionScript=$scriptDir/QTL_interaction.multivariant.py


#=========================#
# Format expression input #
#=========================#

tmpExpressionFile=$(mktemp --suffix=".expression.txt")

zcat $expressionBed | cut -f 4- | awk -v FS="\t" -v OFS="," 'NR==1 {$1="Name"} {$1=$1; print}' | gzip -c > $tmpExpressionFile


#=====================================#
# Run single-variant interaction test #
#=====================================#

$interactionScript $eQTLFile $inVCF $tmpExpressionFile $populationTable $outFileSingle --covFile $covFile --maf_threshold $mafThreshold


#====================================#
# Run multi-variant interaction test #
#====================================#

$multivariantInteractionScript $eQTLFile $inVCF $tmpExpressionFile $populationTable $outFileMulti --covFile $covFile --maf_threshold $mafThreshold


#==========#
# Clean-up #
#==========#

rm $tmpExpressionFile
