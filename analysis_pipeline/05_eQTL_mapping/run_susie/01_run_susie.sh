#!/usr/bin/bash

#==========#
# Get args #
#==========#

geneID=$1 # This script is meant to be run with a single gene. It can be parallelized using something like SLURM
inVCF=$2 # The same VCF used with FastQTL
invNormTMMBed=$3 # The same TSS Bed file used with FastQTL
covFile=$4 # The same covariates file used with FastQTL
nomResultsFile=$5 # The nominal pass output from FastQTL
outDir=$6

window=1e6 # window up- and downstream of TSS to test for associations 
mafThreshold=0.01 # MAF threshold to test for associations


#=============#
# Set scripts #
#=============#

# Get the directory of this script
thisDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

scriptDir=$thisDir/scripts

susieScript=$scriptDir/run_susie.R


#======================#
# Get variants to test #
#======================#

tmpPrefix=$(mktemp)
tmpGeneBed=$tmpPrefix.gene.bed
tmpVariantsPrefix=$tmpPrefix.variants

zgrep -m1 $geneID $invNormTMMBed | awk -v window=$window 'BEGIN{OFS="\t";} {print $1, ($2-window>0)?$2-window:0, $3+window}' > $tmpGeneBed
bcftools view --regions-overlap pos -R $tmpGeneBed --min-af ${mafThreshold}:minor $inVCF > $tmpVariantsPrefix.vcf
vcftools --gzvcf $tmpVariantsPrefix.vcf --012 --out $tmpVariantsPrefix 2> /dev/null
grep ^[^#] $tmpVariantsPrefix.vcf | cut -f 3 > $tmpVariantsPrefix.012.snps

# Clean-up unneeded files
rm $tmpPrefix
rm $tmpGeneBed
rm $tmpVariantsPrefix.vcf
rm $tmpVariantsPrefix.012.pos


#============================#
# Extract expression of gene #
#============================#

tmpExpressionBed=$tmpPrefix.expression.bed
zcat $invNormTMMBed | head -n1 > $tmpExpressionBed
zgrep -m1 $geneID $invNormTMMBed >> $tmpExpressionBed


#=========================#
# Extract nominal results #
#=========================#

tmpNominalResults=$tmpPrefix.nomResults.txt

zcat $nomResultsFile | head -n 1 > $tmpNominalResults
zgrep $geneID $nomResultsFile >> $tmpNominalResults


#===========#
# Run SuSiE #
#===========#

$susieScript --geneID $geneID \
			 --gtPrefix $tmpVariantsPrefix \
			 --expressionBed $tmpExpressionBed \
			 --nominalResults $tmpNominalResults \
			 --covFile $covFile \
			 --outDir $outDir


#==========#
# Clean up #
#==========#

rm $tmpPrefix*
