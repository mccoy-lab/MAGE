#!/usr/bin/bash

#==========#
# Get args #
#==========#

geneID=$1 # This script is meant to be run with a single gene. It can be parallelized using something like SLURM
inVCF=$2 # The same VCF used with FastQTL
normedRatiosBed=$3 # The same filtered, normalized splicing ratios bed file used with FastQTL
covFile=$4 # The same covariates file used with FastQTL
phenGroupsFile=$5
nomResultsFile=$6 # The nominal pass output from FastQTL, indexed with tabix
outDir=$7

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

tmpGenePrefix=$(mktemp)
tmpGeneBed=$tmpGenePrefix.gene.bed
tmpVariantsPrefix=$tmpGenePrefix.variants

firstIntron=`cat $phenGroupsFile | awk -v geneID=$geneID '$2==geneID {print $1}' | head -n 1`

zgrep -m1 $firstIntron $normedRatiosBed | awk -v window=$window 'BEGIN{OFS="\t";} {print $1, ($2-window>0)?$2-window:0, $3+window}' > $tmpGeneBed
bcftools view -R $tmpGeneBed --min-af ${mafThreshold}:minor $inVCF > $tmpVariantsPrefix.vcf
vcftools --gzvcf $tmpVariantsPrefix.vcf --012 --out $tmpVariantsPrefix 2> /dev/null
grep ^[^#] $tmpVariantsPrefix.vcf | cut -f 3 > $tmpVariantsPrefix.012.snps

# Clean-up unneeded files
rm $tmpGenePrefix
rm $tmpVariantsPrefix.vcf
rm $tmpVariantsPrefix.012.pos


#==========================================#
# Run SuSiE for each intron of target gene #
#==========================================#

for intronID in `cat $phenGroupsFile | awk -v geneID=$geneID '$2==geneID {print $1}'`; do

	tmpIntronPrefix=$(mktemp)

	echo $tmpIntronPrefix

	#=========================#
	# Extract usage of intron #
	#=========================#

	tmpUsageBed=$tmpIntronPrefix.usage.bed
	zcat $normedRatiosBed | head -n1 > $tmpUsageBed
	tabix -R $tmpGeneBed $normedRatiosBed | grep -m1 $intronID $normedRatiosBed >> $tmpUsageBed

	#=========================#
	# Extract nominal results #
	#=========================#

	tmpNominalResults=$tmpIntronPrefix.nomResults.txt

	zcat $nomResultsFile | head -n 1 | cut -f 4- > $tmpNominalResults
	tabix -R $tmpGeneBed $nomResultsFile | grep $intronID | cut -f 4- >> $tmpNominalResults

	#===========#
	# Run SuSiE #
	#===========#

	$susieScript --intronID $intronID \
				 --gtPrefix $tmpVariantsPrefix \
				 --expressionBed $tmpUsageBed \
				 --nominalResults $tmpNominalResults \
				 --covFile $covFile \
				 --outDir $outDir

	#=======================#
	# Clean up intron stuff #
	#=======================#

	rm $tmpIntronPrefix*

done


#=====================#
# Clean up gene stuff #
#=====================#

rm $tmpGenePrefix*
