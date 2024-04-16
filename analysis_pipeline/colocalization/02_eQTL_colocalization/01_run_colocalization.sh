#!/usr/bin/bash


#==========#
# Get args #
#==========#

sentinelVariant=$1 # Should be formatted chrom:pos:ref:alt (e.g. chr19:16092494:C:A). Position is expected to be 1-based
sumstatsFile=$2 # Formatted summary statistics for the corresponding trait
ldPrefix=$3 # Prefix (including path) of GWAS LD matrix (and snp list) for this sentinel variant
numSamples=$4 # Number of samples included in the GWAS
traitType=$5 # GWAS trait type. "quant" for quantitative traits, "cc" for case/control traits

eQTL_VCF=$6 # VCF file used for eQTL mapping
eQTL_permutationResultsFile=$7 # Tabixed eQTL results from FastQTL permutation pass 
eQTL_nominalResultsFile=$8 # Tabixed eQTL results from FastQTL nominal pass
eQTL_covFile=$9 # Covariates file used for eQTL mapping

outPrefix=${10} # Prefix (including path of output files)

# Global variables
colocWindow=500000 # +/- 500kb around each sentinel variant will be tested for colocalization


#=============#
# Set scripts #
#=============#

# Get the directory of this script
thisDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
scriptDir=$thisDir/scripts

colocScript=$scriptDir/run_coloc.R


#===========================#
# Get sentinel variant info #
#===========================#

chrom=$(echo $sentinelVariant | cut -d':' -f1)
pos=$(echo $sentinelVariant | cut -d':' -f2)
ref=$(echo $sentinelVariant | cut -d':' -f3)
alt=$(echo $sentinelVariant | cut -d':' -f4)


#================================#
# Get colocalization region info #
#================================#

colocStart=$(( $pos-$colocWindow ))
colocStart=$(( $colocStart > 1 ? $colocStart : 1 )) # Nothing smaller than 1
colocEnd=$(( $pos+$colocWindow-1 ))

printf "\nRunning colocalization in region $chrom:$colocStart-$colocEnd\n"


#======================#
# Check if files exist #
#======================#

if [[ ! -f $ldPrefix.ld ]]; then
	printf "\tLD matrix file does not exist. Stopping\n"
	printf "\nAll done.\n"
	exit 0
fi


#===================#
# Set up temp files #
#===================#

tmpPrefix=$(mktemp)


#=============================#
# Get GWAS sumstats in region #
#=============================#

tmpSumstatsFile=$tmpPrefix.$sentinelVariant.sumstats.tsv

printf "\t* Getting GWAS sumstats in this region...\n"
zcat $sumstatsFile | head -n1 > $tmpSumstatsFile
tabix $sumstatsFile $chrom:$colocStart-$colocEnd >> $tmpSumstatsFile


#============#
# Get eGenes #
#============#

tmpPermResults=$tmpPrefix.$sentinelVariant.permutationResults.txt

printf "\t* Looking for eGenes in region $chrom:$colocStart-$colocEnd... "

zcat $eQTL_permutationResultsFile | head -n 1 > $tmpPermResults
tabix $eQTL_permutationResultsFile $chrom:$colocStart-$colocEnd | awk '$22 <= 0.05' >> $tmpPermResults # eGenes at 5% FDR

n_eGenes=$(cat $tmpPermResults | tail -n +2 | wc -l)
printf "$n_eGenes found.\n"

if [[ $n_eGenes == 0 ]]; then
	rm $tmpPrefix*
	printf "All done.\n"
	exit 0
fi

grepString=$(cat $tmpPermResults | tail -n +2 | cut -f 3 | tr '\n' '|' | sed '$ s/.$//') # list of eGenes used for subsetting nominal results
rm $tmpPermResults


#=====================#
# Get nominal results #
#=====================#

tmpNominalResults=$tmpPrefix.$sentinelVariant.nominalResults.txt

printf "\t* Grabbing eQTL nominal results for these genes... "

zcat $eQTL_nominalResultsFile | head -n1 > $tmpNominalResults
tabix $eQTL_nominalResultsFile $chrom:$colocStart-$colocEnd | grep -E "$grepString" >> $tmpNominalResults

printf "Done.\n"


#========================#
# Get eQTL genotype info #
#========================#

tmpVariantsPrefix=$tmpPrefix.$sentinelVariant.variants

printf "\t* Getting formatted genotypes for eQTLs... "

bcftools view -r $chrom:$colocStart-$colocEnd --min-af=0.01:minor --regions-overlap=pos $eQTL_VCF > $tmpVariantsPrefix.vcf
vcftools --gzvcf $tmpVariantsPrefix.vcf --012 --out $tmpVariantsPrefix 2> /dev/null
grep ^[^#] $tmpVariantsPrefix.vcf | awk '{ print $1"_"$2"_"$4"_"$5 }' > $tmpVariantsPrefix.012.snps

printf "Done.\n"

rm $tmpVariantsPrefix.vcf
rm $tmpVariantsPrefix.012.pos


#====================#
# Run colocalization #
#====================#

printf "\t* Running coloc-susie in this region...\n"

$colocScript --gwas_sumstatsFile $tmpSumstatsFile \
			 --gwas_ldPrefix $ldPrefix \
			 --gwas_sampleSize $numSamples \
			 --gwas_type $traitType \
			 --eQTL_nominalResultsFile $tmpNominalResults \
			 --eQTL_gtPrefix $tmpVariantsPrefix \
			 --eQTL_covFile $eQTL_covFile \
			 --outPrefix $outPrefix \
			 --overlap_vars

printf "Done.\n"


#==========#
# Clean-up #
#==========#

rm $tmpPrefix*
printf "All done.\n"
