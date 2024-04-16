#!/usr/bin/bash


#==========#
# Get args #
#==========#

sentinelVariant=$1 # Should be formatted chrom:pos:ref:alt (e.g. chr19:16092494:C:A). Position is expected to be 1-based
sumstatsFile=$2 # Formatted summary statistics for the corresponding trait
ldPrefix=$3 # Prefix (including path) of GWAS LD matrix (and snp list) for this sentinel variant
numSamples=$4 # Number of samples included in the GWAS
traitType=$5 # GWAS trait type. "quant" for quantitative traits, "cc" for case/control traits

sQTL_VCF=$6 # VCF file used for sQTL mapping
sQTL_permutationResultsFile=$7 # Tabixed sQTL results from FastQTL permutation pass 
sQTL_nominalResultsFile=$8 # Tabixed sQTL results from FastQTL nominal pass
sQTL_covFile=$9 # Covariates file used for sQTL mapping

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
# Get sGenes #
#============#

tmpPermResults=$tmpPrefix.$sentinelVariant.permutationResults.txt

printf "\t* Looking for sGenes in region $chrom:$colocStart-$colocEnd... "

zcat $sQTL_permutationResultsFile | head -n 1 > $tmpPermResults
tabix $sQTL_permutationResultsFile $chrom:$colocStart-$colocEnd | awk '$25 <= 0.05' >> $tmpPermResults # sGenes at 5% FDR

n_sGenes=$(cat $tmpPermResults | tail -n +2 | wc -l)
printf "$n_sGenes found.\n"

if [[ $n_sGenes == 0 ]]; then
	rm $tmpPrefix*
	printf "All done.\n"
	exit 0
fi

grepString=$(cat $tmpPermResults | tail -n +2 | cut -f 3 | tr '\n' '|' | sed '$ s/.$//') # list of sGenes used for subsetting nominal results
rm $tmpPermResults


#=====================#
# Get nominal results #
#=====================#

tmpNominalResults=$tmpPrefix.$sentinelVariant.nominalResults.txt

printf "\t* Grabbing sQTL nominal results for these genes... "

zcat $sQTL_nominalResultsFile | head -n1 > $tmpNominalResults
tabix $sQTL_nominalResultsFile $chrom:$colocStart-$colocEnd | grep -E "$grepString" >> $tmpNominalResults

printf "Done.\n"


#========================#
# Get sQTL genotype info #
#========================#

tmpVariantsPrefix=$tmpPrefix.$sentinelVariant.variants

printf "\t* Getting formatted genotypes for sQTLs... "

bcftools view -r $chrom:$colocStart-$colocEnd --min-af=0.01:minor --regions-overlap=pos $sQTL_VCF > $tmpVariantsPrefix.vcf
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
			 --sQTL_nominalResultsFile $tmpNominalResults \
			 --sQTL_gtPrefix $tmpVariantsPrefix \
			 --sQTL_covFile $sQTL_covFile \
			 --outPrefix $outPrefix \
			 --overlap_vars

printf "Done.\n"


#==========#
# Clean-up #
#==========#

rm $tmpPrefix*
printf "All done.\n"
