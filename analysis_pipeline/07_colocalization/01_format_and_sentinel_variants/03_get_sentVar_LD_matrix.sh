#!/usr/bin/bash

#==========#
# Get args #
#==========#

sentinelVariant=$1 # Should be formatted chrom:pos:ref:alt (e.g. chr19:16092494:C:A). Position is expected to be 1-based
sumstatsFile=$2 # Formatted summary statistics for the corresponding trait
inVCF=$3 # VCF to be generating LD matrix from
selectedSampleList=$4 # File with samples to be subset from VCF (if not all samples were used for each trait GWAS)
outPrefix=$5 # Prefix (including path) of output LD matrix and SNP list files


# Global variables
colocWindow=500000 # +/- 500kb around each sentinel variant will be tested for colocalization


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

printf "\nGetting LD matrix in the region $chrom:$colocStart-$colocEnd\n"


#===================#
# Set up temp files #
#===================#

tmpPrefix=$(mktemp)


#============#
# Susbet VCF #
#============#

tmpVCF=$tmpPrefix.subset.vcf

printf "\t- Susbetting VCF to selected samples...\n"
bcftools view --force-samples --regions-overlap=pos -S $selectedSampleList -o $tmpVCF $inVCF $chrom:$colocStart-$colocEnd > /dev/null 2>&1

nSamplesOut=$(grep -m1 "#CHROM" $tmpVCF | cut -f 10- | tr -s '\t' '\n' | wc -l)
printf "\t\t* $nSamplesOut/$nSamplesIn samples found in VCF\n"


#=======================#
# Get sumstats variants #
#=======================#

tmpSnpListFile=$tmpPrefix.sumstats.snps

printf "\t- Getting list of variants in this region with genotypes AND sumstats...\n"
tabix $sumstatsFile $chrom:$colocStart-$colocEnd | awk '{print $1":"$2":"$3":"$4}'> $tmpSnpListFile


#===============#
# Get LD matrix #
#===============#

printf "\t- Getting LD matrix for these variants...\n"
plink --vcf $tmpVCF --keep-allele-order --extract $tmpSnpListFile --r square --make-bed --out $outPrefix > /dev/null 2>&1


#====================#
# Get final SNP list #
#====================#

outSnpList=$outPrefix.snps

printf "chr\tpos\tref\talt\tvariantID\n" > $outSnpList
cat $tmpPrefix.bim | awk -v OFS='\t' '{print "chr"$1, $4, $6, $5, $2}' >> $outSnpList


#==========#
# Clean-up #
#==========#

rm $tmpPrefix*
printf "All done.\n"
