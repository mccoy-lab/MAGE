#!/usr/bin/bash

#==========#
# Get args #
#==========#

leadHitSummaryFile=$1 # Lead hit output from SuSiE
variantInfoFile=$2 # Tab-separated file with position info for SuSiE lead eQTLs.
				   # Must have header row. First column must be the variantID
				   # as it appears in the SuSiE output. Must also contain at least
				   # a "CHROM" column and a (1-based) "POS" column
inVCF=$3 # Same VCF file used for eQTL mapping
expressionBed=$4 # TSS BED with DESeq2-normalized logged pseudocounts
covFile=$5 # Covariate file (same one used for FastQTL and SuSiE)
out_eQTLs=$6 # File of eQTLs, used as input for aFC-n, and the interaction test
out_aFCnFile=$7 # File to write aFCn results to

threads=8


#=============#
# Set scripts #
#=============#

# Get the directory of this script
thisDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

scriptDir=$thisDir/scripts

getQTLsScript=$scriptDir/get_qtls.py
correctCovariatesScript=$scriptDir/correct_covariates.ci95_wGT.py


#===================#
# Format eQTL input #
#===================#

$getQTLsScript $leadHitSummaryFile $variantInfoFile $out_eQTLs


#=================================================#
# Regress covariates out from logged pseudocounts #
#=================================================#

tmpPrefix=$(mktemp)
tmp_correctedExpression=$tmpPrefix.expression.covariate-corrected.ci95_wGT.bed

$correctCovariatesScript \
--vcf $inVCF \
--pheno $expressionBed \
--qtl $out_eQTLs \
--cov $covFile \
--out $tmp_correctedExpression

bgzip $tmp_correctedExpression
tabix -p bed $tmp_correctedExpression.gz


#=========================#
# Format expression input #
#=========================#

tmpExpressionFile=$tmpPrefix.expression.txt

zcat $tmp_correctedExpression | cut -f 4- | awk -v FS="\t" -v OFS="," 'NR==1 {$1="Name"} {$1=$1; print}' | gzip -c > $tmpExpressionFile


#===========#
# Run aFC-n #
#===========#

afcn.py --vcf $inVCF \
		--expr $tmpExpressionFile \
		--eqtl $out_eQTLs \
		--conf \
		--output $out_aFCnFile \
		--nthreads $threads


#==========#
# Clean up #
#==========#

rm $tmpPrefix*
