#!/usr/bin/bash

#==========#
# Get args #
#==========#

intronOutFileList=$1 # Four-column tab-delimted file with intronID, geneID, path to "cs_results.txt" file, path to "_full_pip_results.txt"
outCSSummary=$2 # sQTL_finemapping.introns.credibleSet_summary.MAGE.v0.1.txt
outFullSNPSummary=$3 # sQTL_finemapping.introns.allAssociations.MAGE.v0.1.txt
outSigSNPSummary=$4 # sQTL_finemapping.introns.significantAssociations.MAGE.v0.1.txt
outGeneSigSNPSummary=$5 # sQTL_finemapping.geneMerged.significantAssociations.MAGE.v0.1.txt
outGeneLeadHitSummary=$5 # sQTL_finemapping.geneMerged.leadVariant.MAGE.v0.1.txt


#=============#
# Set scripts #
#=============#

# Get the directory of this script
thisDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

scriptDir=$thisDir/scripts

collateCSResultsScript=$scriptDir/collate_cs_results.py
collateSNPResultsScript=$scriptDir/collate_snp_results.py
getGeneLevelResultsScript=$scriptDir/get_gene_level_results.py


#==================================#
# Get summary of all credible sets #
#==================================#

$collateCSResultsScript $intronOutFileList $outCSSummary


#=======================#
# Get SNP-level summary #
#=======================#

# Because the full SuSiE results are likely to be huge, it's better not to load it
# all into memory at once. As such, output files will need to be gzipped after
# finished writing
full_gzip_check=False
sig_gzip_check=False

if [[ $outFullSNPSummary == *.gz ]]; then
	outFullSNPSummary=${outFullSNPSummary%".gz"}
	full_gzip_check=True
fi

if [[ $outSigSNPSummary == *.gz ]]; then
	outSigSNPSummary=${outSigSNPSummary%".gz"}
	sig_gzip_check=True
fi

# Actually collect the results
$collateSNPResultsScript $intronOutFileList $outFullSNPSummary $outSigSNPSummary

if [[ $full_gzip_check == "True" ]]; then
	gzip $outFullSNPSummary
fi

if [[ $sig_gzip_check == "True" ]]; then
	gzip $outSigSNPSummary
fi


#========================#
# Get gene level results #
#========================#

$getGeneLevelResultsScript $outSigSNPSummary $outCSSummary $outGeneSigSNPSummary $outGeneLeadHitSummary
