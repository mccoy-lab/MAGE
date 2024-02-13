#!/usr/bin/bash

#==========#
# Get args #
#==========#

geneOutFileList=$1 # Three-column tab-delimted file with geneID, path to "cs_results.txt" file, path to "_full_pip_results.txt"
outCSSummary=$2 # eQTL_finemapping.credibleSet_summary.MAGE.v0.1.txt
outFullSNPSummary=$3 # eQTL_finemapping.allAssociations.MAGE.v0.1.txt
outSigSNPSummary=$4 # eQTL_finemapping.significantAssociations.MAGE.v0.1.txt
outLeadHitSummary=$5 # eQTL_finemapping.leadVariant_aFCn.MAGE.v0.1.txt


#=============#
# Set scripts #
#=============#

# Get the directory of this script
thisDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

scriptDir=$thisDir/scripts

collateCSResultsScript=$scriptDir/collate_cs_results.py
collateSNPResultsScript=$scriptDir/collate_snp_results.py
leadHitScript=$scriptDir/cs_lead_hit.py


#==================================#
# Get summary of all credible sets #
#==================================#

$collateCSResultsScript $geneOutFileList $outCSSummary


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
$collateSNPResultsScript $geneOutFileList $outFullSNPSummary $outSigSNPSummary

if [[ $full_gzip_check == "True" ]]; then
	gzip $outFullSNPSummary
fi

if [[ $sig_gzip_check == "True" ]]; then
	gzip $outSigSNPSummary
fi


#===================================================#
# Get best hits (highest PIP) for each credible set #
#===================================================#

$leadHitScript $outSigSNPSummary $outLeadHitSummary
