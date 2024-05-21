#!/usr/bin/env python3

import pandas as pd
import argparse

#============================#
# Format eQTL input for aFCn #
#============================#

def main():

	parser = argparse.ArgumentParser(description="Format SuSiE lead eQTLs to be used with aFC-n")
	parser.add_argument('susieLeadHitsFile', help='Lead eQTL output file from SuSiE')
	parser.add_argument('variantInfoFile', help='Tab-separated file with position info for SuSiE lead eQTLs. Must have header row. First column must be the variantID as it appears in the SuSiE output. Must also contain at least a "CHROM" column and a (1-based) "POS" column')
	parser.add_argument('outFile', help="File to write formatted eQTL data to")

	args = parser.parse_args()

	# Read in data
	susieLeadHits = pd.read_csv(args.susieLeadHitsFile, header=0, sep='\t')
	variantInfo = pd.read_csv(args.variantInfoFile, header=0, index_col=0, sep='\t')

	# Create output data
	outData = susieLeadHits[['geneID', 'variantID']].copy().rename(columns={'geneID': 'gene_id', 'variantID': 'variant_id'})
	outData['variant_chr'] = outData['variant_id'].map(variantInfo['CHROM'])
	outData['variant_pos'] = outData['variant_id'].map(variantInfo['POS'])

	# Write to output
	outData.to_csv(args.outFile, header=True, index=False, sep='\t')


if __name__ == "__main__":
	main()
