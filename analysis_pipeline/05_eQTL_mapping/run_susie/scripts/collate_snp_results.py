#!/usr/bin/env python3

import pandas as pd
import argparse

#===========================#
# Collect SNP-level results #
#===========================#

def main():

	parser = argparse.ArgumentParser(description='Collate snp-level eQTL finemapping results into a file describing all results, and a file describing onlt fine-mapped variants')
	parser.add_argument('geneResultsMap', help='Map of geneID to output credible set and full results files')
	parser.add_argument('outAllFile', help='Name output table file with ALL results (cannot be gzipped)')
	parser.add_argument('outSigFile', help='Name output table file with only significant associations after finemapping (cannot be gzipped)')
	args = parser.parse_args()

	geneFilesDF = pd.read_csv(args.geneResultsMap, sep='\t', header=0, index_col=0)

	i = 0
	for geneID in geneFilesDF.index:
		gene_fullResultsFile = geneFilesDF.loc[geneID, 'fullResultsFile']
		gene_fullResultsDF = pd.read_csv(gene_fullResultsFile, sep='\t', header=0)
		gene_fullResultsDF.insert(loc=0, column='geneID', value=geneID)
		gene_fullResultsDF.rename(columns={'variant_id':'variantID',
										   'pip':'variantPIP',
										   'credible_set':'variantCredibleSet'}, inplace=True)
		gene_sigResultsDF = gene_fullResultsDF[~gene_fullResultsDF['variantCredibleSet'].isna()]

		if i == 0:
			gene_fullResultsDF.to_csv(args.outAllFile, sep='\t', header=True, index=False, na_rep='NA')
			gene_sigResultsDF.to_csv(args.outSigFile, sep='\t', header=True, index=False, na_rep='NA')

		else:
			gene_fullResultsDF.to_csv(args.outAllFile, mode='a', sep='\t', header=False, index=False, na_rep='NA')
			gene_sigResultsDF.to_csv(args.outSigFile, mode='a', sep='\t', header=False, index=False, na_rep='NA')

		i+=1


if __name__=='__main__':
	main()
