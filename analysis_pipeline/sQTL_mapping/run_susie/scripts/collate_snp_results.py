#!/usr/bin/env python3

import pandas as pd
import argparse

#===========================#
# Collect SNP-level results #
#===========================#

def main():

	parser = argparse.ArgumentParser(description='Collate snp-level sQTL finemapping results into a file describing all results, and a file describing only fine-mapped variants')
	parser.add_argument('intronResultsMap', help='Map of intronID to geneID and output credible set and full results files')
	parser.add_argument('outAllFile', help='Name output table file with ALL results (cannot be gzipped)')
	parser.add_argument('outSigFile', help='Name output table file with only significant associations after finemapping (cannot be gzipped)')
	args = parser.parse_args()

	intronFilesDF = pd.read_csv(args.intronResultsMap, sep='\t', header=0, index_col=0)

	i = 0
	for intronID in intronFilesDF.index:
		geneID = intronFilesDF.loc[intronID, 'geneID']
		intron_fullResultsFile = intronFilesDF.loc[intronID, 'fullResultsFile']
		intron_fullResultsDF = pd.read_csv(intron_fullResultsFile, sep='\t', header=0)
		intron_fullResultsDF.insert(loc=0, column='geneID', value=intronID)
		intron_fullResultsDF.insert(loc=0, column='intronID', value=intronID)
		intron_fullResultsDF.rename(columns={'variant_id':'variantID',
										   'pip':'variantPIP',
										   'credible_set':'variantCredibleSet'}, inplace=True)
		intron_sigResultsDF = intron_fullResultsDF[~intron_fullResultsDF['variantCredibleSet'].isna()]

		if i == 0:
			intron_fullResultsDF.to_csv(args.outAllFile, sep='\t', header=True, index=False, na_rep='NA')
			intron_sigResultsDF.to_csv(args.outSigFile, sep='\t', header=True, index=False, na_rep='NA')

		else:
			intron_fullResultsDF.to_csv(args.outAllFile, mode='a', sep='\t', header=False, index=False, na_rep='NA')
			intron_sigResultsDF.to_csv(args.outSigFile, mode='a', sep='\t', header=False, index=False, na_rep='NA')

		i+=1


if __name__=='__main__':
	main()
