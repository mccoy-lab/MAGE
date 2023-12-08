#!/usr/bin/env python3

import pandas as pd
import argparse

#====================#
# Collect CS results #
#====================#

def main():

	parser = argparse.ArgumentParser(description='Collate credible set-level sQTL finemapping results into a single file describing results for each credible set across all introns')
	parser.add_argument('intronResultsMap', help='Map of intronID to geneID and output credible set and full results files')
	parser.add_argument('outFile', help='Name output table file')
	args = parser.parse_args()

	csResultsDF = pd.DataFrame({'intronID': pd.Series(dtype=object),
								'geneID': pd.Series(dtype=object),
								'credible_set': pd.Series(dtype=object),
								'cs_nvars': pd.Series(dtype=object),
								'cs_coverage': pd.Series(dtype=object),
								'cs_min_corr': pd.Series(dtype=object),
								'cs_mean_corr': pd.Series(dtype=object),
								'cs_median_corr': pd.Series(dtype=object)})

	intronFilesDF = pd.read_csv(args.intronResultsMap, sep='\t', header=0, index_col=0)

	for intronID in intronFilesDF.index:
		geneID = intronFilesDF.loc[intronID, 'geneID']
		intron_csResultsFile = intronFilesDF.loc[intronID, 'csResultsFile']
		intron_csResultsDF = pd.read_csv(intron_csResultsFile, sep='\t', header=0)
		if len(intron_csResultsDF) == 0:
			continue
		intron_csResultsDF.insert(loc=0, column='geneID', value=geneID)
		intron_csResultsDF.insert(loc=0, column='intronID', value=intronID)
		csResultsDF = pd.concat([csResultsDF, intron_csResultsDF], axis=0)

	csResultsDF.rename(columns={'credible_set':'credibleSet',
								'cs_nvars':'numVariants',
								'cs_coverage':'coverage',
								'cs_min_corr':'minCorr',
								'cs_mean_corr':'meanCorr',
								'cs_median_corr':'medianCorr'
								}, inplace=True)

	csResultsDF.to_csv(args.outFile, sep='\t', header=True, index=False)


if __name__=='__main__':
	main()
