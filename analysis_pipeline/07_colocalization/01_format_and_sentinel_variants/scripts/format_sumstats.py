#!/usr/bin/env python3

import argparse
import pandas as pd


#=================#
# Format sumstats #
#=================#

def main():

	parser = argparse.ArgumentParser(description="Format GWAS Catalog summary statistics")
	parser.add_argument('inFile', help='Raw GWAS Catalog summary statistics file')
	parser.add_argument('outFile', help='Output formatted GWAS catalog summary statistcs file')

	args = parser.parse_args()

	# Read in summary stats
	usecols = ['hm_chrom', 'hm_pos', 'hm_other_allele', 'hm_effect_allele', 'hm_rsid', 'hm_effect_allele_frequency', 'hm_beta', 'standard_error', 'p_value']
	sumstatsDF = pd.read_csv(args.inFile, header=0, sep='\t', low_memory=False, usecols=usecols+['hm_variant_id'])

	# Drop nan rows
	sumstatsDF = sumstatsDF.loc[~sumstatsDF['hm_variant_id'].isna(), :]

	# Drop duplicate rows
	sumstatsDF.drop_duplicates(subset='hm_variant_id', keep=False, inplace=True)

	# Rename chromosomes
	sumstatsDF['hm_chrom'] = sumstatsDF['hm_chrom'].apply(lambda i: f'chr{i}')

	# Rename and order columns
	col_map = {'hm_chrom': 'chr',
			   'hm_pos': 'pos',
			   'hm_other_allele': 'ref',
			   'hm_effect_allele': 'alt',
			   'hm_rsid': 'source_rsID',
			   'hm_effect_allele_frequency': 'af',
			   'hm_beta': 'beta',
			   'standard_error': 'se',
			   'p_value': 'pval'}
	ordered_cols = [col_map[col] for col in usecols]
	sumstatsDF = sumstatsDF.rename(columns=col_map).loc[:, ordered_cols]
	sumstatsDF['pos'] = sumstatsDF['pos'].astype(pd.Int64Dtype())

	# Limit to autosomes + chrX
	chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX']
	sumstatsDF = sumstatsDF.loc[sumstatsDF['chr'].isin(chroms)]

	# Order output
	def chromNum(chrom):
		if chrom == 'chrX':
			return 23
		return int(chrom[3:])
	sumstatsDF['chrNum'] = sumstatsDF['chr'].apply(chromNum)
	sumstatsDF.sort_values(by=['chrNum', 'pos', 'ref', 'alt'], inplace=True)
	sumstatsDF.drop(columns='chrNum', inplace=True)

	# Write output
	sumstatsDF.to_csv(args.outFile, header=True, index=False, sep='\t')


if __name__ == "__main__":
	main()
