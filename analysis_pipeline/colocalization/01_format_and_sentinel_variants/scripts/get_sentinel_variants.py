#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse


#============================#
# Get GWAS sentinel variants #
#============================#

def main():

	parser = argparse.ArgumentParser(description='Get sentinel variants from formatted GWAS sumstats')
	parser.add_argument('sumstatsFile', type=str, help='Formatted GWAS summary statistics file')
	parser.add_argument('windowSize', type=int, help='+/- window size around each sentinel variant that cannot contain another sentinel variant')
	parser.add_argument('min_pval', type=float, help='Minimum p-value to be considered for sentinel variants')
	parser.add_argument('outFile', type=str, help='File to write sentinel variants results to')

	args = parser.parse_args()

	# Read in summary statistics
	sumstats = pd.read_csv(args.sumstatsFile, header=0, sep='\t', low_memory=False)

	# Get significant variants
	sorted_sumstats = sumstats.loc[sumstats['pval'] <= args.min_pval].sort_values(by='pval', ascending=True)

	# Get sentinel variants
	sentinel_variants = set()

	while len(sorted_sumstats) > 0:

		top_hit = sorted_sumstats.index[0]
		top_hit_chrom, top_hit_pos = sorted_sumstats.loc[top_hit, ['chr', 'pos']]

		variants_to_drop = sorted_sumstats.loc[(sorted_sumstats['chr']==top_hit_chrom) & (sorted_sumstats['pos'] >= top_hit_pos-args.windowSize) & (sorted_sumstats['pos'] < top_hit_pos+args.windowSize), :].index.to_list()
		sorted_sumstats.drop(variants_to_drop, inplace=True)

		sentinel_variants.add(top_hit)

	# Subset data and write to output
	sumstatsSentinel = sumstats.loc[sumstats.index.isin(sentinel_variants), :]
	sumstatsSentinel.to_csv(args.outFile, header=True, index=False, sep='\t', na_rep='NA')


if __name__ == "__main__":
	main()
