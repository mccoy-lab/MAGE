#!/usr/bin/env python3

import qtl.io
import qtl.norm
import pandas as pd
import argparse

#==========================#
# Inverse normal transform #
#==========================#

def main():

	parser = argparse.ArgumentParser(description='Generate inverse normal transformed TMM values for eQTL analyses')
	parser.add_argument('-t', '--tmm', required=True, help='Table of expression TMM values with genes as rows, and samples as columns')
	parser.add_argument('-o', '--out', required=True, help='Name output table file')
	args = parser.parse_args()

	# Read-in tmm values
	tmm_df = pd.read_csv(args.tmm, sep='\t', header=0, index_col=0)

	# Inverse normalize TMM values
	norm_df = qtl.norm.inverse_normal_transform(tmm_df)

	# Write to output file
	norm_df.to_csv(args.out, header=True, index=True, sep='\t')


if __name__=='__main__':
	main()
