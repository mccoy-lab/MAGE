#!/usr/bin/env python3

import argparse
import pandas as pd

#================================#
# Select lead eQTLs from each CS #
#================================#

def main():

	parser = argparse.ArgumentParser(description='Get the lead association for each credible set in an input Susie results file')
	parser.add_argument('susieResults', help="Output of significant associations from SuSiE.")
	parser.add_argument('outFile', help="File to write lead hit results to")
	args = parser.parse_args()

	credibleSetResults = pd.read_csv(args.susieResults, sep='\t', header=0).sort_values(by=['geneID', 'variantCredibleSet', 'variantPIP'], ascending=[True, True, False]).reset_index(drop=True)
	leadHitResults = credibleSetResults.drop_duplicates(subset=['geneID', 'variantCredibleSet'], keep='first')

	leadHitResults.to_csv(args.outFile, sep='\t', header=True, index=False)


if __name__=='__main__':
	main()
