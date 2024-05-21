#!/usr/bin/env python3

import pandas as pd
import argparse
import gzip

#==================#
# Define functions #
#==================#

def openfile(filename, mode='rt'):
	if filename.endswith('.gz'):
		return gzip.open(filename, mode) 
	else:
		return open(filename, mode)


#=========================#
# Get significant results #
#=========================#

def main():

	parser = argparse.ArgumentParser(description='Using q-values from permutation pass, identify eGenes and subset significant nominal results')
	parser.add_argument('-n', '--nominalResults', required=True, help='Nominal pass eQTL results from FastQTL')
	parser.add_argument('-p', '--permResults', required=True, help='Permutation pass eQTL results from FastQTL')
	parser.add_argument('-e', '--out_eGene', required=True, help='File to write list of eGenes to')
	parser.add_argument('-s', '--out_sigResults', required=True, help='File to write significant nominal pass eQTL results to (will not be gzipped)')
	parser.add_argument('-f', '--FDR', required=False, default=0.05, type=float, help='False discovery rate to select eGenes (default 0.05)')

	args = parser.parse_args()

	# Read in permutation results
	permResultsDF = pd.read_csv(args.permResults, header=0, index_col=0, sep='\t')

	# Get eGenes
	eGenes = permResultsDF.loc[permResultsDF['qval']<=args.FDR].index.tolist()
	with openfile(args.out_eGene, 'wt') as out_fs:
		out_fs.write('\n'.join(eGenes) + '\n')

	# Get significant nominal associations
	eGeneThresholds = {eGene: permResultsDF.loc[eGene, 'pval_nominal_threshold'] for eGene in eGenes}
	
	with openfile(args.out_sigResults, 'wt') as out_fs, openfile(args.nominalResults, 'rt') as in_fs:
		for i, line in enumerate(in_fs):
			fields = line.rstrip('\n').split('\t')
			if i == 0:
				geneID_col = fields.index('gene_id')
				pval_col = fields.index('pval_nominal')
				out_fs.write(line)
				continue
			fields = line.rstrip('\n').split('\t')
			geneID = fields[geneID_col]
			pval = float(fields[pval_col])
			if geneID in eGeneThresholds and pval < eGeneThresholds[geneID]:
				out_fs.write(line)


if __name__ == '__main__':
	main()
