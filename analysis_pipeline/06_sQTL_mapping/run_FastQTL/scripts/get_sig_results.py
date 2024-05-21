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

	parser = argparse.ArgumentParser(description='Using q-values from permutation pass, identify sGenes and subset significant nominal results')
	parser.add_argument('-n', '--nominalResults', required=True, help='Nominal pass sQTL results from FastQTL')
	parser.add_argument('-p', '--permResults', required=True, help='Permutation pass sQTL results from FastQTL')
	parser.add_argument('-g', '--phenGroups', required=True, help='File containing phenotype groups used for sQTL mapping')
	parser.add_argument('-e', '--out_sGene', required=True, help='File to write list of sGenes to')
	parser.add_argument('-s', '--out_sigResults', required=True, help='File to write significant nominal pass sQTL results to (will not be gzipped)')
	parser.add_argument('-f', '--FDR', required=False, default=0.05, type=float, help='False discovery rate to select sGenes (default 0.05)')

	args = parser.parse_args()

	# Read in permutation results
	permResultsDF = pd.read_csv(args.permResults, header=0, index_col=17, sep='\t')

	# Get sGenes
	sGenes = permResultsDF.loc[permResultsDF['qval']<=args.FDR].index.tolist()
	with openfile(args.out_sGene, 'wt') as out_fs:
		out_fs.write('\n'.join(sGenes) + '\n')
	
	# Get phenotype groups
	with open(args.phenGroups, 'rt') as in_fs:
		intronGroups = {line.rstrip('\n').split('\t')[0]: line.rstrip('\n').split('\t')[1] for line in in_fs}

	# Get significant nominal associations
	sGeneThresholds = {sGene: permResultsDF.loc[sGene, 'pval_nominal_threshold'] for sGene in sGenes}

	with openfile(args.out_sigResults, 'wt') as out_fs, openfile(args.nominalResults, 'rt') as in_fs:
		for i, line in enumerate(in_fs):
			fields = line.rstrip('\n').split('\t')
			if i == 0:
				intronID_col = fields.index('gene_id')
				pval_col = fields.index('pval_nominal')
				out_fs.write(line)
				continue
			fields = line.rstrip('\n').split('\t')
			intronID = fields[intronID_col]
			pval = float(fields[pval_col])
			if intronGroups[intronID] in sGeneThresholds and pval < sGeneThresholds[intronGroups[intronID]]:
				out_fs.write(line)


if __name__ == '__main__':
	main()
