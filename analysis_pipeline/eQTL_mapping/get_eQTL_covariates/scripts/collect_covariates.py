#!/usr/bin/env python3

import pandas as pd
import argparse
import numpy as np

#===============================#
# Collect and format covariates #
#===============================#

def main():

	parser = argparse.ArgumentParser(description="Collect and format covariates for use with FastQTL")
	parser.add_argument('-e', '--eigenvecs', required=True, help='Genotype PC ".eigenvec" file from Plink')
	parser.add_argument('-p', '--peerX', required=True, help='Output "X.csv" file from PEER')
	parser.add_argument('-a', '--addCov', required=False, help="Tab-separated file with additional covariates. Should have a header row. First column should be sample IDs, additional columns are additional covariates to include.")
	parser.add_argument('-o', '--out', required=True, help="Name of file to write collected formatted covariates to.")

	args = parser.parse_args()

	# Read in Genotype PCs
	fullDF = pd.read_csv(args.eigenvecs, delim_whitespace=True, header=None).T.drop(index=1)
	fullDF.columns = fullDF.iloc[0]
	fullDF.drop(fullDF.index[0], inplace=True)
	fullDF.columns.name = None
	fullDF.insert(0, 'id', [f'PC{i+1}' for i in range(fullDF.shape[0])])
	fullDF.set_index(keys='id', drop=True, inplace=True)

	# Add PEER covariates
	peerArray = np.loadtxt(args.peerX, delimiter=',', dtype=float).T
	peer_labels = [f'PEER{i+1}' for i in range(peerArray.shape[0])]

	for i, label in enumerate(peer_labels):
		fullDF.loc[label, :] = peerArray[i, :]

	# Add additional covariates
	addCovDF = pd.read_csv(args.addCov, header=0, sep='\t').T
	addCovDF.columns = addCovDF.iloc[0]
	addCovDF.drop(addCovDF.index[0], inplace=True)
	addCovDF.columns.name = None
	addCovDF.index.name = 'id'
	fullDF = pd.concat([fullDF, addCovDF], axis=0)

	# Output
	fullDF.to_csv(args.out, sep='\t', header=True, index=True)


if __name__=='__main__':
	main()
