#!/usr/bin/env python3

import gzip
import argparse
import pandas as pd

#==================#
# Define functions #
#==================#

def openfile(filename, mode='rt'):
	if filename.endswith('.gz'):
		return gzip.open(filename, mode) 
	else:
		return open(filename, mode)


#===================#
# Add location info #
#===================#

def main():

	parser =  argparse.ArgumentParser(description='Add intron gene coordinates (in bed format) to nominal results file from FastQTL, to allow tabix indexing')
	parser.add_argument('inFile',  help='Nominal results file from FastQTL')
	parser.add_argument('geneBed', help='Bed file with gene coordinates')
	parser.add_argument('outFile', help='Output file with bed format intron gene coordinates')
	
	args = parser.parse_args()

	geneCoords = {}
	for i, line in enumerate(openfile(args.geneBed, 'rt')):
		if i == 0:
			continue
		fields = line.rstrip('\n').split('\t')
		geneCoords[fields[3]] = tuple(fields[:3])

	with openfile(args.outFile, 'w') as out_fs, openfile(args.inFile, 'rt') as in_fs:
		for i, line in enumerate(in_fs):
			if i == 0:
				out_fs.write(f'geneChrom\tgeneStart\tgeneEnd\t{line}')
				continue
			intronID = line.split('\t')[0]
			chrom, chromStart, chromEnd = geneCoords[intronID]
			out_fs.write(f'{chrom}\t{chromStart}\t{chromEnd}\t{line}')


if __name__ == '__main__':
	main()
