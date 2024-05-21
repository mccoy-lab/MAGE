#!/usr/bin/env python3

import gzip
import argparse

#==================#
# Define functions #
#==================#

def openfile(filename, mode='r'):
	if filename.endswith('.gz'):
		return gzip.open(filename, mode) 
	else:
		return open(filename, mode)


#==========================#
# Map transcripts to genes #
#==========================#

def main():

	parser = argparse.ArgumentParser(description='From an input transcript GTF file (e.g. "gencode.v38.annotation.gtf.gz"), generate an output transcript-to-gene mapping file')
	parser.add_argument('-g', '--gtf', required=True, help='Annotation GTF file with transcript annotations. Transcripts should be labelled "transcript_id", and genes should be labeled "gene_id" in the attribute field')
	parser.add_argument('-o', '--out', required=True, help='Output comma-separated file mapping transcripts to genes')
	args = parser.parse_args()

	with openfile(args.gtf, 'rt') as inFS, openfile(args.out, 'wt') as outFS:
		outFS.write('TXNAME,GENEID\n')
		for line in inFS:
			if line.startswith('#'):
				continue
			fields = line.strip('\n').split('\t')
			if fields[2] != 'transcript':
				continue
			txInfo = dict((x.split(' ')[0], x.split(' ')[1].strip('"')) if len(x.split(' '))==2 else (x, None) for x in fields[8].split('; '))
			outFS.write(f'{txInfo["transcript_id"]},{txInfo["gene_id"]}\n')

if __name__ == "__main__":
	main()
