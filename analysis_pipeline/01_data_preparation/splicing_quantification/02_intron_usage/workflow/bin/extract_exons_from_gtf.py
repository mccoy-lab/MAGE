#!/usr/bin/env python3

import sys
import gzip
import argparse

def is_gzipped(path):
	with open(path, "rb") as f:
		return f.read(2) == b'\x1f\x8b'

if __name__=='__main__':

	parser = argparse.ArgumentParser(description="Extract exons from an input GTF file for use with Leafcutter's \"map_clusters_to_genes.R\" tool.")
	parser.add_argument('inGTF', help='The input GTF file.')
	parser.add_argument('-o', '--outputFile', help='The output file to write the exon descriptions to. If not included, will write to stdout.')
	args = parser.parse_args()

	gtfFile = args.inGTF
	open_fn = gzip.open if is_gzipped(gtfFile) else open
	out_fs = sys.stdout if not args.outputFile else open(args.outputFile, 'w')

	with open_fn(gtfFile, 'rt') as in_fs:
		out_fs.write('chr\tstart\tend\tstrand\tgene_id\tgene_name\n')

		seen_exons = set()
		for line in in_fs:
			if line.startswith('#'):
				continue
			fields = line.strip('\n').split('\t')
			featureType = fields[2]
			if featureType != 'exon':
				continue
			chrom = fields[0]
			start = fields[3]
			end = fields[4]
			strand = fields[6]
			info = {x.split(' ')[0]: x.split(' ')[1].strip('"') for x in fields[8].split('; ')}
			gene_id = info['gene_id']
			gene_name = info['gene_name']
			exon = (chrom, start, end, strand, gene_id, gene_name)
			if exon in seen_exons:
				continue
			seen_exons.add(exon)
			out_fs.write(f'{chrom}\t{start}\t{end}\t{strand}\t{gene_id}\t{gene_name}\n')
