#!/usr/bin/env python3

import argparse
import gzip
import sys


#==================================#
# Convert haploid calls to diploid #
#==================================#

def is_gzipped(path):
	with open(path, "rb") as f:
		return f.read(2) == b'\x1f\x8b'

def main():

	# Get arguments 
	parser = argparse.ArgumentParser(description='Add variant IDs to VCF file with missing IDs. These are in the format of the 1KGP NYGC VCF')
	parser.add_argument('inVCF', help='VCF with (at least some) missing variant IDs')
	args = parser.parse_args()

	# Convert genotypes
	open_fn = gzip.open if is_gzipped(args.inVCF) else open

	hap_to_dip = {'0':'0|0', '1':'1|1', '.':'.|.'}

	seen_pos = set()

	with open_fn(args.inVCF, 'rt') as in_fs:
		for i, line in enumerate(in_fs):
			if line.startswith('#'):
				sys.stdout.write(line)
				continue
			row = line.rstrip('\n').split('\t')
			chrom = row[0][3:]
			pos = int(row[1])
			snp_pos = f'{chrom}:{pos}'
			snp_id = row[2]
			ref = row[3]
			alt = row[4]

			if snp_id == '.':
				while True:
					if snp_pos in seen_pos:
						pos += 1
						snp_pos = f'{chrom}:{pos}'
					else:
						seen_pos.add(snp_pos)
						row[2] = f'{chrom}:{pos}:{ref[:10]}:{alt[:10]}'
						break

			out_row = '\t'.join(row) + '\n'
			sys.stdout.write(out_row)


if __name__ == '__main__':
	main()
