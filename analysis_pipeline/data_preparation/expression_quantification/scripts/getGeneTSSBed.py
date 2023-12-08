#!/usr/bin/env python3

import qtl.io
import argparse

#========================#
# Convert GTF to TSS Bed #
#========================#

def main():

	parser = argparse.ArgumentParser(description='Convert GTF to TSS Bed')
	parser.add_argument('-g', '--gtf', required=True, help='Annotation GTF file with transcript annotations. Transcripts should be labelled "transcript_id", and genes should be labeled "gene_id" in the attribute field')
	parser.add_argument('-o', '--out', required=True, help='Output TSS Bed')
	args = parser.parse_args()

	tss_bed_df = qtl.io.gtf_to_tss_bed(args.gtf, feature='gene')
	tss_bed_df.to_csv(args.out, header=False, index=False, sep='\t')

if __name__=='__main__':
	main()
