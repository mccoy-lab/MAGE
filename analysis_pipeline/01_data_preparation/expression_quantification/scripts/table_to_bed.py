#!/usr/bin/env python3

import qtl.io
import qtl.norm
import pandas as pd
import argparse

#==========================#
# Convert table to TSS Bed #
#==========================#

def main():

	parser = argparse.ArgumentParser(description='Generate expression TSS BED files')
	parser.add_argument('-i', '--inputTable', required=True, help='Table of expression values with genes as rows, and samples as columns')
	parser.add_argument('-g', '--gtf', required=True, help='Annotation GTF file with transcript annotations. Transcripts should be labelled "transcript_id", and genes should be labeled "gene_id" in the attribute field')
	parser.add_argument('-o', '--out', required=True, help='Output TSS bed file')
	args = parser.parse_args()

	# Read-in expression values
	expression_df = pd.read_csv(args.inputTable, sep='\t', header=0, index_col=0)
		
	bed_template_df = qtl.io.gtf_to_tss_bed(args.gtf, feature='gene')
	bed_df = pd.merge(bed_template_df, expression_df, left_index=True, right_index=True)
	bed_df = bed_df.groupby('chr', sort=False, group_keys=False).apply(lambda x: x.sort_values('start'))
	bed_df.rename(columns={'chr': '#Chr', 'gene_id': 'ID'}, inplace=True)
	bed_df.to_csv(args.out, header=True, index=False, sep='\t')


if __name__=='__main__':
	main()
