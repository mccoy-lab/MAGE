#!/usr/bin/env python3

from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
import pandas as pd
import numpy as np
import argparse
import pysam
import gzip
import warnings
import sys
warnings.filterwarnings("ignore", category=DeprecationWarning) 


#==================#
# Define functions #
#==================#

# Parsing functions are from my fork of the afcn tool: https://github.com/dtaylo95/aFCn

def read_eqtls(eqtl_filename):
	
	'''Used to read in the eqtl file into a pandas dataframe'''

	#open the eqtl file
	eqtl_file = pd.read_csv(eqtl_filename, sep='\t')
	
	# #drop gene versions
	# eqtl_file['gene_id_clean'] = eqtl_file.gene_id.str.split('.').str[0]
	
	# Clean up gene and variant IDs
	gene_map = {gene_id: f'g{i}' for i, gene_id in enumerate(eqtl_file.gene_id.unique())}
	eqtl_file['gene_id_clean'] = eqtl_file.gene_id.map(gene_map)

	variant_map = {variant_id: f'v{i}' for i, variant_id in enumerate(eqtl_file.variant_id.unique())}
	eqtl_file['variant_id_clean'] = eqtl_file.variant_id.map(variant_map)

	# if is_variant_id_format_valid(eqtl_file) == False:

	#     raise Exception('''Variant IDs must not begin with a digit, 
	#             reformat your vcf and EQTL matrix''')
 

	return eqtl_file

def read_expressions(expressions_filename, eqtl_file):

	'''Used to read the relevant expressions into a pandas dataframe'''
	
	gene_map = eqtl_file.loc[:, ['gene_id', 'gene_id_clean']].set_index('gene_id')['gene_id_clean'].to_dict()
	expr_df = []
	expr_df_cols = []
	#open the expression file, get the expressions for the genes we have eqtls for
	with gzip.open(expressions_filename) as f:

		for line in f:

			line = line.decode()
			line = line.replace('\n','')

			expr_data = line.split(',')

			if expr_data[0] == 'Name':

				expr_df_cols = expr_data


			else:

				#try splitting the gene by version number, in 
				#case our user forgot
				thisgene = expr_data[0]
				
				if thisgene in gene_map:
					expr_df.append(expr_data)

	expr_dataframe = pd.DataFrame(expr_df, columns=expr_df_cols)

	#if no genes were found
	if expr_dataframe.empty:
		raise Exception("No matching genes found between the EQTL and expressions files")

	#now do the gene version correction once and for all if needed
	expr_dataframe['Name'] = expr_dataframe['Name'].map(gene_map)
	
	return expr_dataframe

def get_vcf_header(vcfName):

	''' This function returns the header of a gzipped or non-gzipped vcf file '''

	#open the gzipped file 
	vcf_stream = gzip.open(vcfName)


	#find the first line that does not begin with "#"
	for line in vcf_stream:

		if isinstance(line, bytes) and not isinstance(line, str):

			line = line.decode()

		if '#CHROM' in line:

			return line.replace('\n','').split('\t')

def read_haplotypes(vcfName, eqtl_file):

	'''Reads in the variants we have eqtls and expressions(genes) for'''
	
	variant_map = eqtl_file.loc[:, ['variant_id', 'variant_id_clean']].set_index('variant_id')['variant_id_clean'].to_dict()

	#match these to the vcf file
	vcf_df = []
	vcf_df_cols = []
	#the eqtl variants
	counter = 0
	tabix_haplotypes = pysam.Tabixfile( vcfName,"r")


	#write the header first
	vcf_df_cols = get_vcf_header(vcfName)

	for variant in variant_map:

		#get coordinates
		chrom = eqtl_file.loc[eqtl_file['variant_id']==variant, ['variant_chr']].iloc[0,0]
		pos = eqtl_file.loc[eqtl_file['variant_id']==variant, ['variant_pos']].iloc[0,0]

		#query vcf using tabix
		try:

			records = tabix_haplotypes.fetch(chrom, int(pos) - 1, int(pos))

			for record in records:
				record_info = record.split('\t')
				if record_info[2] == variant:
					vcf_df.append(record_info)


		except:

			#variant not found
			continue

	vcf_dataframe = pd.DataFrame(vcf_df, columns=vcf_df_cols).drop_duplicates()

	# Use clean variant ids
	vcf_dataframe['ID'] = vcf_dataframe.ID.map(variant_map)

	#do a quick check and raise an exception if no eqtls were found
	if vcf_dataframe.empty:

		raise Exception("No matching EQTLs were found in VCF file - check variant ID formatting")

	
	return vcf_dataframe

def read_pop_table(tableName):
	pop_df = pd.read_csv(tableName, header=0, index_col=0, sep='\t')
	pop_df.columns.values[0] = 'population'
	return pop_df

def read_covariates(covName):
	cov_df = pd.read_csv(covName, header=0, sep='\t')
	cov_df.columns.values[0] = 'id'
	return cov_df

def gt_to_ac(gt_str):
	try:
		return int(gt_str[0])+int(gt_str[2])
	except:
		return np.nan

def get_MAF(gt_vector):
	AF = np.sum(gt_vector)/(2*len(gt_vector))
	return min(AF, 1-AF)

def openfile(filename, mode='rt'):
	if filename.endswith('.gz'):
		return gzip.open(filename, mode) 
	else:
		return open(filename, mode)


#=========================================================#
# Run eQTL interaction test (single-variant associations) #
#=========================================================#

def main():

	parser = argparse.ArgumentParser(description='Discover significant genotype-by-population interactions in an input set of QTL associations')
	parser.add_argument('associationTable', help="Four-column tab-separated file with QTL associations to test: 1) geneID, 2) variantID, 3) variantChrom, 4) variantPos")
	parser.add_argument('inVCF', help='Tabix-indexed VCF with sample genotypes')
	parser.add_argument('expressionFile', help='Comma-separated file with expression values for each sample')
	parser.add_argument('populationTable', help='Tab-separated table of sample populations/continental groups')
	parser.add_argument('outFile', help="File to write formatted eQTL data to")
	parser.add_argument('-c', '--covFile', help='Tab-separated table of covariates across samples')
	parser.add_argument('-m', '--maf_threshold', default=0.00, type=float, help='Only consider variants with MAF > threshold in at least two populations')

	args = parser.parse_args()


	#==============#
	# Parse inputs #
	#==============#

	print('\nReading inputs...\n\t1. eQTL table... ', end='', flush=True)

	eqtl_df = read_eqtls(args.associationTable)
	print('Done.\n\t2. Expression table... ', end='', flush=True)

	expr_df = read_expressions(args.expressionFile, eqtl_df)
	expr_df = expr_df.set_index('Name', drop=True).astype(float)
	print('Done.\n\t3. Genotypes from VCF... ', end='', flush=True)

	vcf_df = read_haplotypes(args.inVCF, eqtl_df)
	vcf_df = vcf_df.set_index('ID', drop=False).iloc[:, 9:].applymap(gt_to_ac)
	print('Done.\n\t4. Population table... ', end='', flush=True)

	pop_df = read_pop_table(args.populationTable)
	pops = sorted(pop_df['population'].unique().tolist())
	print('Done.')

	cov_list = []
	if args.covFile is not None:
		print('\t5. Covariates... ', end='', flush=True)
		cov_df = read_covariates(args.covFile)
		cov_df = cov_df.set_index('id', drop=True)
		cov_list = cov_df.index.tolist()
		print('Done.')

	#=========================#
	# Create design dataframe #
	#=========================#

	print('\nPreparing design matrix... ', end='', flush=True)
	full_design_df = pd.concat([expr_df.T, vcf_df.T, pop_df], axis=1)

	if args.covFile is not None:
		full_design_df = pd.concat([full_design_df, cov_df.T], axis=1)
		for cov in cov_list:
			try:
				full_design_df[cov] = full_design_df[cov].astype(float)
			except:
				pass
	print('Done.')


	#================#
	# Run regression #
	#================#

	interaction_results = pd.DataFrame({'gene_id': pd.Series(dtype=str),
										'variant_id': pd.Series(dtype=str),
										'variant_chr': pd.Series(dtype=str),
										'variant_pos': pd.Series(dtype=int),
										'interaction_pval': pd.Series(dtype=float)})

	nAssociations = eqtl_df.shape[0]

	print(f'\nRunning single-variant interaction model for {nAssociations} gene-eQTL associations', flush=True)
	print(f'\t0/{nAssociations} tests completed', end='\r', flush=True)

	for i, row in eqtl_df.iterrows():
		geneID = row['gene_id']
		variantID = row['variant_id']
		variantChr = row['variant_chr']
		variantPos = int(row['variant_pos'])
		geneID_clean = row['gene_id_clean']
		variantID_clean = row['variant_id_clean']
		
		variant_AFs = [get_MAF(full_design_df.loc[full_design_df['population']==pop, variantID_clean]) for pop in pops]
		if min(sorted(variant_AFs)[-2:]) < args.maf_threshold:
			interaction_results.loc[i, :] = [geneID, variantID, variantChr, variantPos, np.nan]
			print(f'\t{i+1}/{nAssociations} tests completed', end='\r', flush=True)
			continue

		formula_1 = f"{geneID_clean} ~ {variantID_clean} + population"
		formula_2 = f"{geneID_clean} ~ {variantID_clean}*population"
		for cov in cov_list:
			formula_1 += f' + {cov}'
			formula_2 += f' + {cov}'

		m1 = ols(formula_1, data=full_design_df).fit()
		m2 = ols(formula_2, data=full_design_df).fit()
		pval = anova_lm(m1, m2).loc[1, 'Pr(>F)']
		interaction_results.loc[i, :] = [geneID, variantID, variantChr, variantPos, pval]
		print(f'\t{i+1}/{nAssociations} tests completed', end='\r', flush=True)

	print(f'\t{i+1}/{nAssociations} tests completed\n\nAll done!\n')
	

	#==============#
	# Write output #
	#==============#

	interaction_results['variant_pos'] = interaction_results['variant_pos'].astype(int)
	interaction_results.to_csv(args.outFile, header=True, index=False, sep='\t', na_rep='NA')


if __name__ == "__main__":
	main()
