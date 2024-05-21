#!/usr/bin/env python3

import argparse
import warnings
import gzip
import pysam
import pandas as pd
import time
import statsmodels.formula.api as smf
import copy

#==================#
# Define functions #
#==================#

def sample_column_map(path, start_col=9, line_key="#CHR"):
	stream_in = gzip.open(path, "r")

	out_map = {}
	for line in stream_in:
		if isinstance(line, bytes) and not isinstance(line, str):
			line = line.decode()
		if line_key in line:
			line = line.rstrip().split("\t")
			for i in range(start_col,len(line)):
				out_map[line[i]] = i

			break

	stream_in.close()

	return(out_map)


def correct_covariates(df_test):

	# correct for covariates
	# add genotype categorical covariates
	genoIDs = [x for x in df_test.columns if x.startswith('geno')]
	n_genos = len(genoIDs)
	for genoID in genoIDs:
		i = int(genoID.split('_')[1])
		cov_homo_ref = [int(x == 0) for x in df_test[genoID]]
		if sum(cov_homo_ref) > 0:
			df_test[f'cov_homo_ref_{i}'] = cov_homo_ref

		cov_homo_alt = [int(x == 2) for x in df_test[genoID]]
		if sum(cov_homo_alt) > 0:
			df_test[f'cov_homo_alt_{i}'] = cov_homo_alt

	cov_ids = [x for x in df_test.columns if "cov_" in x]

	# convert categorical covariates to n-1 binary covariates
	new_cols = {}
	drop_cols = []

	for xcov in cov_ids:
		if df_test.dtypes[xcov] == object:
			values = list(set(df_test[xcov]))[1:]
			for xval in values:
				xname = xcov+"_"+xval
				new_cols[xname] = [int(x == xval) for x in df_test[xcov]]

			drop_cols.append(xcov)

	df_test.drop(drop_cols,axis=1,inplace=True)
	for xcov in new_cols.keys():
		df_test[xcov] = new_cols[xcov]
	cov_ids = [x for x in df_test.columns if "cov_" in x]

	# NOTE any variable that is a string will be treated as categorical - this is the same functionality as FASTQTL, so good
	# see: http://statsmodels.sourceforge.net/devel/example_formulas.html

	xformula = "pheno ~ "+"+".join(cov_ids)
	result = smf.ols(formula=xformula, data=df_test).fit()

	# use only significant (95% CI doesn't overlap 0) covariates to correct expression values
	# do not include intercept or genotypes in correction

	drop_covs = []
	for xcov in list(result.params.index):
		if xcov.startswith("cov_homo_ref") or xcov.startswith("cov_homo_alt"):
			continue
		if xcov in df_test.columns:
			coefficient = result.params.loc[xcov]
			# pval = result.params.loc[xcov]
			# if pval >= 0.01:
			# 	drop_covs.append(xcov) 
			upper_ci = result.conf_int(0.05).loc[xcov][1]
			lower_ci = result.conf_int(0.05).loc[xcov][0]
			if (lower_ci <= 0 and upper_ci >= 0):
				drop_covs.append(xcov)

	# drop insignificant covariates
	df_test.drop(drop_covs, axis=1, inplace=True)
	cov_ids = [x for x in df_test.columns if "cov_" in x]

	# redo regression without insignificant covs
	if len(cov_ids) > 0:
		xformula = "pheno ~ "+"+".join(cov_ids)
		result = smf.ols(formula=xformula, data=df_test).fit()

		df_test_corrected = copy.deepcopy(df_test)
		for xcov in list(result.params.index):
			coefficient = result.params.loc[xcov]
			if xcov == "Intercept" or xcov.startswith("cov_homo_ref") or xcov.startswith("cov_homo_alt"):
				df_test_corrected[xcov] = [0] * len(df_test_corrected.index)
			else:
				df_test_corrected[xcov] = [x * coefficient for x in df_test_corrected[xcov]]

		# add residual to dataframe
		df_test_corrected['pheno_cor'] = [row['pheno'] - sum(row[1+n_genos:len(row)]) for index, row in df_test_corrected.iterrows()]

	else:
		# if none of the covariates are significant then just leave the values as is
		df_test_corrected = copy.deepcopy(df_test)
		df_test_corrected['pheno_cor'] = df_test_corrected['pheno']


	return(df_test_corrected)


#========================#
# Regress out covariates #
#========================#

def main():

	parser = argparse.ArgumentParser(description="From input gene expression quantifications, regress out all covariates whose effect size 95% CI exludes 0")
	parser.add_argument("--vcf", required=True, help="Genotype VCF")
	parser.add_argument("--pheno", required=True, help="Phenotype file (bed format)")
	parser.add_argument("--qtl", required=True, help="File containing QTL to calculate allelic fold change for. Should contain tab separated columns 'gene_id' with phenotype (gene) IDs and 'variant_id' with SNP IDs. Optionally can include the columns 'variant_chr' and 'variant_pos', which will facilitate tabix retrieval of genotypes, greatly reducing runtime.")
	parser.add_argument("--cov", required=True, help="Covariates file (in the same format used by FastQTL)")
	parser.add_argument("--output", required=True, help="Output file")
	args = parser.parse_args()

	# # disable warnings
	# warnings.filterwarnings("ignore")


	#=========================#
	# Get sample map from VCF #
	#=========================#

	vcf_map = sample_column_map(args.vcf)
	tabix_vcf = pysam.Tabixfile(args.vcf,"r")


	#=================#
	# Load Covariates #
	#=================#

	df_cov = pd.read_csv(args.cov, sep="\t", index_col=False)
	if "ID" in df_cov.columns:
		cov_id_col = "ID"
	elif "id" in df_cov.columns:
		cov_id_col = "id"
	else:
		print("Could not find covariate ID column in covariates column. Please ensure that it is either labeled 'ID' or 'id'")
		quit()

	trans_df_cov = df_cov.T
	trans_df_cov.columns = trans_df_cov.iloc[0]
	trans_df_cov = trans_df_cov.iloc[1:]
	trans_df_cov.columns.name = None
	trans_df_cov.index.rename('id', inplace=True)
	trans_df_cov.rename(columns={x: f'cov_{x}' for x in trans_df_cov.columns}, inplace=True)


	#=====================#
	# Load phenotype data #
	#=====================#

	pheno_map = sample_column_map(args.pheno, line_key="#", start_col=4)
	tabix_pheno = pysam.Tabixfile(args.pheno, "r")


	#======================#
	# Load in QTLs to test #
	#======================#

	df_qtl = pd.read_csv(args.qtl, sep="\t", index_col=False)


	#========================#
	# Load SNP position info #
	#========================#

	set_esnp = set(df_qtl['variant_id'].tolist())
	dict_esnp = {}

	for index, row in df_qtl.iterrows():
		dict_esnp[row['variant_id']] = [row['variant_chr'],int(row['variant_pos'])]


	#====================#
	# Get gene positions #
	#====================#

	set_epheno = set(df_qtl['gene_id'].tolist())
	dict_ephenotype = {}
	orderedGenes = []

	with gzip.open(args.pheno, "rt") as stream_in:
		for line in stream_in:
			if line[0:1] != "#":
				columns = line.rstrip().split("\t")
				#Chr    start   end     ID
				if columns[3] in set_epheno:
					orderedGenes.append(columns[3])
					dict_ephenotype[columns[3]] = [columns[0],int(columns[1]), int(columns[2])]


	#====================#
	# Correct covariates #
	#====================#

	stream_out = open(args.output, "w")
	
	with gzip.open(args.pheno, 'rt') as stream_in:
		for line in stream_in:
			orderedSamples = line.strip('\n').split('\t')[4:]
			stream_out.write(line)
			break

	nGenes = len(df_qtl['gene_id'].unique())
	completed = 0
	start_time = time.time()

	df_qtl['geneOrder'] = df_qtl['gene_id'].apply(lambda x: orderedGenes.index(x))
	df_qtl = df_qtl.sort_values(by='geneOrder')
	df_qtl.drop(columns=['geneOrder'], inplace=True)

	for geneID, gene_df_qtl in df_qtl.groupby('gene_id', sort=False):

		ephenotype = dict_ephenotype[geneID]
		records = tabix_pheno.fetch(ephenotype[0], ephenotype[1]-1, ephenotype[1]+1)
		dict_pheno = {}
		
		for record in records:
			cols = record.rstrip().split("\t")
			if cols[3] == geneID:
				for sample in pheno_map.keys():
					dict_pheno[sample] = float(cols[pheno_map[sample]])
		
		df_test = pd.Series(dict_pheno).rename(index='pheno').to_frame()
		df_test.index.rename('sampleID', inplace=True)

		dict_snpInfo = {}

		for index, row in gene_df_qtl.iterrows():
			
			dict_snpInfo[f'geno_{index}'] = row.tolist()[1:]
			esnp = dict_esnp[row['variant_id']]
			records = tabix_vcf.fetch(esnp[0], esnp[1]-1, esnp[1])
			dict_geno = {}

			for record in records:
				cols = record.rstrip().split("\t")
				if cols[2] == row['variant_id']:
					gt_index = cols[8].split(":").index("GT")
					snp_found = 1
					for sample in pheno_map.keys():
						sample_col = cols[vcf_map[sample]]
						dict_geno[sample] = sample_col.split(":")[gt_index]

			out_dict_geno = {}

			for sample in dict_geno.keys():
				if "." not in dict_geno[sample]:	# only include samples w/ complete genotype data (no '.')
					out_dict_geno[sample] = dict_geno[sample].count("1")
				else:
					out_dict_geno[sample] = np.nan

			df_test[f'geno_{index}'] = pd.Series(out_dict_geno)

		df_test = pd.concat([df_test, trans_df_cov], join='inner', axis=1)
		df_test = df_test.apply(pd.to_numeric, errors='ignore')

		# Do the actual correction
		df_test_corrected = correct_covariates(df_test)
		orderedCorrectExpression = [str(x) for x in df_test_corrected.loc[orderedSamples, 'pheno_cor'].tolist()]
		stream_out.write(f'{ephenotype[0]}\t{ephenotype[1]}\t{ephenotype[2]}\t{geneID}\t' + '\t'.join(orderedCorrectExpression) + '\n')

		completed += 1
		if completed % 10 == 0:
			stream_out.flush()
			print("     COMPLETED %d of %d = %f in %d seconds"%(
			    completed, nGenes,
			    float(completed)/float(nGenes),
			    time.time()-start_time))


	stream_out.close()


if __name__ == "__main__":
	main()
