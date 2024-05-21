#!/usr/bin/env python3

from copy import deepcopy
import pandas as pd
import argparse
import os


#===========================#
# Collate sentinel variants #
#===========================#

def main():

	parser = argparse.ArgumentParser(description='Collate coloc results')
	parser.add_argument('mapFile', type=str, help='Four-column tab-separated file with trait accession, trait description, sentinel variant, coloc results prefix')
	parser.add_argument('outPrefix', type=str, help='Prefix (including path) of output files')

	args = parser.parse_args()

	#-------------------#
	# Set up output DFs #
	#-------------------#

	fullResultsDF = pd.DataFrame({'GWAS_catalog_accessionID': pd.Series(dtype='str'),
								  'traitDescription' : pd.Series(dtype='str'),
								  'GWAS_sentinelVariant' : pd.Series(dtype='str'),
								  'eQTL_gene_ensemblID': pd.Series(dtype='str'),
								  'eQTL_gene_geneSymbol': pd.Series(dtype='str'),
								  'eQTL_variantChrom': pd.Series(dtype='str'),
								  'eQTL_variantPosition': pd.Series(dtype=int),
								  'eQTL_variantRef': pd.Series(dtype='str'),
								  'eQTL_variantAlt': pd.Series(dtype='str'),
								  'eQTL_variant_kgpID': pd.Series(dtype='str'),
								  'eQTL_variant_rsID': pd.Series(dtype='str'),
								  'eQTL_geneCS': pd.Series(dtype='str'),
								  'GWAS_variantChrom': pd.Series(dtype='str'),
								  'GWAS_variantPosition': pd.Series(dtype=int),
								  'GWAS_variantRef': pd.Series(dtype='str'),
								  'GWAS_variantAlt': pd.Series(dtype='str'),
								  'GWAS_variant_kgpID': pd.Series(dtype='str'),
								  'GWAS_variant_rsID': pd.Series(dtype='str'),
								  'GWAS_CS': pd.Series(dtype='str'),
								  'PP.H0': pd.Series(dtype=float),
								  'PP.H1': pd.Series(dtype=float),
								  'PP.H2': pd.Series(dtype=float),
								  'PP.H3': pd.Series(dtype=float),
								  'PP.H4': pd.Series(dtype=float)})

	GWAS_CS_DF = pd.DataFrame({'GWAS_catalog_accessionID': pd.Series(dtype='str'),
							   'traitDescription' : pd.Series(dtype='str'),
							   'GWAS_sentinelVariant' : pd.Series(dtype='str'),
							   'variantChrom' : pd.Series(dtype='str'),
							   'variantPosition' : pd.Series(dtype=int),
							   'variantRef' : pd.Series(dtype='str'),
							   'variantAlt' : pd.Series(dtype='str'),
							   'variant_kgpID' : pd.Series(dtype='str'),
							   'variant_rsID' : pd.Series(dtype='str'),
							   'variantPIP' : pd.Series(dtype=float),
							   'variantCredibleSet' : pd.Series(dtype='str'),
							   })

	eQTL_CS_DF = pd.DataFrame({'GWAS_catalog_accessionID': pd.Series(dtype='str'),
							   'traitDescription' : pd.Series(dtype='str'),
							   'GWAS_sentinelVariant' : pd.Series(dtype='str'),
							   'variantChrom' : pd.Series(dtype='str'),
							   'variantPosition' : pd.Series(dtype=int),
							   'variantRef' : pd.Series(dtype='str'),
							   'variantAlt' : pd.Series(dtype='str'),
							   'variant_kgpID' : pd.Series(dtype='str'),
							   'variant_rsID' : pd.Series(dtype='str'),
							   'ensemblID': pd.Series(dtype='str'),
							   'geneSymbol': pd.Series(dtype='str'),
							   'variantPIP' : pd.Series(dtype=float),
							   'variantCredibleSet' : pd.Series(dtype='str'),
							   })


	#------------------#
	# Get full results #
	#------------------#

	with open(args.mapFile) as in_fs:
		for line in in_fs:
			fields = line.rstrip('\n').split('\t')
			traitID = fields[0]
			traitDescription = fields[1]
			sentinelVariant = fields[2]
			fpath = fields[3] + ".coloc.fullResults.txt"
			
			if os.path.exists(fpath):
				signalFullResultsDF = pd.read_csv(fpath, header=0, sep='\t')
				signalFullResultsDF['GWAS_catalog_accessionID'] = traitID
				signalFullResultsDF['traitDescription'] = traitDescription
				signalFullResultsDF['GWAS_sentinelVariant'] = sentinelVariant
				signalFullResultsDF = signalFullResultsDF.loc[:, fullResultsDF.columns]

				fullResultsDF = pd.concat([fullResultsDF, signalFullResultsDF], axis=0)

	# Output full results
	fullResultsDF = fullResultsDF.sort_values(by=['GWAS_catalog_accessionID'])
	fullResultsDF.to_csv(f'{args.outPrefix}.full_results.txt', header=True, index=False, sep='\t', na_rep='NA')

	# Output moderate results
	moderateResultsDF = fullResultsDF.loc[fullResultsDF['PP.H4'] >= 0.5, :]
	moderateResultsDF.to_csv(f'{args.outPrefix}.moderate_results.txt', header=True, index=False, sep='\t', na_rep='NA')


	#------------------------#
	# Get GWAS credible sets #
	#------------------------#

	with open(args.mapFile) as in_fs:
		for line in in_fs:
			fields = line.rstrip('\n').split('\t')
			traitID = fields[0]
			traitDescription = fields[1]
			sentinelVariant = fields[2]
			fpath = fields[3] + ".coloc.GWAS_CS.txt"

			if os.path.exists(fpath):
				signalGWAS_CS_DF = pd.read_csv(fpath, header=0, sep='\t')
				signalGWAS_CS_DF['GWAS_catalog_accessionID'] = traitID
				signalGWAS_CS_DF['traitDescription'] = traitDescription
				signalGWAS_CS_DF['GWAS_sentinelVariant'] = sentinelVariant
				signalGWAS_CS_DF = signalGWAS_CS_DF.loc[:, GWAS_CS_DF.columns]

				GWAS_CS_DF = pd.concat([GWAS_CS_DF, signalGWAS_CS_DF], axis=0)

	GWAS_CS_DF = GWAS_CS_DF.sort_values(by=['GWAS_catalog_accessionID'])
	GWAS_CS_DF.to_csv(f'{args.outPrefix}.moderate_GWAS_CS.txt', header=True, index=False, sep='\t', na_rep='NA')


	#------------------------#
	# Get eQTL credible sets #
	#------------------------#

	with open(args.mapFile) as in_fs:
		for line in in_fs:
			fields = line.rstrip('\n').split('\t')
			traitID = fields[0]
			traitDescription = fields[1]
			sentinelVariant = fields[2]
			fpath = fields[3] + ".coloc.eQTL_CS.txt"

			if os.path.exists(fpath):
				signaleQTL_CS_DF = pd.read_csv(fpath, header=0, sep='\t')
				signaleQTL_CS_DF['GWAS_catalog_accessionID'] = traitID
				signaleQTL_CS_DF['traitDescription'] = traitDescription
				signaleQTL_CS_DF['GWAS_sentinelVariant'] = sentinelVariant
				signaleQTL_CS_DF = signaleQTL_CS_DF.loc[:, eQTL_CS_DF.columns]

				eQTL_CS_DF = pd.concat([eQTL_CS_DF, signaleQTL_CS_DF], axis=0)

	eQTL_CS_DF = eQTL_CS_DF.sort_values(by=['GWAS_catalog_accessionID'])
	eQTL_CS_DF.to_csv(f'{args.outPrefix}.moderate_eQTL_CS.txt', header=True, index=False, sep='\t', na_rep='NA')


if __name__ == "__main__":
	main()
