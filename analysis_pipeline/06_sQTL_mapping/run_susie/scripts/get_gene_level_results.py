#!/usr/bin/env python3

import argparse
import pandas as pd
from copy import deepcopy

#==================#
# Define functions #
#==================#

def merge_sets(list_of_sets):
	pooled = deepcopy(list_of_sets)
	merging = True
	while merging:
		merging=False
		for i,group in enumerate(pooled):
			matching_group = next((g for g in pooled[i+1:] if g.intersection(group)),None)
			if not matching_group: continue
			group.update(matching_group)
			pooled.remove(matching_group)
			merging = True
	return pooled


#========================#
# Get gene-level results #
#========================#

def main():

	parser = argparse.ArgumentParser(description='Merge overlapping credible sets into a gene-wise summary from intron-wise sQTL fine-mapping results')
	parser.add_argument('intronResults', help="Tab-separated file with significant sQTL fine-mapping associations, one row for each intron-variant association")
	parser.add_argument('credibleSetResults', help="Tab-separated file with intron-level credible set summary, one row for each credible set")
	parser.add_argument('outGeneResults', help='File name to write output gene-wise merged credible set results to (all included variants)')
	parser.add_argument('outLeadResults', help='File name to write out gene-wise lead sVariant results to (one variant per merged credible set)')
	args = parser.parse_args()

	#==============#
	# Read in data #
	#==============#

	# All intron-level credible set results
	intronResults = pd.read_csv(args.intronResults, sep='\t', header=0)
	intronResults['intronCredibleSet'] = intronResults['intronID'] + ',' + intronResults['variantCredibleSet']

	# Intron-level credible set summary
	credibleSetResults = pd.read_csv(args.credibleSetResults, sep='\t', header=0)
	credibleSetResults['intronCredibleSet'] = credibleSetResults['intronID'] + ',' + credibleSetResults['credibleSet']
	credibleSetResults.set_index('intronCredibleSet', inplace=True)


	#===================#
	# Create output DFs #
	#===================#

	independentSetsResults = pd.DataFrame({'geneID': pd.Series(dtype='str'),
										   'variantID': pd.Series(dtype='str'),
										   'mergedCredibleSet': pd.Series(dtype='str'),
										   'intronIDs': pd.Series(dtype='str')})

	independentSetLeadResults = pd.DataFrame({'geneID': pd.Series(dtype='str'),
											  'variantID': pd.Series(dtype='str'),
											  'mergedCredibleSet': pd.Series(dtype='str'),
											  'intronID': pd.Series(dtype='str'),
											  'intronCredibleSet': pd.Series(dtype='str'),
											  'intronCS_coverage': pd.Series(dtype=float),
											  'intronCS_variantPIP': pd.Series(dtype=float),
											  'mergedIntronSets': pd.Series(dtype=str)})

	#===========#
	# Grab data #
	#===========#

	numGenes = len(intronResults['geneID'].unique())
	for gi, geneID in enumerate(intronResults['geneID'].unique()):
		geneResults = intronResults.loc[intronResults['geneID']==geneID]
		geneCredibleSets = []
		for intronID in geneResults['intronID'].unique():
			geneIntronResults = geneResults.loc[geneResults['intronID']==intronID]
			for credibleSet in geneIntronResults['variantCredibleSet'].unique():
				geneCredibleSets.append(set(geneIntronResults.loc[geneIntronResults['variantCredibleSet']==credibleSet, 'variantID']))
		mergedGeneCredibleSets = merge_sets(geneCredibleSets)

		for i, mergedSet in enumerate(mergedGeneCredibleSets):
			mergedSetID = f'M{i+1}'
			mergedSetIntrons = set()
			for variantID in sorted(mergedSet, key=lambda x: int(x.split(':')[1])):
				intronIDs = set(geneResults.loc[geneResults['variantID']==variantID, 'intronCredibleSet'].tolist())
				mergedSetIntrons |= intronIDs
				intronIDs_string = ";".join(sorted(list(intronIDs)))
				independentSetsResults.loc[len(independentSetsResults)] = [geneID, variantID, mergedSetID, intronIDs_string]
			mergedSetIntrons = sorted(list(mergedSetIntrons))
			mergedSetIntrons_string = ";".join(mergedSetIntrons)
			sortedIntronCredibleSets = credibleSetResults.loc[mergedSetIntrons, :].sort_values(by='coverage')
			leadIntronID = sortedIntronCredibleSets['intronID'][-1]
			leadIntronCredibleSet = sortedIntronCredibleSets['credibleSet'][-1]
			leadIntronCredibleSetCoverage = sortedIntronCredibleSets['coverage'][-1]
			sortedCredibleSetVariants = geneResults.loc[geneResults['intronCredibleSet']==leadIntronID+','+leadIntronCredibleSet, :].sort_values(by='variantPIP')
			leadVariantID = sortedCredibleSetVariants['variantID'].iloc[-1]
			leadVariantPIP = sortedCredibleSetVariants['variantPIP'].iloc[-1]
			independentSetLeadResults.loc[len(independentSetLeadResults)] = [geneID, leadVariantID, mergedSetID, leadIntronID, leadIntronCredibleSet, leadIntronCredibleSetCoverage, leadVariantPIP, mergedSetIntrons_string]
		
		if gi % 10 == 0:
			print(f'{gi}/{numGenes}')


	#=================#
	# Write to output #
	#=================#

	independentSetsResults.to_csv(args.outGeneResults, sep='\t', header=True, index=False)
	independentSetLeadResults.to_csv(args.outLeadResults, sep='\t', header=True, index=False)


if __name__ == '__main__':
	main()
