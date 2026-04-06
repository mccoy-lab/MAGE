#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import subprocess
import os
import gzip
import contextlib
from datetime import datetime
import tempfile
import shutil
import glob
from sklearn.decomposition import PCA
import qtl.io
import sys

def is_gzipped(path):
	with open(path, "rb") as f:
		return f.read(2) == b'\x1f\x8b'

class Unbuffered(object):
	def __init__(self, stream):
		self.stream = stream

	def write(self, data):
		self.stream.write(data)
		self.stream.flush()

	def writelines(self, datas):
		self.stream.writelines(datas)
		self.stream.flush()

	def __getattr__(self, attr):
		return getattr(self.stream, attr)

sys.stdout = Unbuffered(sys.stdout)

if __name__=='__main__':

	parser = argparse.ArgumentParser(description='Filter Leafcutter output, prepare for FastQTL')
	parser.add_argument('leafcutter_output_dir', help='Path to Leafcutter output directory')
	parser.add_argument('leafcutter_output_prefix', help='Prefix of Leafcutter output files (i.e. everything before `_perind_numers.counts.gz` and `_perind.counts.gz`)')
	parser.add_argument('sample_participant_lookup', help='Lookup table linking samples to participants')
	parser.add_argument('genes_gtf', help='GTF file with gene annotations')
	parser.add_argument('cluster_to_gene', help='File containing map from clusterID to geneID')
	parser.add_argument('--keep_samples', default=None, type=str, help='File containing the samples to keep in the output dataframe. If this argument is not used, all samples will be included')
	parser.add_argument('--pct_cutoff', default=50, type=int, help='The minimum percent of (kept) samples that introns must have counts > 0 in. E.g. "50" for a 50% threshold')
	parser.add_argument('-l', '--leafcutter_prep_script', default='./prepare_phenotype_table.py', help='Path to leafcutter `prepare_phenotype_table.py` helper script')
	parser.add_argument('-o', '--output_dir', default='.', help='Output directory for filtered counts file.')
	parser.add_argument('-a', '--autosomes_chrX', action='store_true', help='If included, will limit output to introns on autosomes + chrX')
	args = parser.parse_args()

	print(f'[{datetime.now().strftime("%b %d %H:%M:%S")}] leafcutter filtering')

	print('  * preparing data')
	counts_df = pd.read_csv(args.leafcutter_output_dir+'/'+args.leafcutter_output_prefix+'_perind.counts.gz', sep='\s+').set_index('chrom')
	col_dict = {i:i.split('.')[0] for i in counts_df.columns}
	counts_df.rename(columns=col_dict, inplace=True)
	sample_participant_lookup_s = pd.read_csv(args.sample_participant_lookup,
										  sep='\t', index_col=0, dtype=str, header=None).squeeze("columns")
	assert counts_df.columns.isin(sample_participant_lookup_s.index).all()
	counts_df.rename(columns=sample_participant_lookup_s, inplace=True)

	if args.keep_samples is not None:
		print('    ** Subsetting to defined samples')
		with open(args.keep_samples) as fs:
			ncols = len(fs.readline().strip().split('\t'))
		if ncols == 1:
			samples = [line.strip() for line in open(args.keep_samples)]
			counts_df = counts_df.loc[:,samples]
		elif ncols == 2:
			samples = {line.strip().split('\t')[0]:line.strip().split('\t')[1] for line in open(args.keep_samples)}
			counts_df = counts_df.loc[:, list(samples.keys())]
			counts_df.rename(columns=samples, inplace=True)
		else:
			raise Exception('"--keep_samples": File must have 1 or 2 columns')

	# Removing introns on chrY
	if args.autosomes_chrX:
		prefilter_clusters = len(set(np.unique(counts_df.index.map(lambda x: x.split(':')[-1]))))
		prefilter_introns = counts_df.shape[0]
		chroms = tuple(f'chr{i}' for i in range(1, 23))+('chrX',)
		autosomes_chrX_mask = counts_df.index.str.startswith(chroms)
		removed_introns = np.size(autosomes_chrX_mask) - np.count_nonzero(autosomes_chrX_mask)
		counts_df = counts_df.loc[autosomes_chrX_mask]
		introns = counts_df.shape[0]
		clusters = len(set(np.unique(counts_df.index.map(lambda x: x.split(':')[-1]))))
		print('    ** dropping {} introns not belonging to autosomes or chrX\n'
			  '       {}/{} introns remain ({}/{} clusters)'.format(
				  removed_introns, introns,
				   prefilter_introns, clusters, prefilter_clusters)
			)



	print('  * filtering counts')
	calculate_frac = lambda x: float(x[0])/float(x[1]) if x[1] > 0 else 0
	get_count = lambda x: float(x[0]) if x[1] > 0 else 0
	frac_df = counts_df.applymap(lambda x: calculate_frac([int(i) for i in x.split('/')]))
	pct_zero = (frac_df == 0).sum(1) / frac_df.shape[1]  # for zero counts, frac is zero
	n_unique = frac_df.apply(lambda x: len(x.unique()), axis=1)
	zscore_df = ((frac_df.T-frac_df.mean(1)) / frac_df.std(1)).T

	# filter out introns with low counts or low complexity
	n = np.floor(frac_df.shape[1]*0.1)
	if n < 10:
		n = 10
	mask = (pct_zero <= (100-args.pct_cutoff)/100) & (n_unique >= n) # <- change this to change percent sample threshold

	# additional filter for low complexity
	ns = zscore_df.shape[1]
	mask2 = ((zscore_df.abs()<0.25).sum(1) >= ns-3) & ((zscore_df.abs() > 6).sum(1) <= 3)
	if np.any(mask & mask2):
		print(f'    ** dropping {np.sum(mask & mask2)} introns with low variation')
	mask = mask & ~mask2

	filtered_counts_df = counts_df.loc[mask].copy()

	cluster_ids = np.unique(counts_df.index.map(lambda x: x.split(':')[-1]))
	filtered_cluster_ids = np.unique(filtered_counts_df.index.map(lambda x: x.split(':')[-1]))
	print('    ** dropping {} introns with counts in fewer than {}% of samples\n'
		  '       {}/{} introns remain ({}/{} clusters)'.format(
			   counts_df.shape[0]-filtered_counts_df.shape[0], args.pct_cutoff, filtered_counts_df.shape[0],
			   counts_df.shape[0], len(filtered_cluster_ids), len(cluster_ids))
		)
	
	temp_counts_file = os.path.join(args.output_dir, f'tmp_{args.leafcutter_output_prefix}_perind.counts.filtered.gz')
	filtered_counts_df.to_csv(temp_counts_file, sep=' ', header=True, index=True, compression='gzip')

	# preparing phenotype table
	subprocess.check_call(
		'python3 '+os.path.join(args.leafcutter_prep_script) \
		+f' {temp_counts_file}' \
		+f' -p 0', shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

	# concatenating chromosome-level BED files
	bed_files = sorted(glob.glob(os.path.join(args.output_dir, '*_perind.counts.filtered.gz.qqnorm_*')))
	bed_df = []
	for f in bed_files:
		bed_df.append(pd.read_csv(f, sep='\t', dtype=str))
	bed_df = pd.concat(bed_df, axis=0)
	bed_df['chr_ix'] = bed_df['#Chr'].str.replace('chr','').str.replace('X','23').str.replace('Y', '24').astype(np.int32)
	for c in ['start', 'end']:
		bed_df[c] = bed_df[c].astype(np.int32)
	bed_df.sort_values(['chr_ix', 'start', 'end'], inplace=True)
	bed_df.drop('chr_ix', axis=1, inplace=True)
	bed_df.rename(columns={'#Chr':'#chr'}, inplace=True)
	
	# writing merged BED
	bed_file = os.path.join(args.output_dir, f'{args.leafcutter_output_prefix}.perind.counts.filtered.qqnorm.bed.gz')
	qtl.io.write_bed(bed_df, bed_file)

	# Removing introns filtered during "prepare_phenotype_table" script
	filtered_introns = bed_df['ID'].to_list()
	double_filtered_counts_df = filtered_counts_df.loc[filtered_introns, :]
	double_filtered_cluster_ids = np.unique(double_filtered_counts_df.index.map(lambda x: x.split(':')[-1]))
	print('    ** dropping {} introns whose cluster is missing in 40% of samples OR with low variation across samples\n'
		  '       {}/{} introns remain ({}/{} clusters)'.format(
			   filtered_counts_df.shape[0]-double_filtered_counts_df.shape[0], double_filtered_counts_df.shape[0],
			   counts_df.shape[0], len(double_filtered_cluster_ids), len(cluster_ids))
		)

	# Removing clusters with only one intron
	one_intron_mask = double_filtered_counts_df.index.map(lambda x: x.split(':')[-1]).duplicated(keep=False)
	final_filtered_counts_df = double_filtered_counts_df.loc[one_intron_mask]
	final_filtered_cluster_ids = set(np.unique(final_filtered_counts_df.index.map(lambda x: x.split(':')[-1])))
	print('    ** dropping {} clusters with only one intron\n'
		  '       {}/{} introns remain ({}/{} clusters)'.format(
			  double_filtered_counts_df.shape[0]-final_filtered_counts_df.shape[0], final_filtered_counts_df.shape[0],
			   counts_df.shape[0], len(final_filtered_cluster_ids), len(cluster_ids))
		)

	print('  * converting cluster coordinates to gene coordinates')
	tss_df = qtl.io.gtf_to_tss_bed(args.genes_gtf, )
	cluster2gene_dict = pd.read_csv(args.cluster_to_gene,
		sep='\t', index_col=0).squeeze("columns").to_dict()

	print('    ** assigning introns to gene mapping(s)')
	n = 0
	gene_bed_df = []
	group_s = {}
	for _,r in bed_df.iterrows():
		s = r['ID'].split(':')
		if s[-1] not in final_filtered_cluster_ids: # filter out clusters with one intron
			continue
		cluster_id = s[0]+':'+s[-1]
		if cluster_id in cluster2gene_dict:
			gene_ids = cluster2gene_dict[cluster_id].split(',')
			for g in gene_ids:
				gi = r['ID']+':'+g
				gene_bed_df.append(tss_df.loc[g, ['chr', 'start', 'end']].tolist() + [gi] + r.iloc[4:].tolist())
				group_s[gi] = g
		else:
			n += 1
	if n > 0:
		print(f'    ** discarded {n} introns without a gene mapping')

	print('  * writing BED files for QTL mapping')
	gene_bed_df = pd.DataFrame(gene_bed_df, columns=bed_df.columns)
	# sort by TSS and gene
	gene_bed_df['gene'] = gene_bed_df['ID'].apply(lambda x: x.split(':')[-1])
	gene_bed_df = gene_bed_df.groupby('#chr', sort=False, group_keys=False).apply(lambda x: x.sort_values(['start', 'gene', 'ID']))
	gene_bed_df.drop('gene', axis=1, inplace=True)
	qtl.io.write_bed(gene_bed_df, os.path.join(args.output_dir, f'{args.leafcutter_output_prefix}_perind.counts.normed.genes.filtered.bed.gz'))
	gene_bed_df[['start', 'end']] = gene_bed_df[['start', 'end']].astype(np.int32)
	gene_bed_df[gene_bed_df.columns[4:]] = gene_bed_df[gene_bed_df.columns[4:]].astype(np.float32)
	pd.Series(group_s).sort_values().to_csv(
		os.path.join(args.output_dir, f'{args.leafcutter_output_prefix}_perind.counts.normed.genes.filtered.phenotype_groups.txt'),
		sep='\t', header=False)

	print('  * writing filtered leafcutter files for Manta')
	final_filtered_counts_df.insert(loc=0, column='cluster', value=final_filtered_counts_df.index.map(lambda x: x.split(':')[0] + ':' + x.split(':')[-1]), allow_duplicates=True)
	final_filtered_counts_df.index.name = 'intron:cluster'
	final_filtered_counts_df = final_filtered_counts_df.reset_index()
	final_filtered_counts_df.sort_values(by=['cluster', 'intron:cluster'], inplace=True)
	
	filtered_counts_file = os.path.join(args.output_dir, f'{args.leafcutter_output_prefix}_perind.counts.filtered.gz')
	final_filtered_counts_df.to_csv(filtered_counts_file, sep='\t', header=True, index=False)

	# delete intermediary files
	files = glob.glob(os.path.join(args.output_dir, f'tmp_{args.leafcutter_output_prefix}_perind.counts.filtered.gz.qqnorm_chr*')) \
		  + glob.glob(os.path.join(args.output_dir, f'tmp_{args.leafcutter_output_prefix}_perind.counts.filtered.gz.phen_chr*'))
	for f in files:
		os.remove(f)
	os.remove(os.path.join(args.output_dir, f'tmp_{args.leafcutter_output_prefix}_perind.counts.filtered.gz_prepare.sh'))
	os.remove(os.path.join(args.output_dir, f'{args.leafcutter_output_prefix}.perind.counts.filtered.qqnorm.bed.gz'))
	os.remove(os.path.join(args.output_dir, f'{args.leafcutter_output_prefix}.perind.counts.filtered.qqnorm.bed.gz.tbi'))
	os.remove(os.path.join(args.output_dir, f'tmp_{args.leafcutter_output_prefix}_perind.counts.filtered.gz.ave'))
	os.remove(temp_counts_file)

	print(f'[{datetime.now().strftime("%b %d %H:%M:%S")}] done')
