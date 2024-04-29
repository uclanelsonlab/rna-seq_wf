#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from os.path import basename
import gc
import resource

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
parser.add_argument('--proband_bam2sj', required=True, help='input file')
parser.add_argument('--chromosome', required=True, help='limit to given chromosome')
parser.add_argument('--reference_bam2sj_hdf5', required=True, help='per-chromosome reference bam2sj hdf5')
parser.add_argument('--sjdb', required=True, help='sjdbList.fromGTF.out.tab')
parser.add_argument('--minimum_cutoff', default=4, type=int, help='minimum read cutoff')
parser.add_argument('--udn_id_key', help='path to udn id key')
parser.add_argument('--output_dir', required=True, help='output dir')
args = parser.parse_args()

chrom = args.chromosome
sjdb_file = args.sjdb
minimum_cutoff = args.minimum_cutoff
sample_name = basename(args.proband_bam2sj)
if sample_name.endswith('.bam2SJ.out.tab.gz'):
  sample_name = sample_name[:-len('.bam2SJ.out.tab.gz')]

# load udn id key
udn_dict = dict()
with open(args.udn_id_key) as f:
  for l in f:
    line_list = l.strip().split('\t')
    udn_dict[line_list[0]] = line_list[1]

# function to get proband id from sample name
tissue_types = ['muscle', 'skin', 'fibroblast', 'blood', 'transdiff-neurons', 'TdF', '2wkDiffNeu', 'NSC', 'fat', 'liver', 'bonemarrow', 'heart']
def get_proband_id(column_name):
  sample_name = column_name.split('.')[0]
  for tissue_type in tissue_types:
    sample_name = sample_name.split('-' + tissue_type)[0]
  if sample_name.startswith('IPH'):
    return sample_name[:8] + 'P'
  elif sample_name.startswith('UDN'):
    return udn_dict.get(sample_name[:9], sample_name)
  else:
    return sample_name

# function to add canonical sjdb annotation to data frame index
def add_canonical_annotations_to_index(input_df):
  global sjdb_file
  column_names = ['chr', 'start', 'end', 'strand']
  column_dtypes = ['str', 'int', 'int', 'str']
  column_dtype_dict = dict(zip(column_names, column_dtypes))
  sjdb_list = pd.read_table(sjdb_file, names=column_names, index_col=[0, 1, 2], dtype=column_dtype_dict)
  sjdb_list['is_canonical'] = True
  sjdb_list = sjdb_list.drop('strand', axis=1) # skip for now
  sjdb_list = sjdb_list[~sjdb_list.index.duplicated(keep='first')]
  canonical_starts = sjdb_list.groupby(['chr', 'start']).any()
  canonical_ends = sjdb_list.groupby(['chr', 'end']).any()
  
  input_df_index = pd.DataFrame(index=input_df.index)
  input_df_with_annotations = input_df_index.join(
    sjdb_list, on=['chr', 'start', 'end']).join(
    canonical_starts, on=['chr', 'start'], rsuffix='_start').join(
    canonical_ends, on=['chr', 'end'], rsuffix='_end').fillna(False)
  input_df_with_annotations = input_df_with_annotations.set_index(['is_canonical', 'is_canonical_start', 'is_canonical_end'], append=True)
  input_df_with_annotations = input_df.reindex(index=input_df_with_annotations.index)
  return input_df_with_annotations

print('Resource checkpoint 1:', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)

# load the input bam2sj file:
column_names = ['chr', 'start', 'end', 'n_unique']
column_dtypes = ['str'] + ['int']*3
column_dtype_dict = dict(zip(column_names, column_dtypes))
df_s = pd.read_table(args.proband_bam2sj, names=column_names, index_col=[0, 1, 2], dtype=column_dtype_dict)
df_s = df_s.query('chr == "' + chrom + '"')
df_s_with_minimum_cutoff = add_canonical_annotations_to_index(df_s).fillna(0).replace(range(1, minimum_cutoff), 0)
df_s_with_minimum_cutoff = df_s_with_minimum_cutoff.loc[(df_s_with_minimum_cutoff != 0).any(axis=1)]

# load the reference file
df = pd.read_hdf(args.reference_bam2sj_hdf5, 'splicing_master_table')
df = df.loc[df.index.intersection(df_s.index)] # https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#indexing-deprecate-loc-reindex-listlike
df_with_minimum_cutoff = add_canonical_annotations_to_index(df).fillna(0).replace(range(1, minimum_cutoff), 0)
df_with_minimum_cutoff.loc[(df_with_minimum_cutoff != 0).any(axis=1)]

print('Resource checkpoint 2:', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
del df
del df_s
gc.collect()

# make ratios
df_start_counts = df_with_minimum_cutoff.query('is_canonical == True').groupby(['chr', 'start']).sum()
df_end_counts = df_with_minimum_cutoff.query('is_canonical == True').groupby(['chr', 'end']).sum()
df_start_ratios = df_with_minimum_cutoff/(df_start_counts.reindex(df_with_minimum_cutoff.index).fillna(0) + df_with_minimum_cutoff)
df_end_ratios = df_with_minimum_cutoff/(df_end_counts.reindex(df_with_minimum_cutoff.swaplevel('start', 'end').index).swaplevel('start', 'end').fillna(0) + df_with_minimum_cutoff)

df_s_start_counts = df_s_with_minimum_cutoff.query('is_canonical == True').groupby(['chr', 'start']).sum().reindex(df_s_with_minimum_cutoff.index).fillna(0)
df_s_end_counts = df_s_with_minimum_cutoff.query('is_canonical == True').groupby(['chr', 'end']).sum().reindex(df_s_with_minimum_cutoff.swaplevel('start', 'end').index).swaplevel('start', 'end').fillna(0)
df_s_start_ratios = df_s_with_minimum_cutoff/(df_s_start_counts + df_s_with_minimum_cutoff)
df_s_end_ratios = df_s_with_minimum_cutoff/(df_s_end_counts + df_s_with_minimum_cutoff)

print('Resource checkpoint 3:', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
del df_start_counts
del df_s_start_counts
gc.collect()

# make tables for start of junction
df_start_family_ratios = df_start_ratios.groupby(lambda x: get_proband_id(x), axis=1).max().fillna(0)
df_start_mean_family_ratio = df_start_family_ratios.mean(axis=1).rename('mean_family_max_ratio')

# get top family ratios
arr = np.argsort(-df_start_family_ratios.values, axis=1)
df_start_family_sorted_column_names = pd.DataFrame(df_start_family_ratios.columns[arr], index=df_start_family_ratios.index)
df_start_family_sorted_ratios = np.sort(-df_start_family_ratios.values, axis=1)
n_top = 4
df_start_family_top_ratios = df_start_family_sorted_column_names.loc[:, list(range(n_top))].copy()
df_start_family_top_ratios.columns = ['top_family_' + str(i) for i in range(n_top)]
for i in range(n_top):
  df_start_family_top_ratios['top_family_ratio_' + str(i)] = -df_start_family_sorted_ratios[:, i]

# get top sample ratios -- this is not as important but it's a kludge to get the number of reads
arr = np.argsort(-df_start_ratios.values, axis=1)
df_start_sorted_column_names = pd.DataFrame(df_start_ratios.columns[arr], index=df_start_ratios.index)
df_start_sorted_ratios = np.sort(-df_start_ratios.values, axis=1)
n_top = 4
df_start_top_ratios = df_start_sorted_column_names.loc[:, list(range(n_top))].copy()
df_start_top_ratios.columns = ['top_sample_' + str(i) for i in range(n_top)]
for i in range(n_top):
  df_start_top_ratios['top_sample_ratio_' + str(i)] = -df_start_sorted_ratios[:, i]

# get reads for top sample ratios
df_start_joined = df_s_start_ratios.rename(columns={'n_unique': 'proband_ratio'}).join(df_s_with_minimum_cutoff.rename(columns={'n_unique': 'proband_reads'})).join(df_start_mean_family_ratio)
df_start_joined['n_family_ratios_ge_0_01'] = (df_start_family_ratios >= 0.01).sum(axis=1)
df_start_joined['n_family_ratios_ge_0_1'] = (df_start_family_ratios >= 0.1).sum(axis=1)
df_start_joined['n_family_ratios_ge_0_2'] = (df_start_family_ratios >= 0.2).sum(axis=1)
df_start_joined['n_family_ratios_ge_0_3'] = (df_start_family_ratios >= 0.3).sum(axis=1)
df_start_joined = df_start_joined.join(df_start_family_top_ratios).join(df_start_top_ratios)
df_start_joined_with_read_count = df_start_joined.join(df_with_minimum_cutoff)
for i in range(n_top):
  # this is a weird hack to fill na values so that the column lookup works
  df_start_joined_with_read_count['top_sample_' + str(i)] = df_start_joined_with_read_count['top_sample_' + str(i)].fillna('top_sample_' + str(i))
  df_start_joined['top_sample_reads_' + str(i)] = df_start_joined_with_read_count.apply(lambda r: r[r['top_sample_' + str(i)]], axis=1)
df_start_joined['side'] = 'start'

# make tables for end of junction
df_end_family_ratios = df_end_ratios.groupby(lambda x: get_proband_id(x), axis=1).max().fillna(0)
df_end_mean_family_ratio = df_end_family_ratios.mean(axis=1).rename('mean_family_max_ratio')

# get top family ratios
arr = np.argsort(-df_end_family_ratios.values, axis=1)
df_end_family_sorted_column_names = pd.DataFrame(df_end_family_ratios.columns[arr], index=df_end_family_ratios.index)
df_end_family_sorted_ratios = np.sort(-df_end_family_ratios.values, axis=1)
n_top = 4
df_end_family_top_ratios = df_end_family_sorted_column_names.loc[:, list(range(n_top))].copy()
df_end_family_top_ratios.columns = ['top_family_' + str(i) for i in range(n_top)]
for i in range(n_top):
  df_end_family_top_ratios['top_family_ratio_' + str(i)] = -df_end_family_sorted_ratios[:, i]

# get top sample ratios -- this is not as important but it's a kludge to get the number of reads
arr = np.argsort(-df_end_ratios.values, axis=1)
df_end_sorted_column_names = pd.DataFrame(df_end_ratios.columns[arr], index=df_end_ratios.index)
df_end_sorted_ratios = np.sort(-df_end_ratios.values, axis=1)
n_top = 4
df_end_top_ratios = df_end_sorted_column_names.loc[:, list(range(n_top))].copy()
df_end_top_ratios.columns = ['top_sample_' + str(i) for i in range(n_top)]
for i in range(n_top):
  df_end_top_ratios['top_sample_ratio_' + str(i)] = -df_end_sorted_ratios[:, i]

# get reads for top sample ratios
df_end_joined = df_s_end_ratios.rename(columns={'n_unique': 'proband_ratio'}).join(df_s_with_minimum_cutoff.rename(columns={'n_unique': 'proband_reads'})).join(df_end_mean_family_ratio)
df_end_joined['n_family_ratios_ge_0_01'] = (df_end_family_ratios >= 0.01).sum(axis=1)
df_end_joined['n_family_ratios_ge_0_1'] = (df_end_family_ratios >= 0.1).sum(axis=1)
df_end_joined['n_family_ratios_ge_0_2'] = (df_end_family_ratios >= 0.2).sum(axis=1)
df_end_joined['n_family_ratios_ge_0_3'] = (df_end_family_ratios >= 0.3).sum(axis=1)
df_end_joined = df_end_joined.join(df_end_family_top_ratios).join(df_end_top_ratios)
df_end_joined_with_read_count = df_end_joined.join(df_with_minimum_cutoff)
for i in range(n_top):
  # this is a weird hack to fill na values so that the column lookup works
  df_end_joined_with_read_count['top_sample_' + str(i)] = df_end_joined_with_read_count['top_sample_' + str(i)].fillna('top_sample_' + str(i))
  df_end_joined['top_sample_reads_' + str(i)] = df_end_joined_with_read_count.apply(lambda r: r[r['top_sample_' + str(i)]], axis=1)
df_end_joined['side'] = 'end'

# assemble data frames to write
df_concat = pd.concat([df_start_joined, df_end_joined], sort=False)
df_concat.to_csv(args.output_dir.rstrip('/') + '/' + sample_name + '_rare_junctions_chr_' + args.chromosome + '.tsv', sep='\t')

# TO DO:
# remove family from controls?
# implement per tissue check?
# include other methods
