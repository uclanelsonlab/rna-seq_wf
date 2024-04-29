#!/usr/bin/env python

import argparse
import pandas as pd
from collections import defaultdict
from os.path import basename

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
parser.add_argument('--input_tsv', required=True, help='tsv')
parser.add_argument('--output_filename', required=True, help='output_filename')
parser.add_argument('--gencode', required=True, help='gencode exons cds only')
args = parser.parse_args()

output_filename = basename(args.input_tsv)
gencode = args.gencode
if output_filename.endswith('.tsv'):
  output_filename = output_filename[:-len('.tsv')]
output_filename += '.xlsx'

# We already pre-filtered using awk because it vastly decreases the memory requirement for Python.

# read input
df_all = pd.read_table(args.input_tsv, index_col=list(range(6)))

# Sheet 1: one_side_canonical
df_at_least_one_side_canonical = df_all.query("(is_canonical_start == True and side == 'start') or (is_canonical_end == True and side == 'end')").copy()
df_at_least_one_side_canonical['sort_value'] = (df_at_least_one_side_canonical['omim_disorders_inheritance'].fillna('').str.replace(',', '') == '')*10**19 + \
  (df_at_least_one_side_canonical['gene'].fillna('').str.replace(',', '') == '')*10**18 + \
  df_at_least_one_side_canonical['n_family_ratios_ge_0_3']*10**15 + \
  df_at_least_one_side_canonical['n_family_ratios_ge_0_2']*10**12 + \
  df_at_least_one_side_canonical['n_family_ratios_ge_0_1']*10**9 + \
  df_at_least_one_side_canonical['n_family_ratios_ge_0_01']*10**6 - \
  df_at_least_one_side_canonical['proband_ratio']*10**3 - \
  df_at_least_one_side_canonical['proband_reads']

# Sheet 2: private_to_family
df_private_to_family = df_all.query('top_family_ratio_1 == 0 and is_canonical_start == False and is_canonical_end == False').copy()

exon_intervals = defaultdict(list)
with open(gencode) as f:
  for l in f:
    chrom, start, end = l.strip('\n').split('\t')[:3]
    start, end = int(start), int(end)
    exon_intervals[chrom].append((start, end))

endpoint_to_exon_list = []
spans_cds_list = []
for index, row in df_private_to_family.iterrows():
  chrom = str(index[0])
  row_start = index[1]
  row_end = index[2]
  min_d = 99999
  spans_cds = False
  for interval_start, interval_end in exon_intervals[chrom]:
    d = min(abs(row_end - interval_start), abs(row_start - interval_end))
    min_d = min(d, min_d)
    if row_start <= interval_end and row_end >= interval_start:
      spans_cds = True
  endpoint_to_exon_list.append(min_d)
  spans_cds_list.append(spans_cds)

df_private_to_family['endpoint_to_exon'] = endpoint_to_exon_list
df_private_to_family['spans_cds'] = spans_cds_list

df_private_to_family['sort_value'] = df_private_to_family['endpoint_to_exon']*10**20 + \
  (df_private_to_family['omim_disorders_inheritance'].fillna('').str.replace(',', '') == '')*10**19 + \
  (df_private_to_family['gene'].fillna('').str.replace(',', '') == '')*10**18 + \
  df_private_to_family['n_family_ratios_ge_0_3']*10**15 + \
  df_private_to_family['n_family_ratios_ge_0_2']*10**12 + \
  df_private_to_family['n_family_ratios_ge_0_1']*10**9 + \
  df_private_to_family['n_family_ratios_ge_0_01']*10**6 - \
  df_private_to_family['proband_ratio']*10**3 - \
  df_private_to_family['proband_reads']

# write data frames
with pd.ExcelWriter(args.output_filename) as writer:
  df_at_least_one_side_canonical.sort_values(by=['sort_value']).drop('sort_value', 1).to_excel(writer, merge_cells=False, sheet_name='at_least_one_side_canonical')
  df_private_to_family.sort_values(by=['sort_value']).drop('sort_value', 1).to_excel(writer, merge_cells=False, sheet_name='private_to_family')
