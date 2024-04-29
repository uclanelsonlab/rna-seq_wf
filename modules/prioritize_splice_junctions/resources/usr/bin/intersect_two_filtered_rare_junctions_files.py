#!/usr/bin/env python

import argparse
import pandas as pd
from os.path import basename

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
parser.add_argument('--input_1', required=True, help='xlsx 1')
parser.add_argument('--input_2', required=True, help='xlsx 2')
parser.add_argument('--output_dir', required=True, help='output dir')
# parser.add_argument('--output_dir', required=True, help='output dir')
args = parser.parse_args()

input_prefix_1 = basename(args.input_1).replace('_rare_junctions_filtered.xlsx', '')
input_prefix_2 = basename(args.input_2).replace('_rare_junctions_filtered.xlsx', '')
output_filename = input_prefix_1 + "_" + input_prefix_2 + "_rare_junctions_intersected.xlsx"

number_of_sheets = min(len(pd.ExcelFile(args.input_1).sheet_names), len(pd.ExcelFile(args.input_1).sheet_names))
output_list = []
for sheet_number in range(number_of_sheets):
  # If Excel has numbers stored as text, Pandas tries to respect that and causes problems. So cast all chrs to str.
  input_1 = pd.read_excel(args.input_1, sheet_name=sheet_number)
  input_1['chr'] = input_1['chr'].astype('str')
  input_1 = input_1.set_index(['chr', 'start', 'end', 'is_canonical', 'is_canonical_start', 'is_canonical_end'])

  input_2 = pd.read_excel(args.input_2, sheet_name=sheet_number)
  input_2['chr'] = input_2['chr'].astype('str')
  input_2 = input_2.set_index(['chr', 'start', 'end', 'is_canonical', 'is_canonical_start', 'is_canonical_end'])

  joined_index = input_1.index.intersection(input_2.index)

  if joined_index.size > 0:
    df_in = pd.concat([input_1.loc[joined_index], input_2.loc[joined_index]])
  
    # sort (borrowed from postprocess.py)
    df_in['sort_value'] = (df_in['omim_disorders_inheritance'].fillna('').str.replace(',', '') == '')*10**19 + \
    (df_in['gene'].fillna('').str.replace(',', '') == '')*10**18 + \
    df_in['n_family_ratios_ge_0_3']*10**15 + \
    df_in['n_family_ratios_ge_0_2']*10**12 + \
    df_in['n_family_ratios_ge_0_1']*10**9 + \
    df_in['n_family_ratios_ge_0_01']*10**6 - \
    df_in['proband_ratio']*10**3 - \
    df_in['proband_reads']

    output_list.append(df_in.copy())
  else:
    output_list.append(None)

# make names if the inputs are inconsistent
output_sheet_names = ['at_least_one_side_canonical', 'private_to_family']

# if output_list != [None, None]:
# write data frames
with pd.ExcelWriter(args.output_dir.rstrip('/') + '/' + output_filename) as writer:
  for sheet_number in range(number_of_sheets):
    if output_list[sheet_number] is not None:
      output_list[sheet_number].sort_values(by=['sort_value']).drop('sort_value', 1).to_excel(writer, merge_cells=False, sheet_name=output_sheet_names[sheet_number])
