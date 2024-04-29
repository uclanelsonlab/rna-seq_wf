#!/usr/bin/env python

import argparse
from os.path import basename
from collections import defaultdict

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
parser.add_argument('--input_filename', required=True, help='input tsv')
parser.add_argument('--output_dir', required=True, help='output dir')
parser.add_argument('--constraint', required=True, help='constraint with interval canonical')
args = parser.parse_args()

output_filename = basename(args.input_filename)
constraint = args.constraint
if output_filename.endswith('.tsv'):
  output_filename = output_filename[:-len('.tsv')]
output_filename += '_with_constraints.tsv'

constraint_dict = defaultdict(list)
with open(constraint) as f:
  for l in f:
    if l.startswith('gene'):
      continue
    line_list = l.strip('\n').split('\t')
    gene, transcript, canonical, obs_lof, exp_lof, oe_lof, oe_lof_lower, oe_lof_upper, obs_mis, exp_mis, oe_mis, oe_mis_lower, oe_mis_upper, obs_syn, exp_syn, oe_syn, oe_syn_lower, oe_syn_upper, lof_z, mis_z, syn_z, pLI, pRec, pNull, gene_issues, chrom, start, end = line_list
    start = int(start)
    end = int(end)
    constraint_dict[chrom].append((start, end, gene, pLI, pRec, mis_z))

def max_with_NA(a, b):
  if 'NA' in (a, b):
    return 'NA'
  return max(float(a), float(b))

c = 0
with open(args.output_dir.rstrip('/') + '/' + output_filename, 'w') as g:
  with open(args.input_filename) as f:
    for l in f:
      line_list = l.strip('\n').split('\t')
      chrom = line_list[0]
      if chrom in ('chr', 'chrom'):
        line_list.extend(['gene', 'pLI', 'pRec', 'mis_z'])
        g.write('\t'.join(line_list) + '\n')

      else:
        line_start = int(line_list[1])
        line_end = int(line_list[2])

        gene_list = []
        pLI_max, pRec_max, mis_z_max = -9, -9, -9

        # if there are n transcripts, make 1 output line
        for start, end, gene, pLI, pRec, mis_z in constraint_dict[chrom]:
          if int(line_start) <= end:
            if int(line_end) >= start:
              gene_list.append(gene)
              pLI_max = max_with_NA(pLI_max, pLI)
              pRec_max = max_with_NA(pRec_max, pRec)
              mis_z_max = max_with_NA(mis_z_max, mis_z)

        if gene_list:
          line_list += [','.join(gene_list), str(pLI_max), str(pRec_max), str(mis_z_max)]
        else:
          line_list += [''] * 4
        g.write('\t'.join(line_list) + '\n')

      c += 1
      if c % 1000 == 0:
        print(c)
