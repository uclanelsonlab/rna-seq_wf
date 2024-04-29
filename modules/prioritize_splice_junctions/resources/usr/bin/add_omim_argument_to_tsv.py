#!/usr/bin/env python

import argparse
from os.path import basename

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
parser.add_argument('--input_filename', required=True, help='input tsv')
parser.add_argument('--omim_argument', required=True, help='omim')
parser.add_argument('--output_dir', required=True, help='output dir')
args = parser.parse_args()

output_filename = basename(args.input_filename)
if output_filename.endswith('.tsv') or output_filename.endswith('.bed'):
  output_filename = output_filename[:-len('.tsv')]
output_filename += '_with_omim.tsv'

omim_list = []

with open(args.omim_argument) as f:
  for l in f:
    chrom, start, end, gene_name = l.strip('\n').split('\t')[:4]
    if chrom == 'Chromosome':
      continue
    start = int(start)
    end = int(end)
    disorders_inheritance = l.strip('\n').split('\t')[20]
    omim_list.append((chrom, start, end, gene_name, disorders_inheritance))

c = 0
with open(args.output_dir.rstrip('/') + '/' + output_filename, 'w') as g:
  with open(args.input_filename) as f:
    for l in f:
      line_list = l.strip('\n').split('\t')
      chrom = line_list[0]
      if l.startswith('#') or chrom in ('chr', 'chrom'):
        line_list.append('omim_genes')
        line_list.append('omim_disorders_inheritance')
        g.write('\t'.join(line_list) + '\n')

      else:
        start = int(line_list[1])
        end = int(line_list[2])
        omim_gene_list = []
        omim_disorders_inheritance_list = []
        for omim_chrom, omim_start, omim_end, omim_gene_name, omim_disorders_inheritance in omim_list:
          if chrom == omim_chrom:
            if int(start) <= int(omim_end):
              if int(end) >= int(omim_start):
                  omim_gene_list.append(omim_gene_name)
                  omim_disorders_inheritance_list.append(omim_disorders_inheritance)
        line_list.append(','.join(omim_gene_list))
        line_list.append(','.join(omim_disorders_inheritance_list))
        g.write('\t'.join(line_list) + '\n')

        c += 1
        if c % 100 == 0:
          print(c)
