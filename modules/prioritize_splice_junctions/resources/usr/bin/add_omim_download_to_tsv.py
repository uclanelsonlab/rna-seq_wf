#!/usr/bin/env python

import argparse
from os.path import basename

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
parser.add_argument('--input_filename', required=True, help='input tsv')
parser.add_argument('--omim_genemap2_lifted_hg19', required=True, help='OMIM genemap path')
parser.add_argument('--output_dir', required=True, help='output dir')
args = parser.parse_args()

output_filename = basename(args.input_filename)
if output_filename.endswith('.tsv') or output_filename.endswith('.bed'):
  output_filename = output_filename[:-len('.tsv')]
output_filename += '_with_omim.tsv'

# Make OMIM list
omim_position_list = []
with open(args.omim_genemap2_lifted_hg19) as f:
  for l in f:
    if l.startswith('#'):
      continue
    line_list = l.strip('\n').split('\t')

    chrom, start, end = line_list[:3]
    if chrom.startswith('chr'):
      chrom = chrom[3:]
    start = int(start)
    end = int(end)

    gene_hgnc_symbol = line_list[11]
    phenotype_list = line_list[15].split(';')
    omim_disorders_inheritance_list = []

    for p in phenotype_list:
      if p == '':
        continue

      disorder_mim = None
      for x in p.split(','):
        if len(x.strip()) == 10 and x.strip()[:6].isdigit():
          if disorder_mim:
            assert False, disorder_mim
          disorder_mim = x.strip()

      if disorder_mim:
        disorder_title = p.split(disorder_mim)[0].strip(', ')
        omim_disorders_inheritance = p.split(disorder_mim)[1].strip(', ')
      else:
        disorder_title = p.strip(', ')
        disorder_mim = ''
        omim_disorders_inheritance = ''

      if omim_disorders_inheritance:
        omim_disorders_inheritance_list.append(omim_disorders_inheritance)

    omim_position_list.append([chrom, start, end, gene_hgnc_symbol, omim_disorders_inheritance_list])

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
        omim_gene_set = set()
        omim_disorders_inheritance_set = set()
        for omim_chrom, omim_start, omim_end, omim_gene_name, omim_disorders_inheritance_list in omim_position_list:
          if chrom == omim_chrom:
            if int(start) <= int(omim_end):
              if int(end) >= int(omim_start):
                  omim_gene_set.add(omim_gene_name)
                  omim_disorders_inheritance_set.update(omim_disorders_inheritance_list)
        line_list.append(','.join(sorted(omim_gene_set)))
        line_list.append(','.join(sorted(omim_disorders_inheritance_set)))
        g.write('\t'.join(line_list) + '\n')

        c += 1
        if c % 100 == 0:
          print(c)
