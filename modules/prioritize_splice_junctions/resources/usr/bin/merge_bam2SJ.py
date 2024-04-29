#!/usr/bin/env python

from collections import Counter

filenames = ['UDN541948-P-muscle.bam2SJ.out.tab', 'UDN541948-P-muscle-2.bam2SJ.out.tab', 'UDN541948-P-muscle-3.bam2SJ.out.tab', 'UDN541948-P-muscle-4.bam2SJ.out.tab']

c = Counter()
for filename in filenames:
  with open(filename) as f:
    for l in f:
      chrom, start_str, end_str, count_str = l.strip('\n').split('\t')
      start = int(start_str)
      end = int(end_str)
      count = int(count_str)
      c[(chrom, start, end)] += count

with open('UDN541948-P-muscle-1-4.merged_manual.bam2SJ.out.tab', 'w') as g:
  for x in c:
    chrom, start, end = x
    start_str = str(start)
    end_str = str(end)
    count_str = str(c[x])
    g.write('\t'.join([chrom, start_str, end_str, count_str]) + '\n')
