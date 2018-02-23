#!/usr/bin/env python
'''
  combine intersect files into one big intersect

  Usage:
    python snv_intersect.py --vcf_names caller1 [callers] --intersects file(s)
'''

import argparse
import collections
import itertools
import logging
import sys
import vcf

import matplotlib
matplotlib.use('Agg')
import venn

def calculate(venns, venn_img, venn_txt, names, max_variants=None):
  total = collections.defaultdict(int)
  grand = 0
  for f in venns:
    lines = [line.strip('\n').split(',') for line in open(f, 'r').readlines()]
    file_total = sum([int(x[1]) for x in lines])
    if max_variants is not None and file_total > max_variants:
      logging.info('skipping %s with %i > %i variants', f, file_total, max_variants)
      continue
    for line in lines:
      total[line[0]] += int(line[1])
      grand += int(line[1])
    logging.info('added %s with variant count %i for total %i', f, file_total, grand)
  
  # venn diagram
  if venn_img is not None and 2 <= len(names) <= 6:
    logging.info('writing venn diagram to %s...', venn_img)
    venn_fn = (None, None, venn.venn2, venn.venn3, venn.venn4, venn.venn5, venn.venn6)
    fig, ax = venn_fn[len(names)](total, names=names)
    matplotlib.pyplot.savefig(venn_img)
    logging.info('writing venn diagram: done')
  else:
    logging.info('skipping venn diagram generation')

  if venn_txt is not None:
    logging.info('writing venn text to %s...', venn_txt)
    open(venn_txt, 'w').write('\n'.join(['{},{}'.format(k, total[k]) for k in sorted(total)]))
    logging.info('writing venn text: done')

def main():
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  parser = argparse.ArgumentParser(description='SNV consensus')
  parser.add_argument('--venn', required=False, help='filename to write venn diagram to')
  parser.add_argument('--venn_txt', required=False, help='filename to write venn intersects to')
  parser.add_argument('--vcf_names', required=True, nargs='+', help='original vcf names')
  parser.add_argument('--max_variants', required=False, type=int, default=1e12, help='original vcf names')
  parser.add_argument('--intersects', nargs='+', help='intersect files')
  args = parser.parse_args()
  calculate(args.intersects, args.venn, args.venn_txt, args.vcf_names, args.max_variants)

if __name__ == '__main__':
  main()
