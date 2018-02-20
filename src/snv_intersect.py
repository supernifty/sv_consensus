#!/usr/bin/env python
'''
  calculates the number of variants that intersect across every combination of vcf files

  Usage:
    python snv_intersect.py file(s)
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

def calculate(vcfs, venn_img):
  vcf_calls = collections.defaultdict(dict)
  who = {}
  for vcf_id, filename in enumerate(vcfs):
    logging.info('processing {}...'.format(filename))
    reader = vcf.Reader(open(filename, 'r'))
    line = 0
    last = (None, None)
    duplicates = 0
    who[filename] = vcf_id
    for line, record in enumerate(reader):
      # has this caller already called this position?
      if record.CHROM == last[0] and record.POS == last[1]:
        duplicates += 1
        continue
      # add to list
      if vcf_id not in vcf_calls[record.CHROM]:
        vcf_calls[record.CHROM][vcf_id] = set()
      vcf_calls[record.CHROM][vcf_id].add(record.POS)

      # remember last for dupes
      last = (record.CHROM, record.POS)
    logging.info('processing {}: {} records. {} duplicates'.format(filename, line, duplicates))

  sys.stdout.write('Sets,Count\n')
  labels = {}
  for set_size in range(1, len(vcfs) + 1):
    for vcf_set in itertools.combinations(vcfs, set_size): # AB AC BC
      # determine intersect of all of these sets, for each chromosome
      total_intersect = 0
      total_difference = 0

      for chrom in vcf_calls:
        if who[vcf_set[0]] not in vcf_calls[chrom]: # chromosome not in first set
          #logging.info('{} has no variants on chromosome {}'.format(who[vcf_set[0]], chrom))
          continue

        # take the intersection of all selected vcfs
        chrom_intersect = vcf_calls[chrom][who[vcf_set[0]]].copy()

        if len(vcf_set) > 1: # additional sets in list
          for intersect_id in range(1, len(vcf_set)):
            if who[vcf_set[intersect_id]] not in vcf_calls[chrom]:
              #logging.info('{} has no variants on chromosome {}'.format(who[vcf_set[0]], chrom))
              chrom_intersect = set()
              continue
            chrom_intersect.intersection_update(vcf_calls[chrom][who[vcf_set[intersect_id]]])

        total_intersect += len(chrom_intersect)

        # remove all non-selected sets
        for excluded in vcfs:
          if excluded not in vcf_set:
            if who[excluded] not in vcf_calls[chrom]:
              continue
            chrom_intersect.difference_update(vcf_calls[chrom][who[excluded]])
        
        total_difference += len(chrom_intersect)

      sys.stdout.write('{},{}\n'.format(' '.join(vcf_set), total_intersect))

      label = []
      for i in vcfs:
        if i in vcf_set:
          label.append('1')
        else:
          label.append('0')
      labels[''.join(label)] = total_difference

  # venn diagram
  if venn_img is not None and 2 <= len(vcfs) <= 6:
    logging.info('writing venn diagram to %s...', venn_img)
    venn_fn = (None, None, venn.venn2, venn.venn3, venn.venn4, venn.venn5, venn.venn6)
    fig, ax = venn_fn[len(vcfs)](labels, names=vcfs)
    matplotlib.pyplot.savefig(venn_img)
    logging.info('writing venn diagram: done')
  else:
    logging.info('skipping venn diagram generation')

def main():
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  parser = argparse.ArgumentParser(description='SNV consensus')
  parser.add_argument('--venn', required=False, help='filename to write venn diagram to')
  parser.add_argument('vcfs', nargs='+', help='input vcfs')
  args = parser.parse_args()
  calculate(args.vcfs, args.venn)

if __name__ == '__main__':
  main()
