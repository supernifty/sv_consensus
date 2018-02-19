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

def calculate(vcfs):
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
  for set_size in range(1, len(vcfs) + 1):
    for vcf_set in itertools.combinations(vcfs, set_size):
      # determine intersect of all of these sets, for each chromosome
      total_intersect = 0

      for chrom in vcf_calls:
        if who[vcf_set[0]] not in vcf_calls[chrom]:
          #logging.info('{} has no variants on chromosome {}'.format(who[vcf_set[0]], chrom))
          continue
        chrom_intersect = vcf_calls[chrom][who[vcf_set[0]]]
        if len(vcf_set) > 1:
          for intersect_id in range(1, len(vcf_set)):
            if who[vcf_set[intersect_id]] not in vcf_calls[chrom]:
              #logging.info('{} has no variants on chromosome {}'.format(who[vcf_set[0]], chrom))
              chrom_intersect = set()
              continue
            chrom_intersect.intersection_update(vcf_calls[chrom][who[vcf_set[intersect_id]]])
        total_intersect += len(chrom_intersect)

      sys.stdout.write('{},{}\n'.format(' '.join(vcf_set), total_intersect))

def main():
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  parser = argparse.ArgumentParser(description='SNV consensus')
  parser.add_argument('vcfs', nargs='+', help='input vcfs')
  args = parser.parse_args()
  calculate(args.vcfs)

if __name__ == '__main__':
  main()
