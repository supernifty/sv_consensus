#!/usr/bin/env python
'''
  counts concordance between vcf files

  Usage:
    python snp_concordance.py file(s)

  Output:
    lists the number of variants with each level of concordance
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
      if record.POS not in vcf_calls[record.CHROM]:
        vcf_calls[record.CHROM][record.POS] = 0 # counter
      vcf_calls[record.CHROM][record.POS] += 1

      # remember last for dupes
      last = (record.CHROM, record.POS)
    logging.info('processing {}: {} records. {} duplicates.'.format(filename, line, duplicates))

  # now count concordance
  counts = collections.defaultdict(int)
  for chrom in vcf_calls:
    for pos in vcf_calls[chrom]:
      counts[vcf_calls[chrom][pos]] += 1

  sys.stdout.write('Concordance,VariantCount\n')
  sys.stdout.write('\n'.join(['{},{}'.format(concordance, counts[concordance]) for concordance in sorted(counts.keys())]))

def main():
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  parser = argparse.ArgumentParser(description='SNV consensus')
  parser.add_argument('vcfs', nargs='+', help='input vcfs')
  args = parser.parse_args()
  calculate(args.vcfs)

if __name__ == '__main__':
  main()
