#!/usr/bin/env python
'''
  write out a vcf with variants that have the specified minimum level of concordance between vcfs

  Usage:
    python snv_concordance_vcf.py [--minimum_concordance n] file(s)
'''

import argparse
import collections
import itertools
import logging
import sys
import vcf

def calculate(vcfs, minimum_concordance):
  vcf_calls = collections.defaultdict(dict)
  who = {}
  written_total = 0
  writer = None
  
  for vcf_id, filename in enumerate(vcfs):
    logging.info('processing {}...'.format(filename))
    reader = vcf.Reader(open(filename, 'r'))
    if writer is None:
      writer = vcf.Writer(sys.stdout, reader)
    line = 0
    last = (None, None)
    duplicates = 0
    written = 0
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
      if vcf_calls[record.CHROM][record.POS] == minimum_concordance: # hit minimum for the first time
        writer.write_record(record)
        written += 1

      # remember last for dupes
      last = (record.CHROM, record.POS)
    written_total += written
    logging.info('processing {}: {} records. {} duplicates. wrote {}'.format(filename, line, duplicates, written))
  logging.info('wrote {} total'.format(written_total))

def main():
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  parser = argparse.ArgumentParser(description='SNV consensus')
  parser.add_argument('--minimum_concordance', type=int, default=2, required=False, help='input vcfs')
  parser.add_argument('vcfs', nargs='+', help='input vcfs')
  args = parser.parse_args()
  calculate(args.vcfs, args.minimum_concordance)

if __name__ == '__main__':
  main()
