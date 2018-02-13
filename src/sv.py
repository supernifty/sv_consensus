#!/usr/bin/env python

import argparse
import collections
import logging
import sys
import vcf

import blist

def write_cluster(chrom, cluster):
  if len(cluster) > 0:
    sys.stdout.write('{}: {}\n'.format(chrom, ' '.join([str(x) for x in cluster])))

def calculate(vcfs, tolerance):
  chroms = collections.defaultdict(blist.sortedlist)
  for filename in vcfs:
    logging.info('processing {}...'.format(filename))
    reader = vcf.Reader(open(filename, 'r'))
    line = 0
    for line, record in enumerate(reader):
      chroms[record.CHROM].add(record.POS)
    logging.info('processing {}: {} records'.format(filename, line))

  # write clusters
  for chrom in chroms:
    current = None
    cluster = []
    for pos in chroms[chrom]:
      if current is not None:
        if pos - current > tolerance: # new cluster
          write_cluster(chrom, cluster) # write old cluster
          cluster = [] # start new cluster
        else: # add to existing
          pass
 
      cluster.append(pos)
      current = pos

    # write what's left
    write_cluster(chrom, cluster)

def main():
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  parser = argparse.ArgumentParser(description='SV consensus')
  parser.add_argument('--tolerance', type=int, default=50, help='maximum distance to still consider the same breakpoint')
  parser.add_argument('vcfs', nargs='+', help='input vcfs')
  args = parser.parse_args()
  calculate(args.vcfs, args.tolerance)

if __name__ == '__main__':
  main()
