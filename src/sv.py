#!/usr/bin/env python

import argparse
import collections
import logging
import sys
import vcf

import blist

def write_cluster_header():
    sys.stdout.write('chrom,size,consensus,contents\n')

def write_cluster(chrom, cluster, vcf_calls):
  if len(cluster) > 0:
    # find all callers in the cluster
    callers = set()
    for pos in cluster:
      for vcf_id, vcf_call in enumerate(vcf_calls):
        if pos in vcf_call[chrom]:
          callers.add(vcf_id)
    sys.stdout.write('{},{},{},{}\n'.format(chrom, max(cluster) - min(cluster), len(callers), ' '.join([str(x) for x in cluster])))

def calculate(vcfs, tolerance):
  all_calls = collections.defaultdict(blist.sortedlist)
  vcf_calls = [collections.defaultdict(set) for _ in range(len(vcfs))]
  who = collections.defaultdict(dict)
  for vcf_id, filename in enumerate(vcfs):
    logging.info('processing {}...'.format(filename))
    reader = vcf.Reader(open(filename, 'r'))
    line = 0
    last = (None, None)
    duplicates = 0
    for line, record in enumerate(reader):
      # has this caller already called this position?
      if record.CHROM == last[0] and record.POS == last[1]:
        duplicates += 1
        continue
      # add to list
      all_calls[record.CHROM].add(record.POS)
      vcf_calls[vcf_id][record.CHROM].add(record.POS)

      # remember last for dupes
      last = (record.CHROM, record.POS)
    logging.info('processing {}: {} records. {} duplicates'.format(filename, line, duplicates))

  # write clusters
  write_cluster_header()

  for chrom in all_calls:
    current = None
    cluster = []
    for pos in all_calls[chrom]:
      if current is not None:
        if pos - current > tolerance: # new cluster
          write_cluster(chrom, cluster, vcf_calls) # write old cluster
          cluster = [] # start new cluster
        else: # add to existing
          pass
 
      cluster.append(pos)
      current = pos

    # write what's left
    write_cluster(chrom, cluster, vcf_calls)

def main():
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  parser = argparse.ArgumentParser(description='SV consensus')
  parser.add_argument('--tolerance', type=int, default=50, help='maximum distance to still consider the same breakpoint')
  parser.add_argument('vcfs', nargs='+', help='input vcfs')
  args = parser.parse_args()
  calculate(args.vcfs, args.tolerance)

if __name__ == '__main__':
  main()
