#!/usr/bin/env python
'''
  Usage:
    python combine_concordance.py file(s)
'''

import argparse
import collections
import itertools
import logging
import sys
import vcf

def calculate(files):
  first = True
  minimum = 99
  maximum = -1
  for f in files:
    values = collections.defaultdict(int)
    with open(f, 'r') as fh:
      for line in fh.readlines()[1:]:
        fields = line.strip('\n').split(',')
        values[int(fields[0])] = fields[1]
        if first:
          minimum = min(minimum, int(fields[0]))
          maximum = max(maximum, int(fields[0]))
    if first:
      sys.stdout.write('Filename,{}\n'.format(','.join([str(x) for x in range(minimum, maximum+1)])))
    sys.stdout.write('{},{}\n'.format(f, ','.join([str(values[x]) for x in range(minimum, maximum+1)])))
    first = False

def main():
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  parser = argparse.ArgumentParser(description='combine concordance files into one file')
  parser.add_argument('files', nargs='+', help='concordance files')
  args = parser.parse_args()
  calculate(args.files)

if __name__ == '__main__':
  main()
