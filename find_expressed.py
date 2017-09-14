#!/usr/bin/env python

import argparse
import collections
import logging
import operator
import sys

import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot

def generate(files, replicates):

  counts = collections.defaultdict(dict)
  time_period = 0

  for i, f in enumerate(files):
    logging.info('processing %s', f)
    if i % replicates == 0:
      time_period += 1

    replicate = i % replicates
  
    for line in open(f, 'r'):
      if line.startswith('#') or line.startswith('Geneid'):
        continue
      fields = line.strip().split('\t')
      if len(fields) < 7:
        continue

      gene = fields[0]
      count = int(fields[6])
      counts[gene]['{},{}'.format(time_period, replicate)] = count

  logging.info('calculating folds')
  folds = {}
  found = 0
  total = 0
  for total, gene in enumerate(counts):
    periods = collections.defaultdict(int)
    for key in counts[gene].keys():
      time_period, replicate = key.split(',')
      periods[int(time_period)] += counts[gene][key]

    #denominator = periods[3] + periods[4] 
    #numerator = periods[1] + periods[6] 
    denominator = periods[4] 
    numerator = periods[2]

    if denominator > 0:
      fold = numerator / denominator
      folds[gene] = fold
      found += 1
    else:
      #logging.debug('zero counts for %s', gene)
      pass
  
  # print sorted folds
  logging.info('writing folds (calculated %i of %i)', found, total)
  for k, v in reversed(sorted(folds.items(), key=operator.itemgetter(1))):
    sys.stdout.write('{},{:.2f}\n'.format(k, v))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='plot gene expression')
  parser.add_argument('--replicates', type=int, default=3, required=True, help='files per group')
  parser.add_argument('--files', nargs='+', help='files containing counts')
  args = parser.parse_args()
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)

  generate(args.files, args.replicates)
