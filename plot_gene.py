#!/usr/bin/env python

import argparse
import collections
import logging

import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot

def generate(gene, files, replicates, output):

  xs = []
  ys = collections.defaultdict(list)

  time_period = 0

  for i, f in enumerate(files):
    logging.info('processing %s', f)
    if i % replicates == 0:
      time_period += 1
      xs.append(time_period)

    replicate = i % replicates
  
    total = 0
    for line in open(f, 'r'):
      if line.startswith('#') or line.startswith('Geneid'):
        continue
      fields = line.strip().split('\t')
      if len(fields) < 7:
        continue
      current = fields[0]
      count = int(fields[6])
      total += count
      if current == gene:
        ys[replicate].append(count)
        logging.debug('added value %i to time %i replicate %i', count, time_period, replicate)

    if gene == 'total':
      ys[replicate].append(total)
  
  # Create traces
  traces = []
  for replicate in ys.keys():
    traces.append(go.Scatter(
      x = xs,
      y = ys[replicate],
      mode = 'lines',
      name = 'Replicate {}'.format(replicate)
    ))
  
  plot({'data':traces, 'layout': {'xaxis': {'title': 'Period'}, 'yaxis': {'title': 'Count'}}}, filename=output)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='plot gene expression')
  parser.add_argument('--gene', required=True, help='which gene to plot (or total)')
  parser.add_argument('--replicates', type=int, default=3, required=True, help='files per group')
  parser.add_argument('--files', nargs='+', help='files containing counts')
  parser.add_argument('--output', default='output.html', help='output file')
  args = parser.parse_args()
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)

  generate(args.gene, args.files, args.replicates, args.output)
