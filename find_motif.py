#!/usr/bin/env python

import argparse
import collections
import logging
import sys
import textwrap

def extract(name, fasta, start, finish):
  sys.stdout.write('>{}\n'.format(name))
  start = max(0, start)
  finish = min(finish, len(fasta))
  sys.stdout.write('\n'.join(textwrap.wrap(fasta[start:finish], 70)))
  sys.stdout.write('\n')
  
def generate(genes, fasta, gff3, offset):
  logging.info('reading fasta...')
  contigs = {}
  contig = None
  for line in open(fasta, 'r'):
    if line.startswith('>'):
      if contig is not None:
        contigs[contig] = ''.join(current)
        logging.info('added contig %s', contig)
      contig = line.strip().split(' ')[0][1:]
      current = []
    else:
      current.append(line.strip())
  contigs[contig] = ''.join(current)

  logging.info('reading features...')
  gene_set = set(genes)
  logging.info(gene_set)
  seen = set()
  transcripts = {}
  # 1	araport11	five_prime_UTR	3631	3759	.	+	.	Parent=transcript:AT1G01010.1
  # 2 araport11 gene  19245591  19248915  . + . ID=gene:AT2G46830;Name=CCA1;biotype=protein_coding;description=Protein CCA1
  for line in open(gff3, 'r'):
    if line.startswith('#'):
      continue
    fields = line.strip().split('\t')
    chromosome = fields[0]
    feature = fields[2]
    if feature == 'gene':
      if 'Name=' not in fields[8]:
        continue
      name = fields[8].split('Name=')[1].split(';')[0] # gene
      if name in gene_set and name not in seen:
        transcript = fields[8].split('ID=gene:')[1].split(';')[0]
        logging.debug('found transcript %s for gene %s', transcript, name)
        transcripts[transcript] = name
    if feature == 'five_prime_UTR':
      transcript = fields[8].split('Parent=transcript:')[1].split('.')[0]
      if transcript in transcripts and transcripts[transcript] in gene_set and transcripts[transcript] not in seen:
        start = int(fields[3])
        finish = int(fields[4])
        strand = fields[6]
        gene = transcripts[transcript]
        logging.debug('extracting sequence for %s', gene)
        if strand == '+':
          extract(gene, contigs[chromosome], start-offset, start)
        else:
          extract(gene, contigs[chromosome], finish, finish+offset)

        seen.add(gene)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='plot gene expression')
  parser.add_argument('--genes', required=True, nargs='+', help='genes to extract')
  parser.add_argument('--fasta', required=True, help='fasta file to extract from')
  parser.add_argument('--gff3', required=True, help='gff file to find features in')
  parser.add_argument('--offset', default=1000, help='gff file to find features in')
  args = parser.parse_args()
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)

  generate(args.genes, args.fasta, args.gff3, args.offset)
