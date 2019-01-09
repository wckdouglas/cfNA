#!/usr/bin/env python

import pandas as pd
from skbio.io import read
from itertools import groupby
from collections import Counter
import re
from sequencing_tools.fastq_tools import reverse_complement


tRNA_PATH = '/stor/work/Lambowitz/ref/hg19_ref/genes'
mature_fa = tRNA_PATH + '/tRNA/hg19-mature-tRNAs.fa'
tRNA_detail = tRNA_PATH + '/tRNA/hg19-tRNAs-confidence-set.ss'
tRNA_fa = tRNA_PATH + '/tRNA.fa'


def tRNA_seq_dict():
    tRNA_dict = {}
    for r in read(mature_fa, 'fasta'):
        description = r.metadata['description'].split(' ')[2].strip(')')
        tRNA_dict[str(r).replace('U','T')+'CCAA'] = description
    return tRNA_dict
    
def anticodon_pos():
    anticodon_dict = {}
    with open(tRNA_detail, 'r') as infile:
        try:
            while True:
                line = next(infile)
                if line.startswith('chr'):
                    tRNA_id = line.split(' ')[0]
                    line = next(infile)
                    pos = re.search('([0-9][0-9]-[0-9][0-9])',line).group(0)
                    anticodon = line.split('Anticodon: ')[1].split(' ')[0]
                    aa = line.split(' ')[1].split('\t')[0]
                    anticodon_dict[tRNA_id] = pos + ',' + anticodon + ',' + aa
        except StopIteration:
            pass
    return anticodon_dict


tRNA_dict = tRNA_seq_dict()
anticodon_dict = anticodon_pos()
rows = []
for record in read(tRNA_fa,'fasta'):
    if record.metadata['id'].startswith('TR'):
        tRNA_id = tRNA_dict[str(record)]
        anticodon = record.metadata['id'].split('-')[1]
        pos, annotated_anticodon, aa = anticodon_dict[tRNA_id].split(',')
        tRNA_length = len(record)
        rows.append((record.metadata['id'], 0, tRNA_length, anticodon, pos, aa, str(record) ))

tablename = tRNA_PATH + '/tRNA/anticodon_annotations.tsv'
df = pd.DataFrame(rows, columns = ['tRNA','start','end','anticodon','anticodon_pos', 'aa', 'seq']) \
    .assign(anticodon_start = lambda d: d.anticodon_pos.str.extract('^([0-9]+)-',expand=False).astype(int))\
    .assign(anticodon_end = lambda d: d.anticodon_pos.str.extract('-([0-9]+)$',expand=False).astype(int)) \
    .assign(predicted_anticodon = lambda d: list(map(lambda x,y,z: x[(y-1):z], d.seq, d.anticodon_start, d.anticodon_end)))\
    .to_csv(tablename, sep='\t', index=False)
print('Written %s' %tablename)

        
