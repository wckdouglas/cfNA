#!/usr/bin/env python

from __future__ import print_function
from Bio import Entrez, SeqIO
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
import sys
import pysam

if len(sys.argv) != 2:
    sys.exit('[usage] python %s <genome_fa>' %sys.argv[0] )

rRNA_gene = ['gi|23898|emb|X12811.1| Human 5S DNA',
    'gi|555853|gb|U13369.1|HSU13369 Human ribosomal DNA complete repeating unit']
Entrez.email = 'wckdouglas@gmail.com'

for rRNA in rRNA_gene:
    id = rRNA.split('|')[1]
    handle = Entrez.efetch(db="nucleotide", id=id, 
                           rettype="fasta", retmode="text")
    record = handle.read()
    print('>' + rRNA +'\n'+ ''.join(record.split('\n')[1:]))



fa = pysam.Fastafile(sys.argv[1])

for chrom, start, end, name, strand in [('chrM',648, 1601, 'MT-RNR1','+'), 
                                        ('chrM', 1671, 3229, 'MT-RNR2','+')]:
    seq = fa.fetch(chrom, start, end)
    print('>{name}\n{seq}'.format(name = name,seq = seq))
