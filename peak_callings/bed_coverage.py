#!/usr/bin/env python

import pyBigWig as pbw
import pysam
import numpy as np
from operator import itemgetter
import pickle
import sys
from tqdm import tqdm

class bed_coverage:
    def __init__(self, inbed, fasta_file, out_bw, strand):
        self.fa  = pysam.Fastafile(fasta_file)
        self.headers = list(zip(self.fa.references, self.fa.lengths))
        self.inbed = pysam.Tabixfile(inbed)
        self.bw = pbw.open(out_bw, 'w')
        self.bw.addHeader(self.headers)
        self.strand = strand
        assert(self.strand in ['+','-'])
        print('Using %s ' %strand)

        self.bed_to_coverage()
        self.bw.close()


    def bed_to_coverage(self):
        self.cov_dict = {}
        for chrom, chrom_length in tqdm(self.headers):
            chrom_cov = np.zeros(chrom_length)
            print('Initialized chromosome %s' %chrom)
        
            # generate coverage
            for line in self.inbed.fetch(chrom, 0, chrom_length):
                fields = line.strip().split('\t')
                chrom, start, end, fraction, strand = itemgetter(0,1,2, 4, 5)(fields)
                if strand == self.strand:
                    chrom_cov[int(start):int(end)] += float(fraction)

            # put in bigwig
            self.bw.addEntries(chrom, 
                          list(range(chrom_cov.shape[0])), 
                          values= chrom_cov.tolist(), 
                          span=1)
            print('Finished chrom: %s' %chrom)





if len(sys.argv) != 5:
    sys.exit('[usage] python %s <tabix file> <fasta file> <out_bigwig> <strand>' %sys.argv[0])
inbed = sys.argv[1]
fasta_file = sys.argv[2]
out_bw = sys.argv[3]
strand = sys.argv[4]

bc = bed_coverage(inbed, fasta_file, out_bw, strand)
bc.make_bigwig()
bc.bed_to_coverage()
print('Written:', out_bw)

