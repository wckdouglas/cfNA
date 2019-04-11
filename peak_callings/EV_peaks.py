#!/usr/bin/env python

import pandas as pd
import pysam
import sys

if len(sys.argv) != 3:
    sys.exit('python %s <peak_coor> <bed.gz file> ' %sys.argv[0])

class peak_count:
    '''
    parse peak coordinates (peak file),
    count tabix reads for each peak
    '''
    def __init__(self, peaks, tabix):
        self.tabix = pysam.Tabixfile(tabix)
        self.peak_df = pd.read_table(peaks) 

    def count_reads(self, chrom, start, end, strand):
        reads = self.tabix.fetch(chrom, start, end)
        read_count = sum(1 for r in reads if r.strip().split('\t')[5] == strand)
        return read_count

    def peak_counting(self):
        return self.peak_df \
            .assign(EV_count = lambda d: list(map(self.count_reads, 
                                                d.chrom, 
                                                d.start, 
                                                d.end, 
                                                d.strand))) 


pc = peak_count(sys.argv[1], sys.argv[2])
pc.peak_counting()\
    .to_csv(sys.stdout, sep='\t', index=False)

        
