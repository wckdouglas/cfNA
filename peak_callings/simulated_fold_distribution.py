#!/usr/bin/env python

import pickle
import glob
import os
import pysam
import logging
from multiprocessing import Pool
import numpy as np
import RNA
from sequencing_tools.fastq_tools import reverse_complement
outdir = '/stor/scratch/Lambowitz/simulated_peaks_long_gene'
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('Peak simulator')

class Fold():
    def __init__(self):
        '''
        take in peak coordinates, return minimum MFE between two strands
        '''
        fa = '/stor/work/Lambowitz/ref/hg19_ref/genome/hg19_genome.fa'
        self.fa = pysam.FastaFile(fa)

    def fold(self,chrom, start, end):
        start, end = int(start), int(end)
        seq_fwd = self.fa.fetch(chrom, start, end)
        seq_rvs = reverse_complement(seq_fwd)
        _, fwd_MFE = RNA.fold(seq_fwd)
        _, rvs_MFE = RNA.fold(seq_rvs)
        return min(fwd_MFE, rvs_MFE)

def processBED(bed):
    fold = Fold()
    MFE = []
    i = bed.split('.')[-2]
    if (int(i)+1) % 20 == 0:
        logger.info('Reading %s' %bed)
    with open(bed) as inbed:
        for line in inbed:
            fields = line.split('\t')
            chrom, start, end = fields[:3]
            MFE.append(fold.fold(chrom, start, end))
    return MFE


def main():
    beds = glob.glob(outdir + '/*bed')
    beds.sort()
    p = Pool(24)
    MFEs = p.map(processBED, beds[:20])
    p.close()
    p.join()

    outfile = outdir + '/peak_folding.pickle'
    with open(outfile, 'wb') as f:
        pickle.dump(np.array(MFEs), f)
    logger.info('Written %s' %outfile)


if __name__ == '__main__':
    main()

