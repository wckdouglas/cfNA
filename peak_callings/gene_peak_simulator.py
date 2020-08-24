#!/usr/bin/env python

import pickle
import glob
import os
import pysam
import random
import pandas as pd
import logging
from multiprocessing import Pool
from scipy.stats import rv_discrete
from functools import partial
from collections import Counter
from operator import itemgetter
import numpy as np
from pybedtools import BedTool, set_tempdir
outdir = '/stor/scratch/Lambowitz/simulated_peaks_all_gene'
if not os.path.isdir(outdir):
    os.mkdir(outdir)
set_tempdir(outdir)
np.random.seed(seed=123)
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('Peak simulator')


class GenePeakGenerator:
    def __init__(self):
        count_file = '/stor/work/Lambowitz/yaojun/Work/cfNA/tgirt_map/bed_files/merged_bed/MACS2/annotated/gene_count.bed'
        count_file = '/stor/work/Lambowitz/yaojun/Work/cfNA/tgirt_map/Counts/all_counts/gene_count.bed'
        self.genes = []
        self.counts = []
        with open(count_file) as counts:
            for line in counts:
                if not line.startswith('.'):
                    fields = line.strip().split('\t')
                    key = ','.join([fields[0],fields[1],fields[2]])
                    self.genes.append(key)
                    self.counts.append(int(fields[-1]))


    def random_genes(self, number_of_peaks = 100):
        random_genes = Counter(random.choices(self.genes, 
                                       weights = self.counts,
                                       k=number_of_peaks))
        
        peaks = []
        for gene, count in random_genes.items():
            chrom, start, end = gene.split(',')
            for pos in np.random.randint(low = int(start), high = int(end), size = count):
                peaks.append((chrom, pos))
        return peaks


class RBPtabix():
    def __init__(self, IDR=False, overlap_cutoff = 11):
        if not IDR:
            bed = '/stor/work/Lambowitz/ref/hg19_ref/genes/filtered_RBP.bed.gz'
        else:
            bed = '/stor/work/Lambowitz/ref/hg19_ref/genes/RBP_IDR.reformated.bed.gz'

        self.bed = pysam.Tabixfile(bed)
        self.overlap_cutoff = overlap_cutoff
    
    def RBP_count(self, chrom, start, end):
        overlapping_base = 0
        for line in self.bed.fetch(chrom, start, end):
            field = line.split('\t')
            rbp_start, rbp_end = itemgetter(1,2)(field)
            overlap = self.__overlap__(start, end, int(rbp_start),int( rbp_end))
            if overlap >= self.overlap_cutoff:
                overlapping_base = max(overlapping_base, overlap)
        return overlapping_base

    def __overlap__(self, s1, e1, s2, e2):
        return min(e1, e2) - max(s1, s2)


def RBP(simulated_peaks, IDR=False):
    if not IDR:
        bed = '/stor/work/Lambowitz/ref/hg19_ref/genes/filtered_RBP.bed.gz'
    else:
        bed = '/stor/work/Lambowitz/ref/hg19_ref/genes/RBP_IDR.reformated.bed.gz'
    rbp_bed = BedTool(bed)
    input_peaks = simulated_peaks.shape[0]
    non_rbp_peak = simulated_peaks\
        .reset_index(drop=True)\
        .pipe(BedTool().from_dataframe)\
        .intersect(b=rbp_bed, v=True) \
        .to_dataframe()  
    
    return  input_peaks - non_rbp_peak.shape[0]

class SizeGenerator:
    def __init__(self, peak_sizes = [100,200,300]):
        '''
        A list of peak sizes

        return a distribution class
        '''

        # make insert distribution
        self.peak_sizes, self.peak_count = np.unique(peak_sizes, return_counts=True)
        px = self.peak_count/self.peak_count.sum()
        self.insert_dist = rv_discrete(name='custm', values=(self.peak_sizes, px))

    def size(self, number_of_peaks=100):
        return self.insert_dist.rvs(size = number_of_peaks)


def simulator(n_simulation, peak_sizes, number_of_peaks, IDR, i):
    '''
    running each chromosome
    '''
    gps = GenePeakGenerator()
    peaks = gps.random_genes(number_of_peaks)
    sizes = SizeGenerator(peak_sizes = peak_sizes)

    simulated_peaks = pd.DataFrame(peaks, columns = ['chrom','start']) \
        .assign(end = lambda d: d.start + sizes.size(number_of_peaks)) \
        .pipe(lambda d: countingRBP(11, d))
    outfile = outdir + '/simulated.%i.bed' %i
    simulated_peaks.to_csv(outdir + '/simulated.%i.bed' %i, 
                           sep='\t', 
                           header=False, 
                           index=False)
    rbp_count = simulated_peaks.query('rbp > 0').shape[0]

    if (i + 1) % (n_simulation // 10)  == 0:
        logger.info('Completed simulation %i with %i RBP' %(i + 1, rbp_count))

    return rbp_count

def countingRBP(overlap_cutoff, bed):
    rbptab = RBPtabix(overlap_cutoff = overlap_cutoff)
    return bed \
        .assign(rbp = lambda d: list(map(rbptab.RBP_count, d.chrom, d.start, d.end) )) 

def rbp_count(overlap_cutoff, bed):
    rbptab = RBPtabix(overlap_cutoff = overlap_cutoff)
    return pd.read_csv(bed, names = ['chrom','start','end'], sep='\t', usecols=[0,1,2]) \
        .assign(rbp = lambda d: list(map(rbptab.RBP_count, d.chrom, d.start, d.end) )) \
        .rbp.values


def main1():
    beds = glob.glob(outdir + '/*bed')
    overlap_cutoff = 0
    func = partial(rbp_count, overlap_cutoff)

    p = Pool(24)
    rbp_counts = p.map(func, beds)
    p.close()
    p.join()

    out_file = outdir + '/rbp_peaks_cutoff.pickle'
    with open(out_file,'wb') as f: 
        pickle.dump(np.array(rbp_counts), f)
    logger.info('Written %s' %out_file)



def main():
    n_simulation = 5000
    IDR = False
    merged_peak = pd.read_excel('Sup_file_061620.xlsx', sheet_name = 'MACS2 peaks (Genomic)')  \
            .pipe(lambda d: d[d['High confidence'] == "Y"])
    RBP_count = merged_peak.pipe(lambda d: d[d['Sense strand RBP']!= '.']).shape[0]
    number_of_peaks = merged_peak.shape[0]
    logger.info('%i peaks with %i RBP binding sites ' %(number_of_peaks, RBP_count))
    peak_sizes = merged_peak['Peak size'].tolist()
    simulator_func = partial(simulator, n_simulation, peak_sizes, number_of_peaks, IDR)

    p = Pool(24)
    result = p.map(simulator_func, range(n_simulation))
    p.close()
    p.join()

    RBP_counts = np.array(result)
    number_of_times_more_RBP = len(RBP_counts[RBP_counts >= RBP_count])
    p_value = number_of_times_more_RBP / n_simulation
    logger.info('RBP enrichment p-value: %.4f' %p_value)

    out_file = outdir + '/rbp_peaks.pickle'
    with open(out_file,'wb') as f: 
        pickle.dump(RBP_counts, f)


if __name__ == '__main__':
    main1()

