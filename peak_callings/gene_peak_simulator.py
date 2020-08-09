#!/usr/bin/env python

import pickle
import os
import pysam
import random
import pandas as pd
import logging
from multiprocessing import Pool
from scipy.stats import rv_discrete
from functools import partial
from collections import Counter
import numpy as np
from pybedtools import BedTool, set_tempdir
outdir = '/stor/scratch/Lambowitz/simulated_peaks_IDR'
if not os.path.isdir(outdir):
    os.mkdir(outdir)
set_tempdir(outdir)
np.random.seed(seed=123)
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('Peak simulator')


class GenePeakGenerator:
    def __init__(self):
        count_file = '/stor/work/Lambowitz/yaojun/Work/cfNA/tgirt_map/bed_files/merged_bed/MACS2/annotated/gene_count.bed'
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
        .assign(end = lambda d: d.start + sizes.size(number_of_peaks))
    outfile = outdir + '/simulated.%i.bed' %i
    simulated_peaks.to_csv(outdir + '/simulated.%i.bed' %i, 
                           sep='\t', 
                           header=False, 
                           index=False)
    rbp_count = RBP(simulated_peaks, IDR=IDR)

    if (i + 1) % (n_simulation // 10)  == 0:
        logger.info('Completed simulation %i with %i RBP' %(i + 1, rbp_count))

    return rbp_count


def main():
    n_simulation = 5000
    IDR = True
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
    main()

