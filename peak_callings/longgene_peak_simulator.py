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
from gene_peak_simulator import GenePeakGenerator, SizeGenerator
outdir = '/stor/scratch/Lambowitz/simulated_peaks_long_gene'
if not os.path.isdir(outdir):
    os.mkdir(outdir)
set_tempdir(outdir)
np.random.seed(seed=123)
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('Peak simulator')


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
    
    if (i+1) % 10 == 0:
        logger.info('Simulated: %i set peaks' %(i+1))

def main():
    n_simulation = 100
    IDR = False
    merged_peak = pd.read_excel('Sup_file_061620.xlsx', sheet_name = 'MACS2 peaks (Genomic)')  \
            .pipe(lambda d: d[d['High confidence'] == "Y"]) \
            .pipe(lambda d: d[d['Peak type'].isin(['RBP', 'Long RNA'])])
    number_of_peaks = merged_peak.shape[0]
    logger.info('%i longRNA peaks with ' %(number_of_peaks))
    peak_sizes = merged_peak['Peak size'].tolist()
    simulator_func = partial(simulator, n_simulation, peak_sizes, number_of_peaks, IDR)

    p = Pool(24)
    result = p.map(simulator_func, range(n_simulation))
    p.close()
    p.join()
    logger.info('Simulated %i set long RNA peaks' %n_simulation)


if __name__ == '__main__':
    main()

