#!/usr/bin/env python

from __future__ import print_function
from builtins import range, zip
import pandas as pd
import numpy as np
from scipy.signal import savgol_filter
import os
import sys
import pyBigWig as pbw
import glob
from functools import partial
import argparse
import pyximport
pyximport.install(setup_args={'include_dirs': np.get_include()})
from call_peak_tools import *
from operator import itemgetter
from multiprocessing import Pool
import re
#from memory_profiler import profile
import gc

def get_opt():
    parser = argparse.ArgumentParser(description='Given a WPS file in bigwig format '+\
                                     'output peak coordinates in  bed file')
    parser.add_argument('-i', '--in_bigwig', help = 'Input bigWig', required=True)
    parser.add_argument('-c', '--control_bigwig', help = 'Control bigWig')
    parser.add_argument('-o', '--out_bed', help = 'Output prefix', required=True)
    parser.add_argument('-s','--strand', help='Which strand?, all peak would have the same strand',
                        choices = ['forward','reverse'], required=True)
    parser.add_argument('--two_pass', help='2 pass peak calling, first pass will be very stringent', 
                        action='store_true')
    parser.add_argument('--threads', type=int, default = 24, help='Number of threads to use (default: 24)')
    args = parser.parse_args()
    return args


def process_bigwig(out_bed, inputWig, controlWig, strand, two_pass, chromosome):
    # get bigwig information
    filename = os.path.basename(inputWig)

    #print message
    print('Running %s for chrom: %s' %(filename, chromosome), file=sys.stderr)

    # read in data
    bw = pbw.open(inputWig)
    length = bw.chroms()[chromosome]
    wps = bw.values(chromosome, 0, length, numpy=True)
    bw.close()
    gc.collect()
    print('[%s] Collected wps' %(filename), file=sys.stderr)
    print('[%s] Start calling peaks: Chrom: %s' %(filename, chromosome), 
            file=sys.stderr)

    peak_count = 0
    temp_bed = out_bed + '.%s_temp' %(chromosome)
    with open(temp_bed, 'w') as peak_bed:
        second_pass = not two_pass # if two pass is True, this is the first pass, only output high quality peak
        peak_count = write_short_peaks(wps, controlWig, peak_bed, chromosome, strand, second_pass = second_pass, peak_count=peak_count)
        ## second_pass = TRUE  out put everything
        ## else, only output high quanlity peak
    print('Written %i peaks to %s' %(peak_count, temp_bed), file=sys.stderr)

    if two_pass:
        ## two pass algorithm, filter out first pass peaks, and do peak calling again
        wps = filter_called_peaks(wps, temp_bed)
        print('[%s] removed 1st pass peak regions' %(filename), file=sys.stderr)
        with open(temp_bed, 'a') as peak_bed:
            peak_count += write_short_peaks(wps, controlWig, peak_bed, chromosome, strand, second_pass = True, peak_count=peak_count)
        print('Written %i peaks to %s after 2nd pass' %(peak_count, temp_bed), file=sys.stderr)
    return temp_bed
    
def get_chroms(bw_file):
    '''
    Get chromosomes and lengths
    '''
    bw = pbw.open(bw_file)
    chromosomes_dict = bw.chroms()
    bw.close()
    regular_chromosome = re.compile('chr[0-9]+$|chr[MXY]$')
    chromosome_dict = {key:value for key,value in chromosomes_dict.items() if regular_chromosome.search(key)}
    return chromosomes_dict

def main():
    args = get_opt()
    bigwig_file = args.in_bigwig
    bigwig_control = args.control_bigwig
    strand = args.strand
    strand = '+' if strand == 'forward' else '-'
    two_pass = args.two_pass

    chromosome_dict = get_chroms(bigwig_file)
    peak_func = partial(process_bigwig, args.out_bed, bigwig_file, bigwig_control, strand, two_pass) 

    p = Pool(args.threads)
    temp_files = p.imap(peak_func, chromosome_dict.keys())
    p.close()
    p.join()

    peak_count = 0
    with open(args.out_bed, 'w') as bed:
        for temp in temp_files:
            for i, line in enumerate(open(temp)):
                bed.write(line)
            os.remove(temp)
            peak_count += i
    print('Written %i peaks to %s' %(peak_count,args.out_bed), file = sys.stderr)
    return 0

if __name__ == '__main__':
    main()
