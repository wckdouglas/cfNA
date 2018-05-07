
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
pyximport.install()
from call_peak_tools import *
from operator import itemgetter
#from memory_profiler import profile
import gc

def get_opt():
    chromosomes = list(range(1,23))
    chromosomes.extend(list('XY'))
    chromosomes = map(lambda x: 'chr' + str(x), chromosomes)
    parser = argparse.ArgumentParser(description='Given a WPS file in bigwig format '+\
                                     'output peak coordinates in  bed file')
    parser.add_argument('-i', '--in_bigwig', help = 'Input bigWig', required=True)
    parser.add_argument('-c', '--chrom',  help = 'Chromosome', required=True)#, choices = chromosomes)
    parser.add_argument('-o', '--out_bed', help = 'Output prefix', required=True)
    parser.add_argument('-l','--length_type', help = 'short or long WPS?', default = 'Long',
                        choices = ['Short','Long'])
    parser.add_argument('-s','--strand', help='Which strand?, all peak would have the same strand',
                        choices = ['forward','reverse'], required=True)
    parser.add_argument('--two_pass', help='2 pass peak calling, first pass will be very stringent', 
                        action='store_true')
    args = parser.parse_args()
    return args


def calling_long_peaks(chromosome, wps, peak_start, peak_end, peak_count, outFile, 
                    strand, peak_size_filter = True, second_pass=False):
    '''
        using peak start and end from wps array,
        find maximum sub array
        and export start and end from maximum subarray
        peak score determine from maximum wps score.
    '''
    adjusted_sub_wps = adjust_peaks(wps, peak_start, peak_end)
    above_median_starts, above_median_ends = find_peak_region(adjusted_sub_wps)

    if len(above_median_starts)>len(above_median_ends):
        above_median_ends = np.append(above_median_ends,len(adjusted_sub_wps))
    if not peak_size_filter:
        above_median_starts, above_median_ends =  pick_peak(above_median_starts, above_median_ends, adjusted_sub_wps)

    for above_median_start, above_median_end in zip(above_median_starts, above_median_ends):
        sub_peak_wps = adjusted_sub_wps[above_median_start:above_median_end]
        nucleosome_start , nucleosome_end = maximum_sum_subarray(sub_peak_wps)


        #adjust coordinate
        nucleosome_start, nucleosome_end = peak_start + above_median_start + np.array([nucleosome_start, nucleosome_end])
        nucleosome_center = int((nucleosome_start + nucleosome_end) /2)
        peak_center = (nucleosome_start + nucleosome_end)/2
        nucleosome_size = abs(nucleosome_end - nucleosome_start)
        if (peak_size_filter and 50 < nucleosome_size  < 150 ) or (not peak_size_filter and nucleosome_size > 5):
            peak_score = wps[nucleosome_start:nucleosome_end].max()
            peak_count += 1
            peak_name = '%s_peak%i' %(chromosome, peak_count)
            line = '\t'.join(map(str,[chromosome, nucleosome_start, nucleosome_end,
                                      peak_name, peak_score, strand, peak_center]))
            print(line, file = outFile)
    return peak_count


def write_long_peaks(wps, out_bed, chromosome, strand, second_pass = False, peak_count = 0 ):
    peaks = peak_iterator(wps)
    peak_regions = merge_and_find_peaks(peaks, tolerance_unprotected = 10)
    for peak_start, peak_end  in peak_regions:
        peak_size = np.abs(peak_end - peak_start)
        if 40 <= peak_size <= 150:
            peak_count = calling_long_peaks(chromosome, wps, peak_start, peak_end, peak_count,
                                            out_bed, strand, peak_size_filter = False, second_pass=False)
        elif 150 < peak_size <= 450:
            peak_count = calling_long_peaks(chromosome, wps, peak_start, peak_end, peak_count,
                                            out_bed, strand, peak_size_filter = True, second_pass=False)
    return peak_count


def process_bigwig(out_bed, inputWig, length_type, strand, two_pass, chromosome):
    # get bigwig information
    filename = os.path.basename(inputWig)

    #print message
    print('Running %s as %s for chrom: %s' %(filename, length_type, chromosome), file=sys.stderr)

    # read in data
    bw = pbw.open(inputWig)
    chrom, length = bw.chroms().items()[0]
    assert chrom == chromosome, 'Wrong chromosomes'
    wps = bw.values(chrom,0,length, numpy=True)
    bw.close()
    if length_type == 'Long':
        wps = adjust_median(wps, window=1000)
        wps = savgol_filter(wps, window_length = 21, polyorder=2)
        peak_caller = partial(write_long_peaks)
    elif length_type == 'Short':
        #wps = adjust_median(wps, window=1000)
        peak_caller = partial(write_short_peaks)
    gc.collect()
    print('[%s] Adjusted wps' %(filename), file=sys.stderr)
    print('[%s] Start calling peaks: Chrom: %s, %s' %(filename, chromosome, length_type), file=sys.stderr)

    peak_count = 0
    with open(out_bed, 'w') as peak_bed:
        second_pass = not two_pass # if two pass is True, this is the first pass, only output high quality peak
        peak_count = peak_caller(wps, peak_bed, chrom, strand, second_pass = second_pass, peak_count=peak_count)
        ## second_pass = TRUE  out put everything
        ## else, only output high quanlity peak
    print('Written %i peaks to %s' %(peak_count, out_bed), file=sys.stderr)

    if two_pass:
        ## two pass algorithm, filter out first pass peaks, and do peak calling again
        wps = filter_called_peaks(wps, out_bed)
        print('[%s] removed 1st pass peak regions' %(filename), file=sys.stderr)
        with open(out_bed, 'a') as peak_bed:
            peak_count += peak_caller(wps, peak_bed, chrom, strand, second_pass = True, peak_count=peak_count)
        print('Written %i peaks to %s after 2nd pass' %(peak_count, out_bed), file=sys.stderr)
    return 0
    

def main():
    args = get_opt()
    bigwig_file = args.in_bigwig
    length_type = args.length_type
    strand = args.strand
    chromosome = args.chrom
    strand = '+' if strand == 'forward' else '-'
    two_pass = args.two_pass
    peak_func = partial(process_bigwig, args.out_bed, bigwig_file, length_type, strand, two_pass) 
    peak_files = peak_func(chromosome)
    return 0

def alt_main():
    bigwig_file = '/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/bed_files/genome_WPS/200_NT_3_S7_R1_001_no_sncRNA_repeats.chrX.Short.forward.bigWig'
    length_type = 'Short'
    strand = 'forward'
    chromosome = 'chrX'
    strand = '+' if strand == 'forward' else '-'
    two_pass = False
    out_bed = '/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/genome_peaks/200_PI_1_S11_R1_001.chr19.Short.forward.bed'
    peak_func = partial(process_bigwig, out_bed, bigwig_file, length_type, strand, two_pass) 
    peak_files = peak_func(chromosome)

if __name__ == '__main__':
#    alt_main()
    main()
