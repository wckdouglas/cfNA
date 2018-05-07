#!/bin/env python

# This is the program for converting bed file to
# bam file and calculate window prediction scores from the selected window size
# input: bed file for alignment, each record is a fragment
#        gtf file
#        fai file

from __future__ import print_function
import numpy as np
import argparse
import pysam
import os
import sys
from sys import stderr, argv
from functools import partial
import pyBigWig as pbw
from numba import jit
import gc

class bigWigFile:
    def __init__(self, filename):
        self.bw_file = pbw.open(filename, 'w')

    def add_chrom(self, chrom_dict):
        self.bw_file.addHeader(list(chrom_dict.items()))

    
    def add_wps(self, chromosome, chrom_array):
        chrom_length = len(chrom_array)
        self.bw_file.addEntries(chromosome, 
                        list(range(chrom_length)), 
                        values=np.array(chrom_array, dtype=np.float64),
                        validate = False, 
                        span=1)

    def close(self):
        self.bw_file.close()


### helper functions
def printMessage(message, sample):
    programname = os.path.basename(argv[0]).split('.')[0]
    print('[%s] %s: %s' %(programname,sample,message), file=sys.stderr)
    return 0


def make_regions(chromosome_length, how_many_bases_to_look_at, half_window):
    start = half_window
    end = start + how_many_bases_to_look_at
    upper_limit = chromosome_length - half_window
    while end < upper_limit:
        yield int(start), int(end)
        start = end
        end = end + how_many_bases_to_look_at
    yield int(start), int(upper_limit)


# get options
def getOpt():
    parser = argparse.ArgumentParser(description='Extract coverage of TSS')
    parser.add_argument('-i','--inFile',help='input bed file (merged paired-end fragments to single fragment bed file)',required=True)
    parser.add_argument('-o','--outprefix',help='output prefix', default='out')
    parser.add_argument('-g','--genome',help='genome file (genome.fa)',required=True)
    parser.add_argument('-w','--window', help='Window size for calculating WPS in memory (default: 10000)', default = 100000, type=int)
    args = parser.parse_args()
    inFile = args.inFile
    outprefix = args.outprefix
    genome = args.genome
    window = args.window
    return inFile, outprefix, genome, window 

@jit()
def push_WPS_to_Array(fields, halfWPSwindow, start, end, window, isize, wpsWindow):
    """
    for a given alignment, compute the regions that can be fully aligned and not
    e.g. [-1, -1 , -1, 1 , 1, 1, 1, 1, 1, -1, -1, -1] for a wps window -f 6 (halfwindow 3 )
    this will be added to the defined transcription start site wps score array after
    adjusting fot the position
    """
    window = end - start
    transcriptAlnWPS = np.zeros(window) # setting the tss window as zeros wps array
    alnWPS = np.zeros(isize + wpsWindow) #adding halfwindow to both side of the alignment
    alnWPS[wpsWindow:-wpsWindow] = 1 # only half window after aln start, the alignment is fully covering the wps window
                                     # and we added half window on the previous line
    alnWPS[alnWPS != 1 ] = -1 #making the alignment wps with ends containing windows
    alnShift = start - (int(fields[1]) - halfWPSwindow) # the distance between alignment start and right side of the wps window:
                                  #  + is left side of the window start, - is to the right
    if alnShift >= 0:
        wps = alnWPS[alnShift:]
        frag_end = len(wps) if len(wps) < window else window
        transcriptAlnWPS[:frag_end] += wps[:frag_end]
    else:
        baseShifted = abs(alnShift)
        frag_end = window if baseShifted + len(alnWPS) > window else baseShifted + len(alnWPS)
        alignedBases = window + alnShift
        wps = alnWPS[:int(alignedBases)]
        transcriptAlnWPS[int(baseShifted):int(frag_end)] += wps
    return transcriptAlnWPS

@jit()
def calculate_WPS(aln_file, chrom, window, wpsWindow, halfWPSwindow, upperBound, lowerBound, start, end):
    '''
    for each gene start site region:
        for i in each position:
            fetch all alignments in the region with the window size (e.g. 120)
            calculate number of alignments fully span and partially mapped
            get wps[i] = number of full map - number of end mapped
    '''
    transcriptWpsForward = np.zeros(end - start)
    transcriptWpsReverse = np.zeros(end - start)
    for aln in aln_file.fetch(reference = chrom, 
                            start = start - halfWPSwindow, 
                            end = end+halfWPSwindow):
        fields = aln.strip().split('\t')
        isize = int(fields[2]) - int(fields[1])
        strand = fields[5]
        if lowerBound <= isize <= upperBound:
            if strand == '-':
                transcriptWpsReverse += push_WPS_to_Array(fields, halfWPSwindow, start, end, window, isize, wpsWindow)
            else:
                transcriptWpsForward += push_WPS_to_Array(fields, halfWPSwindow, start, end, window, isize, wpsWindow)
    return transcriptWpsForward, transcriptWpsReverse


def extract_aln(bam, window, wpsWindow, halfWPSwindow, upperBound,
        lowerBound, chrom, chromSize, lenType, samplename):
    '''
    adding up wps track for all genes
    '''
    chromArrayForward = np.zeros(chromSize)
    chromArrayReverse = np.zeros(chromSize)
    cal_wps = partial(calculate_WPS, bam, chrom,window, wpsWindow, halfWPSwindow, upperBound, lowerBound)
    for start, end in make_regions(chromSize, window, halfWPSwindow):
        forwardWps, reverseWps = cal_wps(start, end)
        chromArrayForward[start:end] += forwardWps
        chromArrayReverse[start:end] += reverseWps

    printMessage('Finished calculating %s WPS for chromosome %s' %(lenType, chrom), samplename)
    return chromArrayForward, chromArrayReverse


def parse_faidx(genome):

    regular_chrom = list(range(1,23))
    regular_chrom = []
    regular_chrom.extend(['X','Y'])
    regular_chrom = ['chr' + str(chrom) for chrom in regular_chrom]

    fa = pysam.Fastafile(genome)
    chroms = {}
    for chrom, chrom_length in zip(fa.references, fa.lengths):
        if chrom in regular_chrom:
            chroms[chrom] = chrom_length
    return chroms


def make_output(outprefix, lenType, chrom_lengths):
    '''
    define output files
    '''
    file_template = '{outprefix}.{strand}.bigWig'

    out_bw = {}
    for strand in ['fwd','rvs']:
        bigwig_name = file_template.format(strand = strand, 
                                    outprefix = outprefix)
        bw = bigWigFile(bigwig_name)
        bw.add_chrom(chrom_lengths)
        out_bw[strand] = bw
    
    return out_bw


def runFile(bed, outprefix, genome, wpsWindow, window, upperBound, lowerBound, lenType, samplename):
    wpsWindow = wpsWindow + 1
    halfWPSwindow = np.divide(wpsWindow,2)
    
    # fetch info and define output
    chrom_lengths = parse_faidx(genome)
    out_bws = make_output(outprefix, lenType, chrom_lengths)

    with pysam.Tabixfile(bed) as tabix_bed:
        '''
        For each chromosome, extract WPS
        '''
        for chromosome in chrom_lengths.keys():
            chrom_size = chrom_lengths[chromosome]
            chromArrayForward, chromArrayReverse = extract_aln(tabix_bed, window, wpsWindow, halfWPSwindow, upperBound,
                    lowerBound, chromosome, chrom_size, lenType, samplename)

            out_bws['fwd'].add_wps(chromosome, chromArrayForward)
            out_bws['rvs'].add_wps(chromosome, chromArrayReverse)
            printMessage('Written %s to BigWig' %(chromosome), samplename)
            
            del chromArrayForward, chromArrayReverse
            gc.collect()

    #close all
    [bw.close() for bw in out_bws.values()]
    return 0

def main(inFile, outprefix, genome, window):
    '''
    main function for controling the work flow
    '''
    samplename = os.path.basename(inFile).split('.')[0]
    printMessage( 'Saving all result to: %s' %outprefix, samplename)

    lowerBound = 15
    upperBound = 100
    lenType = 'Short (%i-%ibp)' %(lowerBound, upperBound)
    wps_window = 5

    runFile(inFile, outprefix, genome, wps_window, window, upperBound, lowerBound, lenType, samplename)
    map(runFile, args)
    return 0

if __name__ == '__main__':
    inFile, outprefix, genome, window = getOpt()
    main(inFile, outprefix, genome, window)
