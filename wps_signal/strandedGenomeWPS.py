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
import subprocess
import os
import sys
from sys import stderr, argv
from functools import partial
from multiprocessing import Pool
from itertools import izip
import pyBigWig as pbw
from numba import jit

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
        yield (start, end)
        start = end
        end = end + how_many_bases_to_look_at
    yield (start, upper_limit)


# get options
def getOpt():
    parser = argparse.ArgumentParser(description='Extract coverage of TSS')
    parser.add_argument('-i','--inFile',help='input bed file (merged paired-end fragments to single fragment bed file)',required=True)
    parser.add_argument('-o','--outprefix',help='output prefix', default='out')
    parser.add_argument('-g','--genome',help='genome file (.fa.fai)',required=True)
    parser.add_argument('-c','--chromosome',help='chromosome name',required=True)
    parser.add_argument('-w','--window', help='Window size for calculating WPS in memory (default: 10000)', default = 10000, type=int)
    args = parser.parse_args()
    inFile = args.inFile
    outprefix = args.outprefix
    genome = args.genome
    window = args.window
    chromosome = args.chromosome
    return inFile, outprefix, genome, window, chromosome

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
    alnShift = start - (long(fields[1]) - halfWPSwindow) # the distance between alignment start and right side of the wps window:
                                  #  + is left side of the window start, - is to the right
    if alnShift >= 0:
        wps = alnWPS[alnShift:]
        end = len(wps) if len(wps) < window else window
        transcriptAlnWPS[:end] += wps[:end]
    else:
        baseShifted = abs(alnShift)
        end = window if baseShifted + len(alnWPS) > window else baseShifted + len(alnWPS)
        alignedBases = window + alnShift
        wps = alnWPS[:alignedBases]
        transcriptAlnWPS[baseShifted:end] += wps
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
        isize = long(fields[2]) - long(fields[1])
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

def writeWig(chrom_array, outputWig, chromosome, samplename):
    outputWig = outputWig.replace('|','_')
    outWig =  pbw.open(outputWig,'w')
    chrom_length = len(chrom_array)
    outWig.addHeader([(chromosome,chrom_length)])
    outWig.addEntries(chromosome, range(chrom_length), values=np.array(chrom_array, dtype=np.float64), span=1)
    outWig.close()
    printMessage('Witten %s' %outputWig, samplename)
    return 0


def parse_faidx(genome):
    chroms = {}
    for line in open(genome):
        fields = line.split('\t')
        chrom, chrom_size = fields[0], fields[1]
        chroms[chrom] = long(chrom_size)
    return chroms


def runFile(arg):
    bed, outprefix, genome, wpsWindow, window, upperBound, lowerBound, lenType, samplename, chromosome = arg
    wpsWindow = wpsWindow + 1
    halfWPSwindow = np.divide(wpsWindow,2)
    output_forward_wig= outprefix + '.' + chromosome+'.'+lenType.split(' ')[0] +'.forward.bigWig'
    output_reverse_wig = output_forward_wig.replace('forward','reverse')
    chromLength = parse_faidx(genome)
    with pysam.Tabixfile(bed) as tabix_bed:
        if chromosome not in chromLength.keys():
            sys.exit('Wrong chromosome name: %s!' %chromosome)
        chromSize = chromLength[chromosome]
        chromArrayForward, chromArrayReverse = extract_aln(tabix_bed, window, wpsWindow, halfWPSwindow, upperBound,
                lowerBound, chromosome, chromSize, lenType, samplename)

    for output_wig, chromArray in izip([output_forward_wig, output_reverse_wig],
                                        [chromArrayForward, chromArrayReverse]):
        writeWig(chromArray, output_wig, chromosome, samplename)
    return 0

def main(inFile, outprefix, genome, window, chromosome):
    '''
    main function for controling the work flow
    '''
    samplename = os.path.basename(inFile).split('.')[0]
    printMessage( 'Saving all result to: %s' %outprefix, samplename)

    iterater =  zip([100, 180],[10, 120],['Short (30-80bp)','Long (120-180bp)'],[5,120])
    args = [(inFile, outprefix, genome, wps_window, window, upperBound, lowerBound, lenType, samplename, chromosome) \
            for upperBound, lowerBound, lenType, wps_window in iterater]
    map(runFile, args)
    return 0

if __name__ == '__main__':
    inFile, outprefix, genome, window, chromosome = getOpt()
    main(inFile, outprefix, genome, window, chromosome)
