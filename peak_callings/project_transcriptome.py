#!/usr/bin/env python


import argparse
import sys
from operator import itemgetter
from sequencing_tools.gene_tools import Bed12Record

def getopt():
    parser = argparse.ArgumentParser(description='To convert bedpe file to bed file')
    parser.add_argument('-i','--input',default='-',help='Input narrowpeak bed')
    parser.add_argument('-o','--output',default='-',help='Output narrowpeak genomic coordinate bed file (default: -)')
    parser.add_argument('-b','--bed12', required=True, help ='bed 12 file for projecting peaks')
    args = parser.parse_args()
    return args


def index_transcripts(bed12):
    transcripts = {}
    with open(bed12, 'r') as bed:
        for t in bed:
            transcripts[t.split('\t')[3]] = Bed12Record(t)
    return transcripts



def main():
    args = getopt()
    fileinput = sys.stdin if args.input == '-' else open(args.input)
    fileoutput=  sys.stdout if args.output == '-' else open(args.output, 'w')
    transcripts = index_transcripts(args.bed12)
    for line in fileinput:
        fields = line.strip().split('\t')
        tid, peak_start, peak_end = itemgetter(0,1,2)(fields)
        fields[5] = transcripts[tid].strand
        peak_start = transcripts[tid].genomic_position(int(peak_start)+1)
        peak_end = transcripts[tid].genomic_position(int(peak_end))
        chrom = transcripts[tid].chrom
        if peak_start > peak_end:
            peak_end, peak_start = peak_start, peak_end
        outline = '{chrom}\t{start}\t{end}\t{fields}'\
            .format(chrom = chrom, start = peak_start,
                    end = peak_end, fields = '\t'.join(fields[3:]))
        print(outline, file = fileoutput)




if __name__ == '__main__':
    main()



