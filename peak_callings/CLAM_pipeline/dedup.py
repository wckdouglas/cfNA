#!/usr/bin/env python

from __future__ import print_function, division
from itertools import groupby
from operator import itemgetter
from sequencing_tools.io_tools import xopen
from collections import defaultdict
import six
import sys
import numpy as np
import argparse


def getopt():
    parser = argparse.ArgumentParser(description = 'adding tag to bam')
    parser.add_argument('-i', '--inbed', required=True, help = 'input bed file')
    parser.add_argument('-o','--outbed', default = '-', help = 'output bed file')
    parser.add_argument('-p','--prefix', default = None, help = 'prefix to read name')
    args = parser.parse_args()
    return args


def group_fragments(line):
    fields = line.strip().split('\t')
    chrom, start, end, strand = itemgetter(0,1,2,5)(fields)
    return chrom, start, end, strand


def main():
    args = getopt()

    in_bed = xopen(args.inbed) if args.inbed != "-" else sys.stdin
    out_bed = xopen(args.outbed, 'w') if args.outbed != "-" else sys.stdout
    if args.prefix:
        read_template = '{chrom}\t{start}\t{end}\t%s_{umi}\t{pseudo_count}\t{strand}' %args.prefix
    else:
        read_template = '{chrom}\t{start}\t{end}\t{umi}\t{pseudo_count}\t{strand}' 
    for (chrom, start, end, strand), lines in groupby(in_bed, group_fragments):
        umis = defaultdict(list)
        for line in lines:
            fields = line.strip().split('\t')
            umi = fields[3].split('_')[0]
            pseudo_count = float(fields[-1])
            umis[umi].append(pseudo_count)
    

        for umi, pseudo_counts in six.iteritems(umis):
            print(read_template \
                .format(chrom = chrom,
                        start = start,
                        end = end,
                        umi = umi,
                        pseudo_count = np.mean(pseudo_counts),
                        strand = strand), file = out_bed)
            


if __name__ == '__main__':
    main()
