#!/usr/bin/env python

import pysam
from bwapy import BwaAligner
import argparse
import logging
import sys

def getopt():
    parser = argparse.ArgumentParser(description = 'Filter read pairs with either read mapping to chrM')
    parser.add_argument('-i', '--inbam', 
                        required=True, default='-', 
                        help = 'input bam file')
    parser.add_argument('-o','--outbam', 
                        default = '-', 
                        help = 'output bam file (defulat: - )')
    parser.add_argument('--filtered_bam', 
                        default = '', 
                        help = 'output filtered bam file (defulat: - )')
    parser.add_argument('-x', '--index',
                        required=True, 
                        help = 'Mitochondrial gnome BWA index')
    args = parser.parse_args()
    return args


def main():
    logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
    logger = logging.getLogger('chrM Filter')
    args = getopt()
    chrM_aligner = BwaAligner(args.index)
    pair = 0
    out_pair =0
    with pysam.Samfile(args.inbam) as inbam:
        with pysam.Samfile(args.outbam, 'wb', template = inbam) as outbam:
            filter_bam = ''
            if args.filtered_bam:
                filter_bam = pysam.Samfile(args.filtered_bam, 'wb', 
                                            template = inbam)
            try:
                while True:
                    read1 = next(inbam)
                    read2 = next(inbam)
                    pair += 1

                    seq1 = read1.get_forward_sequence()
                    seq2 = read2.get_forward_sequence()

                    align1 = chrM_aligner.align_seq(seq1)
                    align2 = chrM_aligner.align_seq(seq2)

                    if align1 or align2:
                        if filter_bam:
                            filter_bam.write(read1)
                            filter_bam.write(read2)
                    else:
                        outbam.write(read1)
                        outbam.write(read2)
                        out_pair += 1
            except StopIteration:
                pass
    if filter_bam:
        filter_bam.close()
    
    logger.info('Read %i pairs, written %i pairs' %(pair, out_pair))


if __name__ == '__main__':
    main()

            
            
