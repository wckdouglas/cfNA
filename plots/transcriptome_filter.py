#!/usr/bin/env python

import mappy as mp
import argparse
import pysam
import numpy as np
import sys
from operator import itemgetter
from collections import defaultdict, Counter
from sequencing_tools.bam_tools import paired_bam, get_strand


class peak_analyzer():
    def __init__(self, bam, idx):
        self.bam = pysam.Samfile(bam)
        self.aligner = mp.Aligner(idx, preset='sr')
    
    def __get_alignments__(self, chrom, start, end):
        '''
        pre-paired alignments from coordinate-sorted bam
        '''
        paired_alignments = defaultdict(lambda: defaultdict(pysam.AlignedRead))
        for aln in self.bam.fetch(str(chrom), int(start), int(end)):
            read = 'read1' if aln.is_read1 else 'read2'
            paired_alignments[aln.query_name][read] = aln
        return paired_alignments

    def __align_seq__(self, seq1, seq2):
        '''
        using minimap to align paired end sequence and return True if concordant pairs are formed
        '''
        alns = self.aligner.map(seq1, seq2)
        paired_alignments = defaultdict(lambda: defaultdict(int))
        for aln in alns:
            paired_alignments[aln.ctg][aln.read_num] = aln
        return paired_alignments
        


    def __filter__(self, paired_alignments):            
        for transcript, transcript_pairs in paired_alignments.items():
            if transcript_pairs[1] and transcript_pairs[2]:
                #if transcript_pairs[1].trans_strand != transcript_pairs[2].trans_strand:
                    return 1
        return 0
        

    def pick_best_transcript(self, transcript_counter):
        max_tc = 0
        max_tr = ''
        for tr, tc in transcript_counter.items():
            if tc > max_tc:
                max_tr = tr
        return max_tr


    def filter_alignments(self, chrom, start, end, strand):
        paired_alignments = self.__get_alignments__(chrom, start, end)
        mapped = 0
        num_pairs = 0
        transcripts = Counter()
        for num_pairs, (query_name, pairs) in enumerate(paired_alignments.items()):
            if pairs['read1'].query_name and pairs['read2'].query_name:
                if  get_strand(pairs['read1']) == get_strand(pairs['read2']) == strand:
                    seq1 = pairs['read1'].get_forward_sequence()
                    seq2 = pairs['read2'].get_forward_sequence()
                    paired = self.__align_seq__(seq1, seq2)
                    mapped += self.__filter__(paired)
                    for tr in paired.keys():
                        transcripts[tr] += 1
        transcript = self.pick_best_transcript(transcripts)

        return mapped, num_pairs + 1,  transcript
        
        


def getopt():
    parser = argparse.ArgumentParser(description = 'MApping fastq to minimap2 index and filter out ')
    parser.add_argument('--idx', required=True,
                        help = 'minimap2 index' )
    parser.add_argument('--bam', required=True,
                        help = 'paired-end alignment bam file')
    parser.add_argument('--bed', default='-', help = 'peak file (default stdin)')

    return parser.parse_args()


def main():
    args = getopt()
    #idx = '/stor/work/Lambowitz/ref/hg19_ref/genes/transcriptome.minimap2_idx'
    #bam = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/merged_bam/dedup/unfragmented.chrM_filter.dedup.bam'
    #chrom = 'chr16' #chr16:31173304
    #start = 31173304
    #end = 31173390
    #strand = '+'
    

    PA = peak_analyzer(args.bam, args.idx)
    infile = sys.stdin if args.bed == "-" else open(args.bed)
    for line in infile:
        line = line.strip()
        fields = line.split('\t')
        chrom, start, end, strand = itemgetter(0,1,2,5)(fields)
        mapped, num_pairs, transcript = PA.filter_alignments(chrom, int(start), int(end), strand)
        fraction = mapped/num_pairs if num_pairs > 0 else 0
        print(line + '%i\t%i\t%.3f\t%s' %(mapped,num_pairs, fraction, transcript))


if __name__ == '__main__':
    main()