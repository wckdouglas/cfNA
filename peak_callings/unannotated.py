#!/usr/bin/env python

import pandas as pd
import pysam
from bwapy import BwaAligner



class chrM_mapper():
    def __init__(self, bam, index):
        self.bam = pysam.Samfile(bam)
        self.aligner = BwaAligner(index, options = '-k 12')

    def __align__(self, seq):
        alignments = self.aligner.align_seq(seq)
        return 1 if alignments else 0

    def run_peak(self, chrom, start, end, peak_strand):
        alns = 0
        chrM_alns = 0
        for aln in self.bam.fetch(chrom, start, end):
            if aln.is_read1:
                aln_strand = '-' if aln.is_reverse  else '+'
                if aln_strand == peak_strand:
                    alns += 1
                    seq = aln.get_forward_sequence()
                    chrM_alns += self.__align__(seq)
        return alns, chrM_alns



project_path = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
peaks = project_path + '/bed_files/merged_bed/MACS2/annotated/unannotated.feather'
bam = project_path + '/merged_bam/unfragmented.bam'
chrM_index = '/stor/work/Lambowitz/ref/hg19/genome/chrM.fa'
peaks = pd.read_feather(peaks)


chrM_tester = chrM_mapper(bam, chrM_index)
for i, row in peaks.iterrows():
    alns, chrM_alns = chrM_tester.run_peak(row['chrom'], row['start'], row['end'], row['strand'])
    print(row['chrom'], row['start'], row['end'], row['strand'], 
        row['antisense_gtype'], row['antisense_gname'],
        alns, chrM_alns, round(chrM_alns/alns*100,3))