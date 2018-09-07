#!/usr/bin/env python

from sequencing_tools.fastq_tools.cutadapt_align import Aligner, locate
from sequencing_tools.fastq_tools import readfq
import glob
import sys


adapter = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTT'
adapter_rc = 'AAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
adapter_length = len(adapter)
aligner = Aligner(reference = adapter,
        max_error_rate = 0.1)
aligner2 = Aligner(reference = adapter_rc,
        max_error_rate = 0.1)

contam = 0
seq_count = 0
fs = glob.glob('/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/Q*/7sl_R1.fq')
cutoff = 3 #nt
for f in fs:
    with open(f) as fq:
        for fq_record in readfq(fq):
            seq_count += 1
            umi = fq_record.seq[:6]
            fq_record.subseq(6, len(fq_record.seq))
            matched = aligner.locate(fq_record.seq)
            matched2 = aligner2.locate(fq_record.seq)
            if matched and matched[1] == adapter_length and matched[4] > cutoff:
                refstart, refstop, querystart, querystop, matches, errors = matched
                if matched2 and matched2[1] == adapter_length and matched2[4] > cutoff:
                    refstart2, refstop2, querystart2, querystop2, matches2, errors2 = matched2
                    print('[{}] {}, {}, {}'.format(umi, fq_record.seq, matched, matched2))
                    print('[%s]%s[%s]' %(fq_record.seq[querystart:querystop],
                                         fq_record.seq[querystop:querystart2],
                                         fq_record.seq[querystart2:]))
                else:
                    print('[{}:{}] {}, {}'.format(fq_record.id, umi, fq_record.seq, matched))
                    print('[%s]%s' %(fq_record.seq[querystart:querystop],fq_record.seq[querystop:]))
                    contam += 1



print('{}/{} ({}%) has adapter > {}nt'.format(contam, seq_count, contam/seq_count*100, cutoff))
