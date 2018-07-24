#!/usr/bin/env python

from sequencing_tools.fastq_tools.cutadapt_align import Aligner
from sequencing_tools.fastq_tools import readfq
import glob


adapter = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
adapter_length = len(adapter)
aligner = Aligner(reference = adapter,
        max_error_rate = 0.1)

contam = 0
seq_count = 0
fs = glob.glob('/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/Qcf*/7sl_R1.fq')
cutoff = 3 #nt
for f in fs:
    with open(f) as fq:
        for fq_record in readfq(fq):
            seq_count += 1
            umi = fq_record.seq[:6]
            fq_record.subseq(6, len(fq_record.seq))
            matched = aligner.locate(fq_record.seq)
            if matched and matched[1] == adapter_length and matched[4] > cutoff:
                refstart, refstop, querystart, querystop, matches, errors = matched
                print('[{}] {}, {}'.format(umi, fq_record.seq, matched))
                print('[%s]%s' %(fq_record.seq[querystart:querystop],fq_record.seq[querystop:]))
                contam += 1

print('{}/{} ({}%) has adapter > {}nt'.format(contam, seq_count, contam/seq_count*100, cutoff))
