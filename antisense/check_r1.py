#!/usr/bin/env python

import glob
import sys
from functools import partial
from sequencing_tools.fastq_tools.cutadapt_align import Aligner, locate
from sequencing_tools.fastq_tools import readfq
from sequencing_tools.io_tools import xopen

def process_read(aligner, adapter_length, nt_cutoff, out_file, fq_record):
    contam = 0
    umi = fq_record.seq[:6]
    fq_record.subseq(6, len(fq_record.seq))
    matched = aligner.locate(fq_record.seq)
    if matched:
        positional_match = (matched[1] == adapter_length or matched[1] == (adapter_length -1))
        match_enough = matched[4] > nt_cutoff
        if positional_match  and match_enough:
            refstart, refstop, querystart, querystop, matches, errors = matched
            print('[{}:{}] {}, {}'.format(fq_record.id, umi, fq_record.seq, matched), file=out_file)
            print('[%s]%s' %(fq_record.seq[querystart:querystop],fq_record.seq[querystop:]), file=out_file)
            contam = 1
    return contam


def find_r2(fqfile, out_file=None, nt_cutoff=6):
    '''
    This module tries to find R2R contaminations in read 1 sequnecese

    1st template switch:

    5' ~~~~~~~~~~~~        (RNA template)
    3' <~~~~~~~~~(N)R2R 5' (new)

    2nd template switch:


    5' R2R------------>R2  3' (new)
    3'    <~~~~~~~~~(N)R2R 5' (bottom strand from last step)


    End product after ligation:

        R1R          UMI         Insert               R2R     
    >>>>>>>>>>>>>>>>>>XXXXXX--(R2)--++++(cDNA)+++++<<<<<<<<<<<<<<<


    usage:
    input:
        fqfile:  filename for Read1 fastq file
        out_file: saved matched reads
        nt_cutoff:  how many nucleotide matches as cutoff
    
    output:
        seq_count:     sequence count
        contam_count:  sequence count with contamination
    '''

    R2R = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTT'
    adapter_length = len(R2R)
    aligner = Aligner(reference = R2R,
            max_error_rate = 0.1)

    contam = 0
    seq_count = 0
    with xopen(fqfile) as fq, open(out_file, 'w') as out:
        read_proccessor = partial(process_read, aligner, adapter_length, nt_cutoff, out)
        for fq_record in readfq(fq):
            seq_count += 1
            contam += read_proccessor(fq_record)
    return seq_count, contam
