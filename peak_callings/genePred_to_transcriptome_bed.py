#!/usr/bin/env python

import fileinput
import re
from operator import itemgetter
import numpy as np

'''
for each bed12 record, find the exon ends and print exon junction in bed
'''

def parse_intervals(ex_ints):
    return [int(ex) for ex in ex_ints.strip(',').split(',') ]

for transcript in fileinput.input():
    pos = 0
    transcript = transcript.strip()
    fields = transcript.split('\t')
    tid, chrom, strand, tss, tes, exon_count, ex_starts, ex_ends = itemgetter(0, 1, 2, 3, 4, 7,8, 9)(fields)

    tss = int(tss)
    ex_starts = parse_intervals(ex_starts)
    ex_ends = parse_intervals(ex_ends)

    assert len(ex_starts) == len(ex_ends) == int(exon_count)

    tx_starts, tx_ends = (ex_starts, ex_ends) if strand == "+" else (ex_starts[::-1], ex_ends[::-1])
    for s, e in zip(ex_starts[:-1], ex_ends[:-1]):
        # remove last exon end
        exon_size = e - s
        exon_end = pos + exon_size
        genome_end = e if strand == '+' else s
        print('{tid}\t{junction_site1}\t{junction_site2}\t{tid}\t0\t+\t{coor}'\
                .format(tid = tid, 
                        junction_site1 = exon_end - 1,
                        junction_site2 = exon_end + 1,
                        coor = chrom + ':' + str(genome_end)))
        pos += exon_size
