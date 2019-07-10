#!/usr/bin/env python


from sequencing_tools.fastq_tools import read_interleaved
import sys


for r1, r2 in read_interleaved(sys.stdin):
    if min(len(r1.seq), len(r2.seq)) > 50:
        print(r1)
        print(r2)



