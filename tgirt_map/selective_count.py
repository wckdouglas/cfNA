#!/usr/bin/env python

import fileinput
from itertools import groupby
import sys
from operator import itemgetter


def extract_fields(x):
    return itemgetter(0,1,2,3)(x.split('\t'))

def compute_overlap_score(read_start, read_end, ref_start, ref_end):
    overlapping_base = min(ref_end, read_end) -  max(ref_start, read_start)
    return overlapping_base/(read_end-read_start) * overlapping_base/(ref_end*ref_start)


def selective(read_start, read_end, lines):
    out_line = ''
    max_overlap = 0
    for line in lines:
        fields = line.split('\t')
        ref_start, ref_end = itemgetter(5,6)(fields)
        if ref_start == ref_end == "-1":
            return line.strip()

        else:
            overlapping_score = compute_overlap_score(int(read_start), int(read_end), 
                                        int(ref_start), int(ref_end))
            if overlapping_score > max_overlap:
                max_overlap = overlapping_score
                out_line = line
    return out_line.strip()


for (chrom, start, end, id), lines in groupby(fileinput.input(), extract_fields):
    out_line = selective(start, end, lines)
    print(out_line)

        

