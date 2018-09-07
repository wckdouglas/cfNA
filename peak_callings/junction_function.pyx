from __future__ import print_function
from builtins import range, map
import pandas as pd
import pysam
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
from collections import Counter
from sequencing_tools.bam_tools import make_regions, check_concordant, get_strand
import os
import sys
import re


cigar_op = re.compile('[MIDNSH]')
cigar_nums = re.compile('[0-9]+')
cpdef set get_junction(aln):

    cdef:
        long parsed = aln.pos
        str op, nt, strand, chrom
        long junc_start, junc_end
        set intron_list = set()

    chrom = aln.reference_name
    strand = get_strand(aln)

    cigar_ops = cigar_op.findall(aln.cigarstring)
    cigar_bases = cigar_nums.findall(aln.cigarstring)

    for op, nt in zip(cigar_ops, cigar_bases):
        if op in ['M','D']: # only M and D is related to genome position
            parsed += int(nt)
        elif op == 'N':
            junc_start = parsed
            parsed += int(nt)
            junc_end = parsed + 1 # 1-pos adjustment to match get_align_pairs
            junction = chrom + ':' + str(junc_start) + '-' + str(junc_end) + '_' + strand
            intron_list.add(junction)
    return intron_list


cpdef set get_splice(AlignedSegment aln):
    cdef:
        bint is_splice
        set junctions = set()

    is_splice = 'N' in aln.cigarstring
    if is_splice and not aln.is_duplicate and not aln.is_supplementary:
        junctions = get_junction(aln)
    return junctions




def find_junction(bam_file):

    cdef:
        AlignedSegment aln
        set aln_introns
        long aln_count = 0
        str intron 

    intron_counter = Counter()
    samplename = os.path.basename(bam_file)
    print('Running %s' %samplename)
    with pysam.Samfile(bam_file, 'r') as bam:
        for aln in bam:
            aln_introns = get_splice(aln)

            if aln_introns:
                for intron in aln_introns:
                    intron_counter[intron] += 1
                    aln_count += 1
                    if aln_count % 10000000 == 0:
                        print('[%s] Parsed %i pairs' %(samplename, aln_count))
    return pd.DataFrame({'intron': list(intron_counter.keys()),
                    'intron_count': list(intron_counter.values())})
