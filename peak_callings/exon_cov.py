#!/usr/bin/env python

import pysam
import numpy as np
import os
from operator import itemgetter
import sys
from pybedtools import BedTool, set_tempdir
from pybedtools.featurefuncs import less_than
from pybedtools.helpers import cleanup
from sequencing_tools.io_tools import xopen
import glob
from multiprocessing import Pool
from functools import partial
import re
from itertools import groupby
set_tempdir('/stor/scratch/Lambowitz/cdw2854')

class exon:
    def __init__(self, bed_line):
        self.fields = bed_line.strip().split('\t')
        self.chrom = self.fields[0]
        self.start = int(self.fields[1])
        self.end = int(self.fields[2])
        self.exon_count = self.fields[4]
        self.strand = self.fields[5]
        self.exon_name = self.fields[3]
        self.exon_size = self.end - self.start
        self.uniform_coverage_score = 0
        self.avg_coverage = 0
        self.extra = ','.join(self.fields[6:])

    def calculate_coverage(self, tabix_file, cutoff = 2):
        coverage_track = np.zeros(self.exon_size)
        try: 
            for fragment in tabix_file.fetch(self.chrom, self.start, self.end):
                f = fragment.strip().split('\t')
                f_start, f_end, f_strand = itemgetter(1,2,5)(f)

                if f_strand == self.strand:
                    adjusted_start = int(f_start) - self.start
                    adjusted_end = int(f_end) - self.start

                    adjusted_start = 0 if adjusted_start < 0 else adjusted_start
                    adjusted_end = self.exon_size - 1 if adjusted_end > self.exon_size - 1 else adjusted_end

                    coverage_track[adjusted_start:adjusted_end] += 1
        except ValueError:
            pass
            
        self.uniform_coverage_score = len(coverage_track[coverage_track >= cutoff])/self.exon_size
        self.avg_coverage = coverage_track.mean()

    def __str__(self):
        template = '{chrom}\t{start}\t{end}\t{exon}\t{score}\t{strand}\t{coverage}\t{info}'
        return template.format(chrom = self.chrom,
                        start = self.start,
                        end = self.end,
                        exon = self.exon_name,
                        score = self.uniform_coverage_score,
                        strand = self.strand,
                        coverage = self.avg_coverage,
                        info = self.extra)
    
class ExonFilter():
    def __init__(self,
            tab_file = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/bed_files/merged_bed/unfragmented.bed.gz',
            exon_file = '/stor/work/Lambowitz/ref/hg19_ref/genes/exons.gencode.bed.gz',
            cutoff=2):
        records = []
        self.high_cov_exons = '/stor/scratch/Lambowitz/cdw2854/high_cov_exon.bed'
        self.exon_file = exon_file
        self.tab_file = tab_file
        if not os.path.isfile(self.high_cov_exons):
            self.initiate(cutoff=cutoff)


    def __groupby__(self, frag):
        return frag.chrom, frag.start, frag.end, frag.name, frag.strand

    def initiate(self, cutoff=2):
        cdef: 
            str chrom, start
        coverages = BedTool(self.exon_file)\
            .coverage(self.tab_file, s=True, d=True, stream=True)
        with open(self.high_cov_exons, 'w') as outbed:
            for (chrom,start,end,name,strand), positions in groupby(coverages, self.__groupby__):
                coverage_track = np.array([position.fields[-1] for position in positions], dtype='int')
                self.uniform_coverage_score = len(coverage_track[coverage_track >= cutoff])/(end-start)
                self.avg_coverage = coverage_track.mean()
                if self.avg_coverage > cutoff and self.uniform_coverage_score > 0.8:
                    outline = '{chrom}\t{start}\t{end}\t{name}\t0\t{strand}'\
                            .format(chrom = chrom,
                                    start = start, 
                                    end = end,
                                    name = name,
                                    strand = strand)
                    print(outline, file = outbed)
    
    def filter(self, bed_df):
        in_count = bed_df.shape[0]
        columns = bed_df.columns
        out_df = BedTool()\
            .from_dataframe(bed_df)\
            .intersect(b = self.high_cov_exons, v=True, f=0.6, s=True)\
            .to_dataframe(names = columns)
        out_count = out_df.shape[0]
        print('Filtered out %i peaks' %(in_count - out_count))

        

        

