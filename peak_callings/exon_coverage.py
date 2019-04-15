#!/usr/bin/env python

import pysam
import numpy as np
import os
from operator import itemgetter
import sys
from pybedtools import BedTool, set_tempdir
import pyBigWig as pbw
import glob
from sequencing_tools.io_tools import xopen
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
            bws = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/bed_files/merged_bed/coverage/unfragmented.{strand}.bigWig',
            exon_file = '/stor/work/Lambowitz/ref/hg19_ref/genes/exons.gencode.bed.gz',
            cutoff=2, 
            force=False):
        records = []
        self.high_cov_exons = '/stor/scratch/Lambowitz/cdw2854/high_cov_exon.bed'
        self.exon_file = exon_file
        self.bws = {strand:pbw.open(bws.format(strand=strand_label)) for strand, strand_label in zip(['-','+'],['rvs','fwd'])}
        if not os.path.isfile(self.high_cov_exons) or force:
            self.initiate(cutoff=cutoff)


    def initiate(self, cutoff=1):
        with xopen(self.exon_file) as exon_lines, \
                open(self.high_cov_exons, 'w') as outbed:
            for exon_line in exon_lines:
                fields = exon_line.split('\t')
                chrom, start, end, name, strand = itemgetter(0,1,2,3,5)(fields)
                if chrom.startswith('chr'):
                    try:
                        bw = self.bws[strand]
                        coverage_track = bw.values(chrom, int(start), int(end),numpy=True)
                    except RuntimeError:
                        print(chrom, start, end)
                        sys.exit()
                    uniform_coverage_score = len(coverage_track[coverage_track >= cutoff])/len(coverage_track)
                    avg_coverage = coverage_track.mean()
                    if avg_coverage > cutoff and uniform_coverage_score > 0.8:
                        outline = '{chrom}\t{start}\t{end}\t{name}\t{uniform_score}\t{strand}\t{avg_score}'\
                            .format(chrom = chrom,
                                    start = start, 
                                    end = end,
                                    name = name,
                                    strand = strand,
                                    uniform_score = uniform_coverage_score,
                                    avg_score = avg_coverage)
                        print(outline, file = outbed)
                    


    
    def filter(self, bed_df, f=0.2):
        in_count = bed_df.shape[0]
        columns = bed_df.columns
        out_df = BedTool()\
            .from_dataframe(bed_df)\
            .intersect(b = self.high_cov_exons, v=True, f=f, s=True)\
            .to_dataframe(names = columns)
        out_count = out_df.shape[0]
        print('Filtered out %i peaks' %(in_count - out_count))
        return out_df

        
def main():
    exon_filter = ExonFilter(force=True, cutoff=1)


if __name__ == '__main__':
    main()


        

