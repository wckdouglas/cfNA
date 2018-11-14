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
REF_PATH = os.environ['REF']
PROJECT_PATH = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/merged_bed'

class exon:
    def __init__(self, bed_line):
        self.fields = bed_line.strip().split('\t')
        self.chrom = self.fields[0]
        self.start = int(self.fields[1])
        self.end = int(self.fields[2])
        self.strand = self.fields[5]
        self.exon_name = self.fields[3]
        self.exon_size = self.end - self.start
        self.coverage_score = 0
        self.avg_coverage = 0
        self.extra = ','.join(self.fields[6:])

    def calculate_coverage(self, tabix_file, cutoff = 3):
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
            
        self.coverage_score = len(coverage_track[coverage_track >= cutoff])/self.exon_size
        self.avg_coverage = coverage_track.mean()
    

    def __str__(self):
        template = '{chrom}\t{start}\t{end}\t{exon}\t{score}\t{strand}\t{coverage}\t{info}'
        return template.format(chrom = self.chrom,
                        start = self.start,
                        end = self.end,
                        exon = self.exon_name,
                        score = self.coverage_score,
                        strand = self.strand,
                        coverage = self.avg_coverage,
                        info = self.extra)


def make_exons(tab_file, cov_exon, exons):
    out_exon_count = 0
    with open(exons, 'r') as exon_records, \
            pysam.Tabixfile(tab_file) as tabix, \
            open(cov_exon, 'w') as out_exon:
        for in_exon_count, exon_record in enumerate(exon_records):
            ex = exon(exon_record)
            ex.calculate_coverage(tabix, cutoff = 3)
            if ex.coverage_score > 0.8 and ex.avg_coverage > 3:
                print(str(ex), file = out_exon)
                out_exon_count += 1
    print('Read %i exons, written %i exons' %(in_exon_count, out_exon_count), 
          file=sys.stdout)


def filter_bed(tab_file, out_prefix, cov_exon):
    
    positive_out = out_prefix + '.fwd.bed'
    negative_out = out_prefix + '.rvs.bed'
    print('Writing:\n1. %s\n2. %s' %(positive_out, negative_out))
    set_tempdir(os.path.dirname(out_prefix))


    bed_filters = [REF_PATH + '/hg19/new_genes/tRNA_yRNA.bed',
        REF_PATH + '/hg19/genome/tRNA.bed',
        REF_PATH + '/hg19/new_genes/sncRNA_x_protein.bed',
        REF_PATH + '/hg19/new_genes/dashr_tRNA.bed',
        REF_PATH + '/hg19/new_genes/rmsk_tRNA.bed',
        REF_PATH + '/hg19/new_genes/rmsk_rRNA.bed',
        REF_PATH + '/hg19/new_genes/rmsk_yRNA.bed',
        REF_PATH + '/hg19/new_genes/refseq_rRNA.bed',
        REF_PATH + '/hg19/new_genes/rRNA_for_bam_filter.bed']

    _filtered = BedTool(tab_file) \
        .intersect(b = bed_filters , v=True)\
        .saveas()
   
    _filtered = _filtered \
        .intersect(b = cov_exon, s=True, v=True) \
        .saveas()

    with open(positive_out, 'w') as pos,\
            open(negative_out, 'w') as neg:
        for fragment_count, fragment in enumerate(_filtered.filter(less_than, 300)):
            out_file = neg if fragment.strand == '-' else pos
            out_file.write(str(fragment))

    print('Output %i fragments' %fragment_count)
    for out in [negative_out, positive_out]:
        os.system('bgzip -f {out}; tabix -f -p bed {out}.gz'.format(out = out))
    return 0


def main():
    if len(sys.argv) != 3:
        sys.exit('[usage] python %s <bed_file> <out_prefix>' %sys.argv[0])

    exons = REF_PATH + '/hg19/new_genes/exons.bed'
    tab_file = sys.argv[1]
    out_prefix =  sys.argv[2]

    prefix = os.path.basename(tab_file).split('.')[0]
    cov_exon = out_prefix + '_exons.bed'

    make_exons(tab_file, cov_exon, exons)
    filter_bed(tab_file, out_prefix, cov_exon)

def main1():
    tab_files = glob.glob(PROJECT_PATH + '/*.bed.gz')
    tab_files = filter(lambda x: re.search('all|unfragmented', os.path.basename(x)), tab_files)
    tab_files = filter(lambda x: 'no_' not in os.path.basename(x), tab_files)

    processor = partial(bed_processor, exons)
    p = Pool(24)
    p.map(processor, tab_files)
    p.close()
    p.join()

    set_tempdir(PROJECT_PATH + '/stranded')
    cleanup()

if __name__ == '__main__':
    main()
