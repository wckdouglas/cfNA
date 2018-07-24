#!/usr/bin/env python

from pybedtools import BedTool
import pandas as pd
import os
import pysam
import re
import numpy as np
from operator import itemgetter
import glob
import sys


def process_broad(broad_peak, bed_path):
    filename = os.path.basename(broad_peak)
    samplename = filename.split('.')[0]
    strand = re.findall('fwd|rvs', filename)[0]

    fragment_file = '%s/%s.%s.bed.gz' %(bed_path, samplename, strand)
    rows = []
    with pysam.Tabixfile(fragment_file) as tab, \
            open(broad_peak,'r') as peaks:
        for peak in peaks:
            peak_fields = peak.strip().split('\t')
            peak_chrom, peak_start, peak_end = itemgetter(0,1,2)(peak_fields)
            peak_start, peak_end = int(peak_start), int(peak_end)
            coverage = np.zeros(peak_end - peak_start)
            for fragments in tab.fetch(peak_chrom, peak_start, peak_end):
                frag_start, frag_end = itemgetter(1,2)(fragments.split('\t'))
                frag_start, frag_end = int(frag_start), int(frag_end)
                frag_start = max(frag_start, peak_start)
                frag_end = min(frag_end, peak_end)
                coverage[(frag_start-peak_start):(frag_end-peak_start)] += 1
            pileup = coverage.max()
            peak_fields.append(pileup)
            rows.append(peak_fields)
    return pd.DataFrame(rows)

def collapse_col(x):
    return ','.join(x)


def make_table(HEPG2=None):
    project_path = '/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/merged_bed'
    bed_path = project_path + '/stranded'
    peak_path = project_path + '/MACS2'
    annotated_path = peak_path + '/annotated'
    
    if not HEPG2:
        out_table = annotated_path + '/unfragmented.annotated_peaks.tsv'
        annotation_file = os.environ['REF'] + '/hg19/new_genes/all_annotation.bed.gz'
    else:
        out_table = annotated_path + '/unfragmented.annotated_peaks_hepg2.tsv'
        annotation_file = os.environ['REF'] + '/hg19/new_genes/all_annotation_hepG2.bed.gz'


    if not os.path.isdir(annotated_path):
        os.mkdir(annotated_path)

    broad_peaks = glob.glob(peak_path + '/unfrag*_peaks.broadPeak')

    bed = pd.concat([process_broad(broad_peak, bed_path) for broad_peak in broad_peaks]) \
        .sort_values([0,1,2])
    inbed = BedTool()\
        .from_dataframe(bed)\
        .intersect(wao=True, b=annotation_file) \
        .to_dataframe(names = ['chrom','start','end',
                                'peakname','score','strand','fc',
                                'log10p','log10q','pileup',
                                'gname','gstrand','gtype','gid'],
                      usecols = [0,1,2,3,4,5,6,7,8,9, 13, 15, 16, 17]) \
        .assign(strand = lambda d: np.where(d.peakname.str.contains('fwd'), '+','-'))\
        .groupby(['chrom','start','end',
                    'peakname','score','strand', 
                    'fc','log10p',
                    'log10q','pileup'], 
                 as_index=False)\
        .agg({'gname': collapse_col,
            'gstrand': collapse_col,
            'gtype': collapse_col,
            'gid': collapse_col}) \
        .sort_values('log10q', ascending=False)\
        .to_csv(out_table, sep='\t', index=False)
    print('Written %s' %out_table)

def main():
    HEPG2s=[True, False]
    [make_table(HEPG2) for HEPG2 in HEPG2s]

if __name__ == '__main__':
    main()
