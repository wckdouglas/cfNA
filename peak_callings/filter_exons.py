#!/usr/bin/env python

import pandas as pd
from pybedtools import BedTool
import pysam


class mRNAFilter():
    def __init__(self):
        ref_path = '/stor/work/Lambowitz/ref/hg19_ref/genes'
        exons = ref_path + '/gencode.exon.bed.gz'
        self.exons = pysam.Tabixfile(exons)
        transcriptom_peaks = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/transcriptome/macs2/unfragmented.fwd_peaks_genomics.narrowPeak.gz'
        self.transcriptome_peaks = pysam.Tabixfile(transcriptom_peaks)

    def search(self, chrom, start, end, attribute = 'exon'):
        if attribute == 'exon':
            it = self.exons
        elif attribute == 'transcriptome':
            it = self.transcriptome_peaks
        return 'yes' if any(it.fetch(chrom, start, end)) else 'no'


m_filter = mRNAFilter()
peaks = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/bed_files/merged_bed/MACS2_unfiltered/annotated/unfragmented.filtered.tsv'
peaks = pd.read_csv(peaks, sep='\t')\
    .query('pileup >= 5 & sample_count >= 5') \
    .assign(is_exon = lambda d: [m_filter.search(chrom, start, end, attribute = 'exon') for chrom, start, end in zip(d.chrom, d.start, d.end)]) \
    .assign(is_transcriptome_peak = lambda d: [m_filter.search(chrom, start, end, attribute = 'transcriptome') for chrom, start, end in zip(d.chrom, d.start, d.end)])

