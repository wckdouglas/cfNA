#!/usr/bin/env python

from pybedtools import BedTool
import pandas as pd
import dask.dataframe as dd
import dask
import os
import pysam
import re
import numpy as np
from operator import itemgetter
import glob
import sys
from functools import partial
import mappy 
from sequencing_tools.fastq_tools import reverse_complement
from multiprocessing import Pool


class PeakClassification():
    '''
    given peak bedfile, select best matches
    1. small RNA,
    2. RBP
    3. repeats
    4. overlap score
    '''
    def __init__(self, annotation_file=None, exon_table=None):
        # rank rna type
        self.preference_RNA = ['misc_RNA','srpRNA','SRP_RNA','vault_RNA','tRNA','snoRNA', 
                        'snRNA','scRNA','miRNA', 'rRNA','piRNA','RBP','SINE','SINE?',
                        'LINE','LINE?','Satellite','Simple_repeat',
                        'Low_complexity','LTR?','LTR', 'DNA','protein_coding','.']#,'RBP']
        self.preference_rank = {rna:i for i, rna in zip(range(100), self.preference_RNA)}
        self.max_rank = max(self.preference_rank.values()) + 1
        self.small_RNA = ['misc_RNA','srpRNA','SRP_RNA',
                    'vault_RNA','tRNA','snoRNA', 
                    'snRNA','scRNA','miRNA']

        # grouping RNA
        self.lrna_regex = 'lincR|protein|pseudo|TR|proces|sense_intr|[tT]elomer'\
                    '|prime|IG|antisen|lncRNA|sense_ov|TEC|RNase_[MP]|non_coding'   
        self.repeats_regex = 'LINE|Satellite|Simple_repeat|SINE|Unknown' \
                        '|Low_complexity|LTR|^DNA$|^DNA\?$|RC|Other'


        #genome filter
        self.annotation_file = annotation_file
        self.exon_table = exon_table


    def rank_type(self, rna):
        '''
        return rna rank
        '''
        try:
            return self.preference_rank[rna]
        except KeyError:
            return self.max_rank + 1


    def __merge_type__(self, x):
        '''
        grouping RNA type
        '''
        if re.search(self.repeats_regex, x):
            return 'Repeats'

        elif x == ".":
            return 'Unannotated'
        
        elif re.search(self.lrna_regex, x):
            return 'Long RNA'
            
        elif re.search('rRNA|rDNA', x):
            return 'rRNA'
        elif re.search('misc|guid|scRN|srpRNA|SRP_RNA|[vV]ault', x):
            return 'misc RNA'
        else:
            return x


    def annotate_peaks(self, bed):
        '''
        bedtools intersect annotation bed file and macs2 broad peaks
        adding overlap score and rna type rank
        '''
        columns = bed.columns.tolist()
        columns.extend(['gchrom','gstart','gend', 'gname','gscore','gstrand',
                        'gtype','gid','overlapped'])
        inbed = BedTool()\
            .from_dataframe(bed)\
            .intersect(wao=True, b=self.annotation_file) \
            .to_dataframe(names = columns) \
            .drop_duplicates() \
            .assign(strand = lambda d: np.where(d.peakname.str.contains('fwd'), '+','-')) \
            .assign(ref_overlap = lambda d: d.overlapped / (d.gend - d.gstart)) \
            .assign(peak_overlap = lambda d: d.overlapped / (d.end - d.start )) \
            .assign(overlap_score = lambda d: d.ref_overlap * d.peak_overlap)  \
            .assign(gtype = lambda d: np.where( (d.gname.str.contains('HY')) & (d.gtype=="scRNA"), 
                                                'misc_RNA',
                                                d.gtype)) \
            .assign(gtype = lambda d: np.where(d.gtype.str.contains('srpRNA'), 
                                                'misc_RNA',
                                                d.gtype)) \
            .assign(gtype = lambda d: np.where(d.gname.isin(["BC200",'7SK','7SL','VTRNA2-1']), 
                                                'misc_RNA', 
                                                d.gtype)) \
            .assign(gtype_rank = lambda d: d.gtype.map(self.rank_type))   \
            .drop(['gscore','peak_summit','gchrom','gstart', 'gend'], axis=1) \
            .fillna(0)
        
        print('Intersected genes ', file= sys.stderr)
        return inbed


    def resolve_annotation(self, inbed):
        '''
        select for greatest overlapped annotation
        '''
        #df = dd.from_pandas(inbed, npartitions=16, sort=True)\
        group_cols = ['chrom','start','end',
                        'peakname','score','is_sense', 
                        'fc','log10p',
                        'log10q','pileup','sample_count']
        df = inbed \
            .assign(is_sense = lambda d: np.where((d.strand == d.gstrand) | (d.gtype.str.contains(self.repeats_regex)), 
                                                'sense', 
                                                'antisense')) \
            .set_index('start')\
            .pipe(dd.from_pandas, npartitions=24) \
            .groupby(group_cols)\
            .apply(self.__select_annotation__, meta = {'gname':'f8',
                                                        'gtype':'f8',
                                                        'strand':'f8',
                                                        'gstrand':'f8'})\
            .compute(scheduler='multiprocessing')\
            .reset_index() \
            .drop_duplicates() \
            .sort_values('log10q', ascending=False)  \
            .assign(gtype = lambda d: d.gtype.map(self.__merge_type__)) \
            .drop('level_11', axis=1)  
    #        .pipe(retype_junctions, exon_table)
        
        sense_df = self.__strand_df__(df, strand = 'sense')
        antisense_df = self.__strand_df__(df, strand = 'antisense')

        print('Resolved genes ', file= sys.stderr)
        return sense_df.merge(antisense_df, how = 'outer').fillna('.') 

    def __strand_df__(self, df, strand = 'sense'):
        '''
        Select gene annotation for each peak in a sense-specific manner
        '''
        condition = 'is_sense == "sense"' if strand == 'sense' else 'is_sense == "antisense"'
        return  df\
            .query(condition) \
            .rename(columns={'gname': strand + '_gname',
                    'gtype': strand + '_gtype' })  \
            .drop(['is_sense','gstrand'],axis=1)


    def __select_annotation__(self, peak_rows):
        '''
        Finiding highest overlap and small rna for annotation
        '''
        max_overlapped = peak_rows\
            .pipe(lambda d: d[d.gtype.isin(self.small_RNA)])
        
        if max_overlapped.shape[0] == 0:
            max_overlapped = peak_rows
        
        max_overlapped = max_overlapped \
            .pipe(lambda d: d[d.overlap_score == d.overlap_score.max()]) \
            .filter(['gname','gtype','gtype_rank','gstrand','strand'])\
            .drop_duplicates() \
            .pipe(lambda d: d[d.gtype_rank == d.gtype_rank.min()])


        if max_overlapped.shape[0] > 1 and max_overlapped.gtype.unique()[0] == 'RBP':
            '''
            output all overlapping RBP
            '''
            max_overlapped = max_overlapped \
                .groupby(['gtype','strand','gstrand'], as_index=False)\
                .agg({'gname': lambda xs: ','.join(xs)})
        else: # or randomly pick one
            max_overlapped = max_overlapped.nsmallest(1, 'gtype_rank')
        
        required_columns = ['gname','gtype','strand', 'gstrand']
        df_dict = {col: max_overlapped[col] for col in required_columns}
        max_overlapped = pd.DataFrame(df_dict)
        return max_overlapped


def process_peaks(bed_path, peak_file):
    filename = os.path.basename(peak_file)
    samplename = filename.replace('.rvs_peaks.narrowPeak','').replace('.fwd_peaks.narrowPeak','')
    strand = re.findall('fwd|rvs', filename)[0]

    fragment_file = '%s/%s.%s.bed.gz' %(bed_path, samplename, strand)
    rows = []
    with pysam.Tabixfile(fragment_file) as tab, \
            open(peak_file,'r') as peaks:
        for peak in peaks:
            peak_fields = peak.strip().split('\t')
            peak_chrom, peak_start, peak_end = itemgetter(0,1,2)(peak_fields)
            peak_start, peak_end = int(peak_start), int(peak_end)
            coverage = np.zeros(peak_end - peak_start)
            sample_count = set()
            for fragments in tab.fetch(peak_chrom, peak_start, peak_end):
                frag_start, frag_end, frag_name = itemgetter(1,2,3)(fragments.split('\t'))
                frag_start, frag_end = int(frag_start), int(frag_end)
                frag_start = max(frag_start, peak_start)
                frag_end = min(frag_end, peak_end)
                coverage[(frag_start-peak_start):(frag_end-peak_start)] += 1
                sample_count.add(frag_name.split(':')[0])
            pileup = coverage.max()
            sample_count = min(len(sample_count),pileup)
            peak_fields.append(pileup)
            peak_fields.append(sample_count)
            rows.append(peak_fields)
    return pd.DataFrame(rows, columns = ['chrom','start','end',
                                    'peakname','score','strand','fc',
                                    'log10p','log10q','peak_summit','pileup','sample_count'])


def main():
    if len(sys.argv) < 3:
        sys.exit('[usage] python %s <out_table> <annotation_file> <bed_path> <exon_Table> [peak1] [peak2]' %sys.argv[0])

    out_table = sys.argv[1]
    annotation_file = sys.argv[2]
    bed_path = sys.argv[3]
    exon_table = sys.argv[4]
    peak_files = sys.argv[5:]
    print('Merging: ', ', '.join(map(os.path.basename, peak_files)))
    peak_classifier = PeakClassification(annotation_file = annotation_file,
                                        exon_table = exon_table)
    peak_processor = partial(process_peaks, bed_path)
    p = Pool(24)
    bed = pd.concat(p.map(peak_processor, peak_files)) \
        .sort_values(['chrom','start','end']) \
        .reset_index(drop=True) 
    p.close()
    p.join()

    #process annotation
    inbed = peak_classifier.annotate_peaks(bed) 
    df = peak_classifier.resolve_annotation(inbed)  
    df.to_csv(out_table, sep='\t', index=False)
    print('Written %s' %out_table)
    assert bed.shape[0] == df.shape[0], 'Peak lost!!!'


def make_table(base_name = 'unfragmented'):
    #test:
    all, base_name = False,  'unfragmented'
    project_path = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/bed_files/merged_bed'
    bed_path = project_path + '/stranded'
    peak_path = project_path + '/MACS2'
    annotated_path = peak_path + '/annotated'
    exon_table = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/merged_bam/unfragmented.spliced_exon.bed.gz'
    out_table = annotated_path + '/%s.annotated_peaks.tsv' %base_name
    annotation_file = os.environ['REF'] + '/hg19_ref/genes/all_annotation.bed.gz'


    if not os.path.isdir(annotated_path):
        os.mkdir(annotated_path)
    filtered = 'filtered' if not all else 'unfiltered'
    peak_files = glob.glob(peak_path + '/%s.%s.*_peaks.narrowPeak' %(base_name, filtered))
    print('Merging: ', ', '.join(map(os.path.basename, peak_files)))

    peak_classifier = PeakClassification(annotation_file = annotation_file,
                                        exon_table = exon_table)
    peak_processor = partial(process_peaks, bed_path)
    p = Pool(24)
    bed = pd.concat(p.map(peak_processor, peak_files)) \
        .sort_values(['chrom','start','end']) \
        .reset_index(drop=True) 
    p.close()
    p.join()

    inbed = peak_classifier.annotate_peaks(bed) 
    df = peak_classifier.resolve_annotation(inbed)  
    df.to_csv(out_table, sep='\t', index=False)
    print('Written %s' %out_table)
    assert bed.shape[0] == df.shape[0], 'Peak lost!!!'



if __name__ == '__main__':
    main()
