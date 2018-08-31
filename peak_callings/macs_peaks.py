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



lrna_regex = 'lincR|protein|pseudo|TR|proces|sense_intr'\
            'prime|IG|antisen|lncRNA|sense_ov|TEC'   
def merge_type(x):
    if re.search('LINE|Satellite|Simple_repeat|SINE|Unknown'
                 '|Low_complexity|LTR|^DNA$|^DNA\?$|RC|Other', x):
        return 'Repeats'

    elif x == ".":
        return 'Unannotated'
    
    elif re.search(lrna_regex, x):
        return 'Long RNA'
        
    elif re.search('rRNA|rDNA', x):
        return 'rRNA'
    elif re.search('misc|guid|scRN|srpRNA', x):
        return 'misc RNA'
    else:
        return x

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


def make_table_old(all=None, base_name = 'unfragmented'):
    project_path = '/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/merged_bed'
    bed_path = project_path + '/stranded'
    peak_path = project_path + '/MACS2'
    annotated_path = peak_path + '/annotated'
    
    if all:
        out_table = annotated_path + '/%s.annotated_peaks.tsv' %base_name
        annotation_file = os.environ['REF'] + '/hg19/new_genes/all_annotation.bed.gz'
    else:
        out_table = annotated_path + '/%s.annotated_peaks_k562.tsv' %base_name
        annotation_file = os.environ['REF'] + '/hg19/new_genes/all_annotation_k562.bed.gz'


    if not os.path.isdir(annotated_path):
        os.mkdir(annotated_path)

    broad_peaks = glob.glob(peak_path + '/%s*_peaks.broadPeak' %base_name)

    bed = pd.concat([process_broad(broad_peak, bed_path) for broad_peak in broad_peaks], sort=False) \
        .sort_values([0,1,2]) \
        .reindex()

    inbed = BedTool()\
        .from_dataframe(bed)\
        .intersect(wao=True, b=annotation_file) \
        .to_dataframe(names = ['chrom','start','end',
                                'peakname','score','strand','fc',
                                'log10p','log10q','pileup','gstart','gend',
                                'gname','gstrand','gtype','gid','overlapped'],
                      usecols = [0,1,2,3,4,5,6,7,8,9,11,12, 13, 15, 16, 17, 18]) \
        .drop_duplicates() \
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


small_RNA = ['misc_RNA','srpRNA','SRP_RNA',
            'vault_RNA','tRNA','snoRNA', 
            'snRNA','scRNA','miRNA']
def select_annotation(peak_rows):
    '''
    Finiding highest overlap and small rna for annotation
    '''
    max_overlapped = peak_rows\
        .pipe(lambda d: d[d.gtype.isin(small_RNA)])
    
    if max_overlapped.shape[0] == 0:
        max_overlapped = peak_rows
       
    max_overlapped = max_overlapped\
        .pipe(lambda d: d[d.overlap_score == d.overlap_score.max()]) 

    
    
    required_columns = ['gname','gtype']
    df_dict = {col: [collapse_col(max_overlapped[col])] for col in required_columns}
    max_overlapped = pd.DataFrame(df_dict)
    return max_overlapped


def strand_df(df, strand = 'sense'):
    '''
    Select gene annotation for each peak in a sense-specific manner
    '''
    condition = 'gstrand == strand' if strand == 'sense' else 'gstrand != strand'
    return  df\
        .query(condition) \
        .rename(columns={'gname': strand + '_gname',
                  'gtype': strand + '_gtype' })  \
        .drop('gstrand',axis=1)
        



def resolve_annotation(inbed):
    '''
    select for greatest overlapped annotation
    '''
    df = dd.from_pandas(inbed, npartitions=16, sort=True)\
        .groupby(['chrom','start','end',
                    'peakname','score','strand', 
                    'fc','log10p',
                    'log10q','pileup','gstrand'])\
        .apply(select_annotation)

    df = df.compute(scheduler='processes',
                      get = dask.multiprocessing.get) \
        .reset_index()\
        .drop_duplicates() \
        .sort_values('log10q', ascending=False) 
    
    sense_df = strand_df(df, strand = 'sense')
    antisense_df = strand_df(df, strand = 'antisense')

    return sense_df.merge(antisense_df, how = 'outer') 


def annotate_peaks(annotation_file, bed):
    '''
    bedtools intersect annotation bed file and macs2 broad peaks
    '''
    inbed = BedTool()\
        .from_dataframe(bed)\
        .intersect(wao=True, b=annotation_file) \
        .to_dataframe(names = ['chrom','start','end',
                                'peakname','score','strand','fc',
                                'log10p','log10q','pileup','gstart','gend',
                                'gname','gstrand','gtype','gid','overlapped'],
                      usecols = [0,1,2,3,4,5,6,7,8,9,11,12, 13, 15, 16, 17, 18]) \
        .drop_duplicates() \
        .assign(strand = lambda d: np.where(d.peakname.str.contains('fwd'), '+','-')) \
        .assign(ref_overlap = lambda d: d.overlapped / (d.gend - d.gstart)) \
        .assign(peak_overlap = lambda d: d.overlapped / (d.end - d.start )) \
        .assign(overlap_score = lambda d: d.ref_overlap * d.peak_overlap)  \
        .assign(gtype = lambda d: np.where( (d.gname.str.contains('HY')) & (d.gtype=="scRNA"), 
                                            'misc RNA',
                                            d.gtype)) \
        .assign(gtype = lambda d: np.where(d.gtype.str.contains('srpRNA'), 
                                             'misc RNA',
                                             d.gtype)) \
        .assign(gtype = lambda d: np.where(d.gname == "BC200", 
                                            'misc RNA', 
                                            d.gtype)) \
        .assign(gtype = lambda d: d.gtype.map(merge_type))
    
    return inbed


def make_table(all=None, base_name = 'unfragmented'):
    project_path = '/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/merged_bed'
    bed_path = project_path + '/stranded'
    peak_path = project_path + '/MACS2'
    annotated_path = peak_path + '/annotated'
    
    if all:
        out_table = annotated_path + '/%s.annotated_peaks.tsv' %base_name
        annotation_file = os.environ['REF'] + '/hg19/new_genes/all_annotation.bed.gz'
    else:
        out_table = annotated_path + '/%s.annotated_peaks_k562.tsv' %base_name
        annotation_file = os.environ['REF'] + '/hg19/new_genes/all_annotation_k562.bed.gz'


    if not os.path.isdir(annotated_path):
        os.mkdir(annotated_path)

    broad_peaks = glob.glob(peak_path + '/%s*_peaks.broadPeak' %base_name)

    bed = pd.concat([process_broad(broad_peak, bed_path) for broad_peak in broad_peaks], sort=False) \
        .sort_values([0,1,2]) \
        .reindex()
    

    inbed = annotate_peaks(annotation_file, bed)

    inbed = resolve_annotation(inbed)\
        .drop('level_11', axis=1)\
        .to_csv(out_table, sep='\t', index=False)
    print('Written %s' %out_table)


def main():
    all_annotations=[True, False]
    base_names = ['unfragmented', 'exonuclease'] 
    for all in all_annotations:
        for base_name in base_names:
            make_table(all=all, base_name = base_name) 

if __name__ == '__main__':
    main()
