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

preference_RNA = ['misc_RNA','srpRNA','SRP_RNA','vault_RNA','tRNA','snoRNA', 
                  'snRNA','scRNA','miRNA', 'rRNA','piRNA','RBP','SINE','SINE?',
                  'LINE','LINE?','Satellite','Simple_repeat',
                 'Low_complexity','LTR?','LTR', 'DNA','protein_coding','.']#,'RBP']
preference_rank = {rna:i for i, rna in zip(range(100), preference_RNA)}
max_rank = max(preference_rank.values()) + 1
def rank_type(rna):
    try:
        return preference_rank[rna]
    except KeyError:
        return max_rank + 1


lrna_regex = 'lincR|protein|pseudo|TR|proces|sense_intr|[tT]elomer'\
            '|prime|IG|antisen|lncRNA|sense_ov|TEC|RNase_[MP]|non_coding'   
repeats_regex = 'LINE|Satellite|Simple_repeat|SINE|Unknown' \
                 '|Low_complexity|LTR|^DNA$|^DNA\?$|RC|Other'
def merge_type(x):
    if re.search(repeats_regex, x):
        return 'Repeats'

    elif x == ".":
        return 'Unannotated'
    
    elif re.search(lrna_regex, x):
        return 'Long RNA'
        
    elif re.search('rRNA|rDNA', x):
        return 'rRNA'
    elif re.search('misc|guid|scRN|srpRNA|SRP_RNA|[vV]ault', x):
        return 'misc RNA'
    else:
        return x


def is_mt(seq, rnr=False):
    is_chrM = 'not_MT'
    chrom_path = '/stor/work/Lambowitz/ref/hg19'
    if rnr:
        genome = chrom_path + '/new_genes/mt_rnr.fa'
    else:
        genome = chrom_path + '/genome/chrM.minimap2_idx'

    aligner = mappy.Aligner(genome,preset='sr')
    if list(aligner.map(seq)):
        is_chrM = 'is_MT'
    return is_chrM

fa = pysam.Fastafile('/stor/work/Lambowitz/ref/hg19/genome/hg19_genome.fa')
def fetch_seq(chrom, start, end, strand):
    seq = fa.fetch(chrom, int(start), int(end))
    seq = seq.upper()
    return seq if strand == "+" else reverse_complement(seq)



def is_junction_exon(junctions, chrom,start, end, strand):
    '''
    junctions are splice sites from bam file
    extended 75 bp from each side
    only exon region.

    this function check if the peak resided on exons that has splicings
    '''
    try:
        fetched = list(junctions.fetch(chrom, start, end))
        if any(_ for _ in fetched if _.split('\t')[5] == strand and int(_.split('\t')[4]) >= 3):
            return 'Spliced_exon:' + fetched[0].split('\t')[-1]
        else:
            return ''
    except ValueError:
        return ''


def retype_junctions(df):
    junctions_tab = '/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/merged_bam'\
                    '/unfragmentd.spliced.tsv.gz'
    junctions = pysam.Tabixfile(junctions_tab)
    
    check_junction = partial(is_junction_exon, junctions)
    return df\
        .assign(is_spliced_exon = lambda d: list(map(check_junction, d.chrom, d.start, d.end, d.strand))) \
        .assign(gtype = lambda d: np.where(d.is_spliced_exon.str.contains("Spliced_exon"), 
                                            'Long RNA',
                                            d.gtype)) \
        .assign(gname = lambda d: np.where(d.is_spliced_exon.str.contains("Spliced_exon"), 
                                            d.is_spliced_exon.str.split(':', expand=True).iloc[:,1],
                                            d.gname)) \
        .drop('is_spliced_exon', axis=1)




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
            .agg({'gname':collapse_col})
    else:
        max_overlapped = max_overlapped.nsmallest(1, 'gtype_rank')
    
    required_columns = ['gname','gtype','strand', 'gstrand']
    df_dict = {col: max_overlapped[col] for col in required_columns}
    max_overlapped = pd.DataFrame(df_dict)
    return max_overlapped


def strand_df(df, strand = 'sense'):
    '''
    Select gene annotation for each peak in a sense-specific manner
    '''
    condition = 'is_sense == "sense"' if strand == 'sense' else 'is_sense == "antisense"'
    return  df\
        .query(condition) \
        .rename(columns={'gname': strand + '_gname',
                  'gtype': strand + '_gtype' })  \
        .drop(['is_sense','gstrand'],axis=1)


def resolve_annotation(inbed):
    '''
    select for greatest overlapped annotation
    '''
    df = dd.from_pandas(inbed, npartitions=16, sort=True)\
        .assign(is_sense = lambda d: np.where((d.strand == d.gstrand) | (d.gtype.str.contains(repeats_regex)), 
                                            'sense', 
                                            'antisense')) \
        .groupby(['chrom','start','end',
                    'peakname','score','is_sense', 
                    'fc','log10p',
                    'log10q','pileup'])\
        .apply(select_annotation, meta={'gname':'f8',
                                        'gtype':'f8',
                                        'strand':'f8',
                                        'gstrand':'f8'}) 

    df = df.compute(scheduler='processes',
                      get = dask.multiprocessing.get) \
        .reset_index()\
        .drop_duplicates() \
        .sort_values('log10q', ascending=False)  \
        .assign(gtype = lambda d: d.gtype.map(merge_type)) \
        .drop('level_10', axis=1)  \
        .pipe(retype_junctions)
    
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
                                            'misc_RNA',
                                            d.gtype)) \
        .assign(gtype = lambda d: np.where(d.gtype.str.contains('srpRNA'), 
                                             'misc_RNA',
                                             d.gtype)) \
        .assign(gtype = lambda d: np.where(d.gname.isin(["BC200",'7SK','7SL','VTRNA2-1']), 
                                            'misc_RNA', 
                                            d.gtype)) \
        .assign(gtype_rank = lambda d: d.gtype.map(rank_type))   \
        .fillna(0)
    
    return inbed


def make_table(all=None, base_name = 'unfragmented'):
    #test:
    # all, base_name = True,  'unfragmented'
    project_path = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/merged_bed'
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
    print('Merging: ', ', '.join(map(os.path.basename, broad_peaks)))

    bed = pd.concat([process_broad(broad_peak, bed_path) for broad_peak in broad_peaks], sort=False) \
        .sort_values([0,1,2]) \
        .reindex() 
    

    inbed = annotate_peaks(annotation_file, bed)  

    df = resolve_annotation(inbed)  
#        .assign(seq = lambda d: list(map(fetch_seq, d.chrom, d.start, d.end, d.strand)))\
#        .assign(chrM = lambda d: d.seq.map(is_mt))
    df.to_csv(out_table, sep='\t', index=False)
    print('Written %s' %out_table)
    assert bed.shape[0] == df.shape[0], 'Peak lost!!!'


def main():
    all_annotations=[True, False]
    base_names = ['unfragmented', 'exonuclease'] 
    for all in all_annotations:
        for base_name in base_names:
            make_table(all=all, base_name = base_name) 

if __name__ == '__main__':
    main()
