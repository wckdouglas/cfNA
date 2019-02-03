#!/usr/bin/env python

import pandas as pd
import glob
import dask.dataframe as dd
from multiprocessing import Pool
import dask
import sys
import re
import os
THREADS = 24

def label_sample(x, salt = False):
    if 'HS' in x:
        return 'High salt (450mM)'
    elif 'Frag' in x:
        return 'Fragmented'
    elif re.search('[-_]sim',x):
        return 'WGS-sim'
    elif re.search('N[aA]|[Aa]lk', x):
        #return 'Alkaline hydrolysis'
        return 'NaOH'
    elif re.search('_L[0-9]+',x):
        return 'Poly(A)-selected'
    elif re.search('[eE]xo|ED|DE', x):
        return 'DNase I + Exo I'
    elif re.search('[aA]ll|[Uu]nt', x):
        return 'Untreated'
    elif re.search('Phos', x):
        return 'DNase I + Phosphatase'
    elif re.search('[Qq][cC][Ff][0-9]+|[uU]nf', x):
        if salt:
            return 'Low salt (200mM)'
        else:
            return 'DNase I'

def read_ref():
    ref_tab = '/stor/work/Lambowitz/ref/hg19_ref/tRNA_fragments/tRNA_fragments.tsv'
    return pd.read_table(ref_tab, skiprows= 6)\
        .pipe(lambda d: d[~d.iloc[:,0].str.startswith('<')]) \
        .rename(columns = {'#MINTbase Unique ID (sequence derived)':'sseqid'})\
        .filter(regex = "sseq|^[53]'-t|^i") \
        .pipe(pd.melt, id_vars = 'sseqid', 
              value_name = 'tRNA', var_name = 'tRF')\
        .pipe(lambda d: d[~pd.isnull(d.tRNA)])


def read_sample(SAMPLE_FOLDER):
    '''
    read blast output
    '''

    col_names = ['qseqid', 'qlen', 'sseqid',
            'slen', 'pident', 'length',
            'mismatch', 'gapopen', 'qstart', 
            'qend', 'sstart', 'send', 'evalue']
    tRF_tab = SAMPLE_FOLDER + '/blast.tRF.tsv'
    samplename = os.path.basename(SAMPLE_FOLDER)
    print('Reading %s' %samplename)
    tRF_df = dd.read_table(tRF_tab, names = col_names, )\
        .repartition(npartitions=THREADS)  \
        .query('slen == qlen ') \
        .groupby('qseqid')\
        .apply(lambda d: d.nlargest(1, 'pident'))\
        .compute(workers=THREADS, scheduler='threads')\
        .reset_index(drop=True) \
        .assign(samplename = samplename)
    return tRF_df


def make_table(project_path, make = False):
    out_table = project_path + '/tRF.feather'
    SAMPLE_FOLDERS = glob.glob(project_path + '/Q*')

    if make:
        p = Pool(THREADS)
        dfs = p.map(read_sample, SAMPLE_FOLDERS)
        p.close()
        p.join()

        pd.concat(dfs, sort=True)\
            .reset_index(drop=True) \
            .to_feather(out_table)
        print('Written %s' %out_table)
    return out_table


def clean_table(table_name):
    return pd.read_feather(table_name) \
        .assign(umi = lambda d: d.qseqid.str.extract('^([ACTG]{6})_', expand=False)) \
        .assign(read_count = 1)\
        .filter(['sseqid','umi','samplename'])\
        .drop_duplicates()\
        .groupby(['sseqid','samplename'], as_index=False)\
        .agg({'umi':'count'}) \
        .merge(read_ref(), on = 'sseqid')\
        .assign(prep = lambda d: d.samplename.map(label_sample))


            
def main():
    project_path = '/stor/work/Lambowitz/cdw2854/cfNA/tsRNA'
    out_table = make_table(project_path, make = False)
    clean_table(out_table).to_feather(out_table)

        
if __name__ == '__main__':
    main()