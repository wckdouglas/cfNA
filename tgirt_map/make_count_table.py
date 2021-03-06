#!/usr/bin/env python

import pandas as pd
import numpy as np
import glob
import os
import re
import sys
from tgirt_map.table_tools import change_gene_type
from multiprocessing import Pool


def read_count_file(file_count, samplename, count_file, count_type, strand, dedup):
    try:
        count_mat = pd.read_csv(count_file, usecols=[3,6,7,8],
                    names=['gene_name','gene_type','gene_id','read_count'],
                    engine='python', sep='\t')  \
            .groupby(['gene_name','gene_type','gene_id'], as_index=False)\
            .sum() \
            .assign(gene_type = lambda d: np.where(d.gene_type == ".", 'No features', d.gene_type))  

    except KeyError:
        print('Error:', count_file)
        sys.exit()
    
    if count_type in ["Repeats", 'reapeats', 'repeats']:
        count_mat = count_mat \
            .assign(gene_name = lambda d: d.gene_type + ':' + d.gene_name) \
            .assign(gene_type = 'Repeats')

    if count_type not in ['sncRNA', 'small_RNA']:
        count_mat = count_mat #\
#            .pipe(lambda d: d[~d.gene_type.isin(['miRNA','misc_RNA','snoRNA','snRNA','tRNA','rRNA','piRNA'])])


    if file_count % 20 == 0:
        print('Parsed %i files' %file_count)

    if count_type == "rRNA_mt":
        count_mat = count_mat.query('gene_type != "No features"')

    return count_mat \
            .query('read_count > 0') \
            .assign(samplename = samplename) \
            .assign(strand = strand)  \
            .assign(dedup = dedup)  

def read_function(args):
    file_count, samplename, count_file, count_type, strand, dedup= args
    print('Running %s' %count_file)
    return read_count_file(file_count, samplename, count_file, count_type, strand, dedup)

def main():
    count_path = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/Counts/all_counts'
    count_types = os.listdir(count_path)
    count_types = map(lambda x: count_path + '/' + x, count_types)
    count_files = []
    for count_type in count_types:
        count_files.extend(glob.glob(count_type + '/*'))

    sample_df = pd.DataFrame({'count_file': count_files}) \
        .assign(count_type = lambda d: list(map(lambda x: x.split('/')[-2], d.count_file))) \
        .assign(strand = lambda d: list(map(lambda x: x.split('.')[-2], d.count_file))) \
        .assign(dedup = lambda d: list(map(lambda x: x.split('.')[-3], d.count_file))) \
        .assign(samplename = lambda d: list(map(lambda x: x.split('/')[-1].split('.')[0], 
                                                d.count_file))) \
        .pipe(lambda d: d[d.samplename.str.contains('[Qq][cC][fF]|[PM]P|genome')])
    sample_df.to_csv('sample.tsv',sep='\t', index=False)
    print(sample_df.head())

    #debug_df = sample_df[sample_df.count_file.str.contains('Qcf12')].query('dedup=="dedup" & strand == "sense"')

    print ('Combining %i files' %sample_df.shape[0])
    iterable = []
    for i, row in sample_df.iterrows():
        iterable.append((i, row['samplename'], row['count_file'],
                        row['count_type'], row['strand'], row['dedup']))

    run_concat = True
    long_tablename = count_path + '/all_counts.tsv'
    spreaded_tablename = count_path + '/spreaded_all_counts.tsv'
    if run_concat:
        p = Pool(24)
        dfs = p.map(read_function, iterable)
        p.close()
        p.join()

        dfs = filter(lambda d: d is not None, dfs)
        concat_df = pd.concat(dfs, axis=0)  \
                .groupby(['samplename','strand','gene_type','gene_name','gene_id', 'dedup'], as_index=False)\
                .agg({'read_count':'sum'})
        concat_df.to_csv(long_tablename, sep = '\t', index=False)
        concat_df.to_feather(long_tablename.replace('.tsv','.feather'))
        print('Written %s' %(long_tablename))

    concat_df = pd.read_table(long_tablename)\
        .assign(samplename = lambda d: d.samplename + ':'+ d.dedup+':' + d.strand)\
        .assign(grouped_type = lambda d: d.gene_type.map(change_gene_type)) \
        .pipe(pd.pivot_table,
            index = ['gene_id','gene_name','gene_type', 'grouped_type'],
            columns = 'samplename', 
            aggfunc=np.sum,
            fill_value = 0,
            values = 'read_count') \
        .reset_index() 
    concat_df.to_csv(spreaded_tablename, sep = '\t', index=False)
    concat_df.to_feather(spreaded_tablename.replace('.tsv','.feather'))
    print('Written %s' %(spreaded_tablename))





if __name__ == '__main__':
    main()
