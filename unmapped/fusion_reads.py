#!/usr/bin/env python

import os
import glob
import pysam
import pandas as pd
from multiprocessing import Pool
import sys


def read_flag_stat(filename):
    info = open(filename,'r').readlines()
    stat_df = {}
    stat_df['read1'] = get_number(get_line(info, 'read1'))
    stat_df['mapped'] = get_number(get_line(info, 'mapped'))
    stat_df['supplementary'] = get_number(get_line(info, 'supplementary'))
    stat_df['proper pair'] = get_number(get_line(info, 'properly paired'))
    stat_df['secondary'] = get_number(get_line(info, 'secondary'))
    return stat_df

def bam_fusion(bam_file):
    fusion_read = set()
    with pysam.Samfile(bam_file,'r') as bam:
        for aln in bam:
            if not aln.is_unmapped and (not aln.is_proper_pair or aln.is_supplementary):
                fusion_read.add(aln.query_name)
    return len(fusion_read)




def read_sample(sample_folder):
    samplename = os.path.basename(sample_folder)
    bam_file = sample_folder + '/unmapped/rescued.bam'
    fusion_count = bam_fusion(bam_file)
    return (samplename, fusion_count)


project_path = '/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map'
sample_folders = glob.glob(project_path + '/*001')

p = Pool(24)
counts = p.map(read_sample, sample_folders)
p.close()
p.join()




pd.DataFrame(counts, columns = ['samplename','fusion_reads']) \
    .sort_values('samplename')\
    .to_csv(sys.stdout, index=False, sep='\t')
