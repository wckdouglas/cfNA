#!/usr/bin/env python

from multiprocessing import Pool
import os
from operator import itemgetter
from collections import defaultdict, Counter
import glob
import pandas as pd
import pysam
from sequencing_tools.bam_tools.poisson_umi_tools import correct_umi_count
from itertools import groupby
import re

def coordinate_grouper(x):
    fields = x.strip().split('\t')
    return itemgetter(0,1,2,5,6)(fields)

class tsRNA_counter():
    def __init__(self, bed, samplename, anticodon_table):
        self.bed = pysam.Tabixfile(bed)
        self.aligned_contig = set(self.bed.contigs)
        self.anticodon_df = pd.read_table(anticodon_table)\
                .pipe(lambda d: d[d.tRNA.isin(self.aligned_contig)])
        self.tsRNA_counter = defaultdict(Counter)
        self.samplename = samplename
        self.splice = re.compile('([0-9]+)S')
        self.aligned = re.compile('([0-9]+)[SM]')
    
    def count_tRNA(self):
        for _, tRNA in self.anticodon_df.iterrows():
            iterator = self.bed.fetch(tRNA['tRNA'])
            for (chrom, start, end, strand, cigar), lines in groupby(iterator, coordinate_grouper):
                clipped = sum(map(int,self.splice.findall(cigar)))
                total_base = sum(map(int, self.aligned.findall(cigar)))

                if clipped/total_base < 0.2 and strand == '+':
                    start, end = int(start), int(end)
                    umis = set()
                    for line in lines:
                        fields = line.split('\t')
                        umi = fields[3].split('_')[0]
                        umis.add(umi)
                    read_count = correct_umi_count(len(umis), umi_nt = 6)
                    
                    start_in_anticodon = tRNA['anticodon_start'] -1 < start < tRNA['anticodon_end']+1
                    end_at3 = end > tRNA['end'] - 5 
                    start_at5 = start <  5 
                    end_in_anticodon = tRNA['anticodon_start'] - 1 < end < tRNA['anticodon_end'] + 1
                    short_frag = end - start < 23

                    if end_at3 and short_frag:
                        self.tsRNA_counter[tRNA['tRNA']]["3' tsRNA"] += read_count

                    elif end_at3 and start_in_anticodon:
                        self.tsRNA_counter[tRNA['tRNA']]["3' half"] += read_count
                    
                    elif start_at5 and end_in_anticodon:
                        self.tsRNA_counter[tRNA['tRNA']]["5' half"] += read_count
                    
                    else:
                        self.tsRNA_counter[tRNA['tRNA']]['Others'] += read_count

    
    def write_table(self):
        dfs = []
        for tRNA, frag_dict in self.tsRNA_counter.items():
            df = pd.DataFrame({'frag_type':list(frag_dict.keys()),
                               'frag_count':list(frag_dict.values())})\
                .assign(samplename =self. samplename)\
                .assign(tRNA = tRNA)
            dfs.append(df)

        if dfs:
            return pd.concat(dfs)\
                .reset_index(drop=True) \
                .merge(self.anticodon_df\
                        .filter(['tRNA','anticodon','aa']),
                    on = 'tRNA')
        else:
            return None


def count_ts(sample_folder):
    samplename = os.path.basename(sample_folder)
    tablename = sample_folder + '/count_temp/tRNA_frag.feather'

    # genomic tRNA
    anticodon_table = '/stor/work/Lambowitz/ref/hg19_ref/genes/tRNA/anticodon_annotations.tsv'
    bed = sample_folder + '/count_temp/small_RNA.all.bed.gz'
    os.system('tabix -f -p bed %s' %bed)

    # mt tRNA
    mt_anticodon_table = '/stor/work/Lambowitz/ref/hg19_ref/genes/mt_tRNA.anticodon_annotations.tsv'
    mt_bed = sample_folder + '/rRNA_mt/mt_tRNA.bed.gz'
    os.system('tabix -f -p bed %s' %mt_bed)

    if True: #not os.path.isfile(tablename):
        print('Running %s' %sample_folder)
        tsRNA = tsRNA_counter(mt_bed, samplename, mt_anticodon_table)
        tsRNA.count_tRNA()
        mt_tab = tsRNA.write_table()

        tsRNA = tsRNA_counter(bed, samplename, anticodon_table)
        tsRNA.count_tRNA()
        tab = tsRNA.write_table()

        if tab is None and mt_tab is None:
            tablename = None
        elif tab is None:
            ts_table = mt_tab
        elif mt_tab is None:
            ts_table = tab
        else:
            ts_table = pd.concat([mt_tab, tab]).reset_index(drop=True)
        ts_table.to_feather(tablename)
        print('Written %s' %tablename)
    return tablename


def main():
    project_path = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
    sample_folders = glob.glob(project_path + '/*001')
    sample_folders = filter(lambda x: not re.search('genome-sim|L[0-9E]+',x), sample_folders)
    sample_folders = filter(lambda x: re.search('[qQ][cC][fF]', x), sample_folders)
    
    p = Pool(24)
    dfs = p.map(count_ts, sample_folders)
    p.close()
    p.join()

    tablename = project_path + '/Counts/tsRNA.feather'
    dfs = filter(lambda df: df is not None, dfs)
    pd.concat(map(pd.read_feather, dfs))\
        .reset_index(drop=True)\
        .to_feather(tablename)
    print('Written %s' %tablename)

if __name__ == '__main__':
    main()



                



