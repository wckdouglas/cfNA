
#!/usr/bin/env python

from __future__ import print_function
from multiprocessing import Pool
import glob
import os
import pandas as pd
import sys
import pysam
from pybedtools import BedTool
import pyximport
from operator import itemgetter
pyximport.install(setup_args={'include_dirs': pysam.get_include()})
import junction_function

def index_table(table_name):
    os.system('cat {tab} | sort -k1,1 -k2,2n -k3,3n '\
              '| bgzip > {tab}.gz'\
              '; tabix -f -p bed {tab}.gz'\
              .format(tab = table_name))
    print('Indexed %s' %table_name)


def add_gene(table_name):
    protein_bed = '/stor/work/Lambowitz/ref/hg19/new_genes/protein.bed.gz'
    BedTool(table_name) \
        .intersect(s=True, wb = True, b = protein_bed) \
        .to_dataframe(usecols = [0,1,2,3,4,5,9])\
        .to_csv(table_name, index=False, header=False, sep='\t')
    print('Intersected protein')


def make_table(out_table, bam_file):
    junction_table = junction_function.find_junction(bam_file) 
    with open(out_table, 'w') as out:
        for i, row in junction_table \
                .assign(chrom = lambda d: d.intron.str.split(':', expand=True).iloc[:,0]) \
                .assign(start = lambda d: d.intron.str\
                            .extract(':([0-9]+)-[0-9]+_[+-]', expand=False).astype(int) - 100) \
                .assign(end = lambda d: d.intron.str\
                            .extract(':[0-9]+-([0-9]+)_[+-]', expand=False).astype(int) + 100) \
                .assign(strand = lambda d: d.intron.str\
                            .extract(':[0-9]+-[0-9]+_([+-])', expand=False)) \
                .filter(['chrom','start','end','intron','intron_count','strand'])\
                .iterrows():

            line_template = '\t'.join([row['chrom'],'{start}',
                                    '{end}', row['intron'],
                                    str(row['intron_count']),
                                    row['strand']])

            start1 = row['start'] - 75
            end1 = row['start'] 

            start2 = row['end']
            end2 = row['end'] + 75


            if start1 > 0:
                print(line_template.format(start = start1, end = end1), file = out)
                print(line_template.format(start = start2, end = end2), file= out)
    print('Written %s' %out_table)


def run(bam_file, out_table):
    make_table(out_table, bam_file)
    add_gene(out_table)
    index_table(out_table)


def main():
    project_path = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/merged_bam'
    bam_file = project_path + '/unfragmented.bam'
    out_table = project_path + '/unfragmentd.spliced.tsv'
    run(bam_file, out_table)

if __name__ == '__main__':
    main()
   
