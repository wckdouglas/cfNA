#!/usr/bin/env python

import os
import glob
import re
import pandas as pd

def filter_bed(run_filter=True):
    REF = os.environ['REF']
    if run_filter:
        filter_command = "| awk '$3-$2 > 10 && $3-$2 < 200' " +\
                '| bedtools intersect -a - -b {REF}/hg19/new_genes/tRNA_yRNA.bed -v '  \
                '| bedtools intersect -a - -b {REF}/hg19/genome/tRNA.bed -v ' \
                '| bedtools intersect -a - -b {REF}/hg19/new_genes/rmsk_tRNA.bed -v ' \
                '| bedtools intersect -a - -b {REF}/hg19/new_genes/rmsk_rRNA.bed -v '\
                '| bedtools intersect -a - -b {REF}/hg19/new_genes/rmsk_yRNA.bed -v '\
                '| bedtools intersect -a - -b {REF}/hg19/new_genes/refseq_rRNA.bed -v '\
                '| bedtools intersect -a - -b {REF}/hg19/new_genes/rRNA_for_bam_filter.bed -v '\
                .format(REF=REF)
    else:
        filter_command =' '
    return filter_command


def output_merge_bed(bed_files, outname):
    '''
    sample_df containing columns:
    bed: full path to bed file

    '''
    if re.search('alkaline', outname):
        filter_command = filter_bed(run_filter=False)
    else:
        filter_command = filter_bed(run_filter=True)
    command = ' zcat {in_beds} '\
            ' {filter_command}' \
            '| sort -k1,1 -k2,2n -k3,3n '\
            '| bgzip '\
            '> {out_bed}'\
            '; tabix -p bed -f {out_bed}'\
            .format(in_beds = ' '.join(bed_files),
                    filter_command = filter_command,
                    out_bed = outname)
    print (command)



project_path = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map_30'
out_path = project_path + '/CLAM/BED_files'
sample_folder = glob.glob(project_path + '/*001') 
sample_folder.sort()

if not os.path.isdir(out_path):
    os.makedirs(out_path)

for regex, label in zip(['Q[Cc][Ff][0-9]+|Exo|[DE][DE]', 'N[aA][0-9]+', 'All'],
                        ['unfragmented','alkaline','untreated']):
    samples = filter(lambda x: re.search(regex, os.path.basename(x).split('no')[0]), sample_folder)
    bed_files = map(lambda x: x + '/Combined/CLAM/pseudo_fragment.bed.gz', samples)
    bed_files = list(bed_files)

    outname = out_path + '/{label}.bed.gz'.format(label = label)
    output_merge_bed(bed_files, outname)

