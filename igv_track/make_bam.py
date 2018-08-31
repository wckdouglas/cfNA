#!/usr/bin/env python

import os
import glob
import re
from picard import run_rna_seq_picard, filter_protein_bam

refflat = '/stor/work/Lambowitz/ref/hg19/new_genes/proteins.refflat'
threads = 6

project_path = '/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/'
merge_path = project_path + '/merged_bam'
filter_path = merge_path + '/filtered_bam'
folders = glob.glob(project_path + '/*001')
folders.sort()

def make_dir(p):
    if not os.path.isdir(p):
        os.mkdir(p)

make_dir(merge_path)
make_dir(filter_path)

for regex, label in zip(['Q[Cc][Ff][0-9]+|[ED][DE]|Exo|HS', 'Frag', 'L[12]','All','N[aA]','ED|DE','HS[123]'],
                        ['unfragmented','fragmented','polyA','untreated', 'alkaline_hydrolysis','exonuclease','high_salt']):
    samples = filter(lambda x: re.search(regex, x), folders)
    if 'poly' in label:
        bam = 'primary.bam'
    else:
        bam = 'primary.marked_duplicate.bam'
    bam_files = list(map(lambda x: x + '/Combined/%s' %bam, samples))


    merged_bam = merge_path + '/' + label + '.bam' 



    if len(bam_files) > 1:
        command = 'sambamba merge -p -t {threads} /dev/stdout  {bamfiles}'\
                    .format(threads = threads,
                            bamfiles = ' '.join(bam_files))
    else:
        command = 'cat {bam_files}'.format(bam_files = bam_files[0])
    command += ' | sambamba sort -p -t {threads} -o {out_bam} /dev/stdin ' \
                .format(threads = threads,
                        out_bam = merged_bam)

    command += filter_protein_bam(merged_bam, filter_path, label, refflat)
    print(command)




