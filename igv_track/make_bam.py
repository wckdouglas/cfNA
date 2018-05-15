#!/usr/bin/env python

import os
import glob
import re
from picard import run_rna_seq_picard

refflat = '/stor/work/Lambowitz/ref/hg19/genome/refFlat.txt'
snc_annotation = os.environ['REF'] + '/hg19/genome/sncRNA_x_protein.sorted.bed.gz'
protein_bed = os.environ['REF'] + '/hg19/genome/protein.bed'

project_path = '/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/'
folders = glob.glob(project_path + '/*001')
folders.sort()


for regex, label in zip(['Q[Cc][Ff][0-9]+|[ED][DE]|Exo|HS', 'Frag', 'L[12]','All','N[aA]'],
                        ['unfragmented','fragmented','polyA','untreated', 'alkaline_hydrolysis']):
    samples = filter(lambda x: re.search(regex, x), folders)
    if 'poly' in label:
        bam = 'primary.bam'
    else:
        bam = 'primary.marked_duplicate.bam'
    bam_files = list(map(lambda x: x + '/Combined/%s' %bam, samples))


    merge_path = project_path + '/merged_bam'
    filter_path = merge_path + '/filtered_bam'
    merged_bam = merge_path + '/' + label + '.bam' 
    filtered_bam = filter_path + '/' + label + '.bam'


    if len(bam_files) > 1:
        command = 'sambamba merge -p -t 6 /dev/stdout  {bamfiles}'.format(bamfiles = ' '.join(bam_files))
    else:
        command = 'cat {bam_files}'.format(bam_files = bam_files[0])
    command += ' | sambamba sort -p -t 6 -o {out_bam} /dev/stdin ' \
                .format(label = label, 
                        out_bam = merged_bam)

    command += '; bedtools intersect -abam {bam} -b {protein} '\
            ' | bedtools intersect -abam - -b {sncRNA} -v '\
            '| samtools view -bF 1024 -F 256 -F 2048 ' \
            ' > {filtered_bam} '\
            .format(bam = merged_bam,
                    sncRNA = snc_annotation,
                    protein = protein_bed,
                    filtered_bam = filtered_bam)
    command += '; '
    command += run_rna_seq_picard(refflat, filter_path, filtered_bam, label,
                                  return_command=True)

    print (command)




