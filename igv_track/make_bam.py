#!/usr/bin/env python

import os
import glob
import re
from picard import run_rna_seq_picard

refflat = '/stor/work/Lambowitz/ref/hg19/genome/refFlat.txt'
snc_annotation = os.environ['REF'] + '/hg19/genome/sncRNA_rRNA_for_bam_filter.bed'
protein_bed = os.environ['REF'] + '/hg19/genome/protein.bed'
plus_bed = os.environ['REF'] + '/hg19/genome/protein_plus.bed'
minus_bed = os.environ['REF'] + '/hg19/genome/protein_minus.bed'

project_path = '/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/'
folders = glob.glob(project_path + '/*001')
folders.sort()


for regex, label in zip(['Q[Cc][Ff][0-9]+|[ED][DE]|Exo|HS', 'Frag', 'L[12]','All','N[aA]','ED|DE','HS[123]'],
                        ['unfragmented','fragmented','polyA','untreated', 'alkaline_hydrolysis','exonuclease','high_salt']):
    samples = filter(lambda x: re.search(regex, x), folders)
    if 'poly' in label:
        bam = 'primary.bam'
    else:
        bam = 'primary.marked_duplicate.bam'
    bam_files = list(map(lambda x: x + '/Combined/%s' %bam, samples))


    merge_path = project_path + '/merged_bam'
    filter_path = merge_path + '/filtered_bam'
    merged_bam = merge_path + '/' + label + '.bam' 
    name_sorted_bam = merge_path + '/' + label + '.name_sort.bam' 
    filtered_bam = filter_path + '/' + label + '.bam'

    plus_bam = filter_path + '/' + label + '.plus.bam' 
    plus_anti_bam = filter_path + '/' + label + '.plus_anti.bam' 
    plus_sense_bam = filter_path + '/' + label + '.plus_sense.bam' 
    minus_bam = filter_path + '/' + label + '.minus.bam'
    minus_anti_bam = filter_path + '/' + label + '.minus_anti.bam'
    minus_sense_bam = filter_path + '/' + label + '.minus_sense.bam'

    protein_sense_bam = filter_path + '/' + label + '.protein_sense.bam'
    protein_anti_bam = filter_path + '/' + label + '.protein_anti.bam'


    if len(bam_files) > 1:
        command = 'sambamba merge -p -t 6 /dev/stdout  {bamfiles}'.format(bamfiles = ' '.join(bam_files))
    else:
        command = 'cat {bam_files}'.format(bam_files = bam_files[0])
    command += ' | sambamba sort -p -t 6 -o {out_bam} /dev/stdin ' \
                .format(label = label, 
                        out_bam = merged_bam)

    filtering_plus = ' awk \'{if (($1~/^@/)|| ($2==99)||($2==147)|| ($2==355)||($2==403)|| ($2==1123)||($2==1171)) print}\''
    filtering_minus = 'awk \'{if (($1~/^@/)|| ($2==83)||($2==163)|| ($2== 339)||($2==419)|| ($2==1187)||($2==1107)) print}\''

    command += '; sambamba sort -n -p -o {name_sorted} {merged_bam}'\
            '; samtools view -h {name_sorted} | {filtering_plus} ' \
            '| samtools view -bS - > {plus_bam} '\
            '; samtools view -h {name_sorted} | {filtering_minus} '\
            '| samtools view -bS - > {minus_bam} '\
            '; bedtools pairtobed -abam {plus_bam} -b {sncRNA} -type neither | bedtools pairtobed -a - -b {plus_bed} > {plus_sense_bam} '\
            '; bedtools pairtobed -abam {minus_bam} -b {sncRNA} -type neither | bedtools pairtobed -a - -b {minus_bed} > {minus_sense_bam} '\
            '; bedtools pairtobed -type neither -abam {minus_bam} -b {minus_bed} '\
            '| bedtools pairtobed -abam - -b {plus_bed} > {minus_anti_bam} '\
            '; bedtools pairtobed -type neither -abam {plus_bam} -b {plus_bed} '\
            '| bedtools pairtobed -abam - -b {minus_bed} > {plus_anti_bam} '\
            '; sambamba merge /dev/stdout {plus_sense_bam} {minus_sense_bam} | samtools view -bF 1024 | sambamba sort {protein_sense_bam} /dev/stdin'\
            '; sambamba merge /dev/stdout  {plus_anti_bam} {minus_anti_bam} | samtools view -bF 1024 | sambamba sort {protein_anti_bam} /dev/stdin '\
            .format(merged_bam = merged_bam,
                    name_sorted = name_sorted_bam,
                    filtering_plus = filtering_plus,
                    filtering_minus = filtering_minus,
                    plus_bam = plus_bam,
                    minus_bam = minus_bam,
                    sncRNA = snc_annotation, 
                    plus_bed = plus_bed,
                    minus_bed = minus_bed,
                    plus_sense_bam = plus_sense_bam,
                    minus_sense_bam = minus_sense_bam,
                    minus_anti_bam = minus_anti_bam,
                    plus_anti_bam = plus_anti_bam,
                    protein_sense_bam = protein_sense_bam,
                    protein_anti_bam = protein_anti_bam)

#    command += '; bedtools intersect -abam {bam} -b {protein} '\
#            ' | bedtools intersect -abam - -b {sncRNA} -v -split '\
#            '| samtools view -bF 1024 -F 256 -F 2048 ' \
#            ' > {filtered_bam} '\
#            .format(bam = merged_bam,
#                    sncRNA = snc_annotation,
#                    protein = protein_bed,
#                    filtered_bam = filtered_bam)
    command += '; '
    command += run_rna_seq_picard(refflat, filter_path, protein_sense_bam, label, return_command=True)
    command += '; '
    command += run_rna_seq_picard(refflat, filter_path, protein_anti_bam, label, return_command=True)
    print(command)




