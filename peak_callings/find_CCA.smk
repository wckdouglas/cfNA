#!/usr/bin/env python

import os

PROJECT_PATH = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/merged_bam/dedup'
OUT_PATH = PROJECT_PATH + '/fold' 
BAM = PROJECT_PATH + '/unfragmented.chrM_filter.dedup.bam'
FQ1 = OUT_PATH + '/{NAME}.1.fq'
FQ2 = OUT_PATH + '/{NAME}.2.fq'
MERGED = OUT_PATH + '/{NAME}.assembled.fastq'
SHAPES = OUT_PATH + '/{NAME}.shapes'
FA = OUT_PATH + '/{NAME}.fa'
ALIGNED_FA = OUT_PATH + '/unknown_tRNA.multi.fa'
coors = ['chr22:24,349,641-24,349,711', 'chr10:71,355,032-71,355,103',
        'chr1:181,392,082-181,392,121', 'chr10:101817588-101817660',
        'chr9:5095156-5095227',
        'chr4:156,379,949-156,380,020',
        'chr2:202,077,805-202,077,873',
        'chr16:20733601-20733674',
        'chr2:125438494-125438563',
        'chr9:81357660-81357728']
names = ['chr22_unknown', 'chr10_unknown',
        'CACNA1E', 'CPN1',
        'JAK2',
        'chr4_unknown',
        'CASP10',
        'THUMPD1',
        'CNTNAP5',
        'chr9_unknown']
names = [c if n.endswith('unknown') else n for c,n in zip(coors, names)]
coors = {n:c for c, n, in zip(coors, names)}


rule all:
    input:
        ALIGNED_FA

rule multi_aligned:
    input:
        expand(FA, NAME = names)

    output:
        ALIGNED_FA

    shell:
        'cat {input} | muscle > {output}'

rule shape:
    input:
        FA = FA

    output:
        SHAPE = SHAPES

    shell:
        'cat {input.FA} '\
        "| RNAshapes -s -e 10 -l -m \'[[][][]]\' "\
        '> {output.SHAPE}'

rule concensus:
    input:
        MERGED

    params:
        NAME = '{NAME}'
    output:
        FA = FA

    shell:
        'cat {input} '\
        '| seqkit seq -m 30 -M 80 '\
        '| seqtk seq -a | muscle  | seqtk seq '\
        '| python muscle_to_concensus.py '\
        "| sed  's/^/\>{params.NAME}\\n/g'  "\
        '> {output} '


rule merge_fq:
    input:
        FQ1 = FQ1,
        FQ2 = FQ2

    params:
        NAME = OUT_PATH + '/{NAME}'

    output:
        MERGED

    shell:
        'pear -f {input.FQ1} -r {input.FQ2} -o {params.NAME}'


rule make_read_fq:
    input:
        BAM 
    
    params:
        COOR = lambda w: coors[w.NAME]

    output:
        FQ1 = FQ1,
        FQ2 = FQ2

    shell:
        "samtools view -b {input} '{params.COOR}' "\
        '| samtools sort -n '\
        '| samtools fastq - ' \
        '| deinterleave_fastq.py -i - -1 {output.FQ1} -2 {output.FQ2}'
