import os
import sys
import pandas as pd
sys.path.append('/stor/home/cdw2854/novel_rna/plots')
from utils import TrnaLookAlike


wildcard_constraints:
    SAMPLENAME = '[A-Z0-9a-z]+',
    COORDINATE = 'chr[0-9XY]+_[0-9]+_[0-9]+'

PROJECT_PATH = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
OUT_PATH = PROJECT_PATH + '/bed_files/bed_files/merged_bed/MACS2/consensus'
BAM = PROJECT_PATH + '/merged_bam/unfragmented.bam'
FQ1 = OUT_PATH +'/fastq/' + '/{COORDINATE}.1.fq'
FQ2 = OUT_PATH +'/fastq/' + '/{COORDINATE}.2.fq'
MERGED = OUT_PATH + '/fastq/' + '/{COORDINATE}.assembled.fastq'
SHAPES = OUT_PATH + '/shapes/' + '/{COORDINATE}.shapes'
COOR_CONSENSUS = OUT_PATH + '/fasta/{COORDINATE}.fa'
FA = OUT_PATH + '/fasta/'+ '/unfragmented.fa'
PEAK_FILE = PROJECT_PATH + '/bed_files/merged_bed/MACS2/annotated/additional_table.csv'
peaks = pd.read_csv(PEAK_FILE)


def get_coordinates(wildcard):
    td = peaks\
        .assign(tn = lambda d: d.chrom + '_' +d.start.astype(str)+ '_' + d.end.astype(str))
    return td.tn.tolist()


def coor_to_coor(wildcard):
    chrom, start, end = wildcard.COORDINATE.split('_')
    return chrom + ':' + start + '-' + end



rule all:
    input:
        FA


rule make_consensus:
    input:
        FA = lambda w: expand(COOR_CONSENSUS,
                    COORDINATE = get_coordinates(w))

    output:
        FA = FA

    shell:
        'cat {input} > {output}'


rule consensus:
    input:
        MERGED

    params:
        NAME = '{COORDINATE}' 

    output:
        FA = COOR_CONSENSUS

    shell:
        'cat {input} '\
        '| seqkit seq -m 30  '\
        '| seqtk seq -a | muscle  | seqtk seq '\
        '| python ~/ngs_qc_plot/muscle_to_concensus.py '\
        "| sed  's/^/\>{params.NAME}\\n/g'  "\
        '> {output} '


rule merge_fq:
    input:
        FQ1 = FQ1,
        FQ2 = FQ2

    params:
        NAME = MERGED.replace('.assembled.fastq','')

    output:
        MERGED

    shell:
        'pear -f {input.FQ1} -r {input.FQ2} -o {params.NAME}'


rule make_read_fq:
    input:
        BAM 
    
    params:
        COOR = lambda w: coor_to_coor(w)

    output:
        FQ1 = FQ1,
        FQ2 = FQ2

    shell:
        "samtools view -b {input} '{params.COOR}' "\
        '| samtools collate -O -  '\
        '| samtools fixmate - - '\
        '| samtools fastq -f0x2 - ' \
        '| deinterleave_fastq.py -i - -1 {output.FQ1} -2 {output.FQ2}'
