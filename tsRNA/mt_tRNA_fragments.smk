import glob
import os
import re

PROJECT_PATH = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
SAMPLE_FOLDERS = glob.glob(PROJECT_PATH + '/*')
SAMPLENAMES = map(os.path.basename, SAMPLE_FOLDERS)
SAMPLENAMES = filter(lambda x: re.search('001$',x), SAMPLENAMES)
SAMPLENAMES = list(SAMPLENAMES)

GENE_REF = os.environ['REF'] + '/hg19_ref/genes/rRNA_mt.bed'
tRNA_REF = os.environ['REF'] + '/hg19_ref/genes/smallRNA'
MT_tRNA_REF = os.environ['REF'] + '/hg19_reg/genes/mt_tRNA.fa'
MT_tRNA_SCAN = os.environ['REF'] + '/hg19_reg/genes/mt_tRNA.tRNA_scan'
SAMPLE_FOLDERS = PROJECT_PATH + '/{SAMPLENAME}'
MT_FOLDER = SAMPLE_FOLDERS + '/rRNA_mt'
MT_BAM = MT_FOLDER + '/aligned.bam'
MT_tRNA_FQ1 = MT_FOLDER + '/mt_tRNA.1.fq'
MT_tRNA_FQ2 = MT_FOLDER + '/mt_tRNA.2.fq'
MT_tRNA_BAM = MT_FOLDER + '/mt_tRNA.bam'
MT_tRNA_BED = MT_FOLDER + '/mt_tRNA.bed.gz'
MT_tRNA_BED_INDEX = MT_tRNA_BED + '.tbi'
THREADS = 1

rule all:
    input:
        expand(MT_tRNA_BED_INDEX, SAMPLENAME = SAMPLENAMES)

rule tRNA_bed_index:
    input:
        MT_tRNA_BED
    
    output:
        MT_tRNA_BED_INDEX
    
    shell:
        'tabix -p bed {input}'

rule tRNA_bed:
    input:
        MT_tRNA_BAM
    
    output:
        MT_tRNA_BED
    
    shell:
        'cat {input} '\
        '| bam_to_bed.py -i - --cigar '\
        '| sort -k1,1 -k2,2n -k3,3n '\
        '| bgzip '\
        '> {output}'

rule map_tRNA:
    input:
        FQ1 = MT_tRNA_FQ1,
        FQ2 = MT_tRNA_FQ2

    threads: THREADS
    params:
        INDEX = tRNA_REF

    output:
        MT_tRNA_BAM

    shell:
        'bowtie2 --very-sensitive-local '\
        '-L 8  --mp 4,2 -N 1 --mm '\
        '--no-mixed --no-discordant --dovetail '\
        '-p {threads} -x {params.INDEX} '\
        '-1 {input.FQ1} -2 {input.FQ2} '\
        '| samtools view -b@ {threads} '\
        '> {output}'


rule extract_fq:
    input:
        MT_BAM
    
    params:
        REF_MODEL = GENE_REF

    output:
        FQ1 = MT_tRNA_FQ1,
        FQ2 = MT_tRNA_FQ2
    
    shell:
        'cat {params.REF_MODEL} '\
        '| grep tRNA --color=no '\
        '| bedtools pairtobed -abam {input} -b - -type both '\
        '| samtools fastq - -1 {output.FQ1} -2 {output.FQ2}'

rule extract_mt_tRNA:
    input:
        tRNA_REF + '.fa'
    
    output:
        MT_tRNA_REF

    shell:
        'cat {input} '\
        '|  seqkit grep -p "MT" -n -r ' \
        '| seqtk seq '\
        '> {output} '


rule find_anticodon:
    input:
        MT_tRNA_REF

    output:
        MT_tRNA_SCAN
        

    shell:
        'tRNAscan-SE -M mammal -o {output} {input}'



