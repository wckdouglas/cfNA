import glob
import os
import re

wildcard_constraints:
    DEDUP_PARAM="total|deduplicated"


PROJECT_PATH= '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
SAMPLE_FOLDER = glob.glob(PROJECT_PATH + '/*001')
SAMPLE_FOLDER = list(filter(lambda x: re.search('[Qq][cC][fF][0-9]+',x), SAMPLE_FOLDER))
DEDUP_PARAM = ['total'] #'deduplicated'
OUT_BAM_TEMPLATE = PROJECT_PATH + "/merged_bam/small_rna/smallRNA.{DEDUP_PARAM}.bam"
SORT_BAM_TEMPLATE = OUT_BAM_TEMPLATE.replace('.bam','.nameSorted.bam')
SUBSAMPLE_BAM_TEMPLATE = OUT_BAM_TEMPLATE.replace('.bam','.subsampled.bam')
MIRNA_BAM_TEMPLATE = OUT_BAM_TEMPLATE.replace('.bam','.miRNA.bam')
VYRNA_BAM_TEMPLATE = OUT_BAM_TEMPLATE.replace('.bam','.vt_yRNA.bam')
REMAP_BAM_TEMPLATE = OUT_BAM_TEMPLATE.replace('.nameSorted.bam','.bowtie2.bam')
FQ1_TEMPLATE = OUT_BAM_TEMPLATE.replace('.bam','_R1.fq.gz')
FQ2_TEMPLATE = FQ1_TEMPLATE.replace('_R1.fq.gz','_R2.fq.gz')
THREADS=12


def select_bam(wildcards):
    if wildcards.DEDUP_PARAM == 'total':
        bam_file = "smallRNA/aligned.bam"
    else:
        bam_file = "smallRNA/aligned.sorted.deduplicated.bam"
    return list(map(lambda S: S + '/' + bam_file, SAMPLE_FOLDER))


rule all:
    input:
        expand(OUT_BAM_TEMPLATE, DEDUP_PARAM = DEDUP_PARAM),
        expand(SUBSAMPLE_BAM_TEMPLATE, DEDUP_PARAM = DEDUP_PARAM),
        expand(MIRNA_BAM_TEMPLATE, DEDUP_PARAM = DEDUP_PARAM),
        expand(VYRNA_BAM_TEMPLATE, DEDUP_PARAM = DEDUP_PARAM)
        #REMAP_BAM_TEMPLATE.format(DEDUP_PARAM = 'total'),


rule sort_bam:
    input:
        BAM_FILE = SORT_BAM_TEMPLATE

    output:
        BAM_FILE = OUT_BAM_TEMPLATE

    params:
        THREADS = THREADS

    shell:
        'cat {input.BAM_FILE}' \
        '| python ~/ngs_qc_plot/bam_viz.py '\
        '| samtools view -b@ {params.THREADS}' \
        '| sambamba sort -t {params.THREADS} --show-progress '\
        '-o {output.BAM_FILE} /dev/stdin'

rule subsample_bam:
    input:
        BAM_FILE = SORT_BAM_TEMPLATE
    
    output:
        BAM_FILE = SUBSAMPLE_BAM_TEMPLATE
    
    params:
        THREADS = THREADS
    
    shell:
        'samtools view -bs 1.005 {input.BAM_FILE}' \
        '| filter_soft_clip.py -i - --pe '\
        '| python ~/ngs_qc_plot/bam_viz.py '\
        '| samtools view -b@ {params.THREADS}' \
        '| sambamba sort -t {params.THREADS} --show-progress '\
        '-o {output.BAM_FILE} /dev/stdin'


rule cat_bam:
    input:
        BAMS = select_bam

    output:
        BAM = SORT_BAM_TEMPLATE

    params:
        THREADS = THREADS

    shell:
        'samtools cat {input.BAMS} '\
        '| sambamba sort -n -o {output.BAM} -t {params.THREADS} /dev/stdin '\


rule vt_yRNA_bam:
    input:
        SORT_BAM_TEMPLATE
    params:
        THREADS = THREADS,
        REF = os.environ['REF'] + '/hg19/new_genes/smallRNA.bed'
    output:
        VYRNA_BAM_TEMPLATE
    shell:
        'cat {params.REF} | grep "vault|RNY" --color=no '\
        '| bedtools intersect -abam {input} -b - '\
        '| sambamba sort -t {params.THREADS} -o {output} /dev/stdin'
        
        

        

rule miRNA_bam:
    input:
        SORT_BAM_TEMPLATE
    params:
        THREADS = THREADS,
        miRNA = os.environ['REF'] + '/hg19/new_genes/smallRNA.bed'
    output:
        MIRNA_BAM_TEMPLATE
    
    shell:
        'cat {params.miRNA} | grep hsa --color=no '\
        '| bedtools intersect -abam {input} -b - '\
        '| sambamba sort -t {params.THREADS} -o {output} /dev/stdin'
