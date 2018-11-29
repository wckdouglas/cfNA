import glob
import os
import re

wildcard_constraints:
    DEDUP_PARAM="total|deduplicated",
    TREATMENT = "[a-zA-Z]+",
    RNA_TYPE = "[a-zRNA_]+"


#VARIABLES
PROJECT_PATH= '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
SAMPLE_FOLDERS = glob.glob(PROJECT_PATH + '/*001')
DEDUP_PARAM = ['total'] #'deduplicated'
OUT_BAM_TEMPLATE = PROJECT_PATH + "/merged_bam/small_rna/{TREATMENT}.{RNA_TYPE}.{DEDUP_PARAM}.bam"
SORT_BAM_TEMPLATE = OUT_BAM_TEMPLATE.replace('.bam','.nameSorted.bam')
SUBSAMPLE_BAM_TEMPLATE = OUT_BAM_TEMPLATE.replace('.bam','.subsampled.bam')
DEDUP_BAM = '{SAMPLE_FOLDER}/{RNA_TYPE}/aligned.sorted.deduplicated.bam'
TOTAL_BAM = '{SAMPLE_FOLDER}/{RNA_TYPE}/aligned.bam'
DEDUP_RG_BAM = DEDUP_BAM.replace('.bam','.add_rg.bam')
TOTAL_RG_BAM = TOTAL_BAM.replace('.bam','.add_rg.bam')
MIRNA_BAM_TEMPLATE = OUT_BAM_TEMPLATE.replace('.bam','.miRNA.bam')
VYRNA_BAM_TEMPLATE = OUT_BAM_TEMPLATE.replace('.bam','.vt_yRNA.bam')
REMAP_BAM_TEMPLATE = OUT_BAM_TEMPLATE.replace('.nameSorted.bam','.bowtie2.bam')
FQ1_TEMPLATE = OUT_BAM_TEMPLATE.replace('.bam','_R1.fq.gz')
FQ2_TEMPLATE = FQ1_TEMPLATE.replace('_R1.fq.gz','_R2.fq.gz')
THREADS=12
RNA_TYPES = ['rRNA_mt','smallRNA']
TREATMENTS = ['unfragmented','phosphatase','fragmented']
REGEXES = ['[Qq][cC][fF][0-9]+','[pP]hos[0-9]+','[fF]rag[0-9]+']

## wildcard functions
REGEX_DICT = {TREATMENT:REGEX for TREATMENT, REGEX in zip(TREATMENTS, REGEXES)}
def select_sample_folders(wildcards):
    REGEX = REGEX_DICT[wildcards.TREATMENT]
    return list(filter(lambda x: re.search(REGEX, x), SAMPLE_FOLDERS))

def select_bam(wildcards):
    bamlist = ''
    if wildcards.DEDUP_PARAM=="total":
        bamlist = expand(TOTAL_RG_BAM, 
                        SAMPLE_FOLDER=select_sample_folders(wildcards), 
                        RNA_TYPE = wildcards.RNA_TYPE) 
    else: 
        bamlist = expand(DEDUP_RG_BAM, 
                        SAMPLE_FOLDER=select_sample_folders(wildcards), 
                        RNA_TYPE = wildcards.RNA_TYPE)
    return bamlist
    
#SHARED commands
RG_COMMAND = 'samtools addreplacerg -r ID:{params.ID} -r SM:{params.ID} '\
        ' -o - {input.BAM} '\
        '| samtools view -b@ {params.THREADS} -F4 '\
        '| samtools sort -@ {params.THREADS} -O bam -o {output.BAM} '

rule all:
    input:
        expand(OUT_BAM_TEMPLATE, DEDUP_PARAM = DEDUP_PARAM, RNA_TYPE = RNA_TYPES, TREATMENT = TREATMENTS),
        expand(SUBSAMPLE_BAM_TEMPLATE, DEDUP_PARAM = DEDUP_PARAM, RNA_TYPE = RNA_TYPES, TREATMENT = TREATMENTS),
        expand(MIRNA_BAM_TEMPLATE, DEDUP_PARAM = DEDUP_PARAM, RNA_TYPE = ['smallRNA'], TREATMENT = TREATMENTS),
        expand(VYRNA_BAM_TEMPLATE, DEDUP_PARAM = DEDUP_PARAM, RNA_TYPE = ['smallRNA'], TREATMENT = TREATMENTS)


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
        THREADS = THREADS,
        FRACTION = lambda w: 1.01 if w.RNA_TYPE == 'smallRNA' else 1.1
    
    shell:
        'samtools view -bs {params.FRACTION} {input.BAM_FILE}' \
        '| filter_soft_clip.py -i - --pe '\
        '| python ~/ngs_qc_plot/bam_viz.py '\
        '| samtools view -b@ {params.THREADS}' \
        '| sambamba sort -t {params.THREADS} --show-progress '\
        '-o {output.BAM_FILE} /dev/stdin'


rule cat_bam:
    input:
        BAMS = lambda w: select_bam(w)

    output:
        BAM = SORT_BAM_TEMPLATE

    params:
        THREADS = THREADS

    shell:
        'sambamba merge -t {params.THREADS} /dev/stdout {input.BAMS} '\
        '| samtools view -b@ {params.THREADS} -F4 '\
        '| sambamba sort -n -t {params.THREADS} -o {output.BAM} /dev/stdin'


rule vt_yRNA_bam:
    input:
        SORT_BAM_TEMPLATE
    params:
        THREADS = THREADS,
        REF = os.environ['REF'] + '/hg19/new_genes/smallRNA.bed'
    output:
        VYRNA_BAM_TEMPLATE
    shell:
        'cat {params.REF} | egrep "vault|RNY" --color=no '\
        '| bedtools intersect -abam {input} -b - '\
        '| python ~/ngs_qc_plot/bam_viz.py '\
        '| samtools view -b@ {params.THREADS}' \
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
        '| python ~/ngs_qc_plot/bam_viz.py '\
        '| samtools view -b@ {params.THREADS}' \
        '| sambamba sort -t {params.THREADS} -o {output} /dev/stdin'


rule add_rg_total:
    input:
        BAM =  TOTAL_BAM 

    params:
        THREADS = THREADS,
        ID = lambda w: os.path.basename(w.SAMPLE_FOLDER)

    output:
        BAM = TOTAL_RG_BAM 

    shell:
        RG_COMMAND

rule add_rg_dedup:
    input:
        BAM =  DEDUP_BAM

    params:
        THREADS = THREADS,
        ID = lambda w: os.path.basename(w.SAMPLE_FOLDER)

    output:
        BAM = DEDUP_RG_BAM

    shell:
        RG_COMMAND
