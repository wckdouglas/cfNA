import glob
import os
import re

wildcard_constraints:
    DEDUP_PARAM="total|deduplicated",
    TREATMENT = "[a-zA-Z0-9_\-]+",
    RNA_TYPE = "[a-zRNA_]+",
    STRAND = 'fwd|rvs',
    FRAG_TYPE = 'five|three'


#VARIABLES
PROJECT_PATH= '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
SAMPLE_FOLDERS = glob.glob(PROJECT_PATH + '/*001')
SMALL_RNA_BED = os.environ['REF'] + '/hg19_ref/genes/smallRNA.bed'
RNA_FAI = os.environ['REF'] + '/hg19_ref/genes/{RNA_TYPE}.viz.genome'
DEDUP_PARAM = ['total'] #'deduplicated'
OUT_PATH = PROJECT_PATH + "/merged_bam/small_rna"
OUT_BAM_TEMPLATE =  OUT_PATH +"/{TREATMENT}.{RNA_TYPE}.{DEDUP_PARAM}.bam"
REV_OUT_BAM_TEMPLATE =  OUT_PATH +"/{TREATMENT}.{RNA_TYPE}.{DEDUP_PARAM}.reverse.bam"
SORT_BAM_TEMPLATE = OUT_BAM_TEMPLATE.replace('.bam','.nameSorted.bam')
SUBSAMPLE_BAM_TEMPLATE = OUT_BAM_TEMPLATE.replace('.bam','.subsampled.bam')
DEDUP_BAM = '{SAMPLE_FOLDER}/{RNA_TYPE}/aligned.sorted.deduplicated.bam'
TOTAL_BAM = '{SAMPLE_FOLDER}/{RNA_TYPE}/aligned.bam'
DEDUP_RG_BAM = DEDUP_BAM.replace('.bam','.add_rg.bam')
TOTAL_RG_BAM = TOTAL_BAM.replace('.bam','.add_rg.bam')
MIRNA_BAM_TEMPLATE = OUT_BAM_TEMPLATE.replace('.bam','.miRNA.bam')
VYRNA_BAM_TEMPLATE = OUT_BAM_TEMPLATE.replace('.bam','.vt_yRNA.bam')
tRNA_FRAG_BAM_TEMPLATE = OUT_BAM_TEMPLATE.replace('.bam','.tRNA_frag.bam')
FILTERED_tRNA_FRAG_BAM_TEMPLATE = OUT_BAM_TEMPLATE.replace('.bam','.{FRAG_TYPE}_tRNA_frag.bam')
BED_TEMPLATE = OUT_BAM_TEMPLATE.replace('.bam','.bed.gz')
BG_TEMPLATE = OUT_BAM_TEMPLATE.replace('.bam','.{STRAND}.bedGraph')
BIGWIG_TEMPLATE = OUT_BAM_TEMPLATE.replace('.bam','.{STRAND}.bigWig')
REMAP_BAM_TEMPLATE = OUT_BAM_TEMPLATE.replace('.nameSorted.bam','.bowtie2.bam')
FQ1_TEMPLATE = OUT_BAM_TEMPLATE.replace('.bam','_R1.fq.gz')
FQ2_TEMPLATE = FQ1_TEMPLATE.replace('_R1.fq.gz','_R2.fq.gz')
THREADS=12
RNA_TYPES = ['rRNA_mt','smallRNA']
TREATMENTS = ['unfragmented','phosphatase','fragmented','HEK293','alkaline_hydrolysis']
STRANDS = ['fwd','rvs']
REGEXES = ['[Qq][cC][fF][0-9]+','[pP]hos[0-9]+','[fF]rag[0-9]+','GC','_[Nn][aA][0-9]+_']

## wildcard functions
REGEX_DICT = {TREATMENT:REGEX for TREATMENT, REGEX in zip(TREATMENTS, REGEXES)}
def select_sample_folders(wildcards):
    REGEX = REGEX_DICT[wildcards.TREATMENT]
    return list(filter(lambda x: re.search(REGEX, x), SAMPLE_FOLDERS))


def cat_bam_command(wildcards):
    if len(select_bam(wildcards)) == 1:
        command = 'samtools sort -@ %i -o ' %THREADS
    else:
        command = 'sambamba merge -t %i ' %THREADS
    return command


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
        expand(REV_OUT_BAM_TEMPLATE, DEDUP_PARAM = DEDUP_PARAM, RNA_TYPE = RNA_TYPES, TREATMENT = TREATMENTS),
        expand(SUBSAMPLE_BAM_TEMPLATE, DEDUP_PARAM = DEDUP_PARAM, RNA_TYPE = RNA_TYPES, TREATMENT = TREATMENTS),
        expand(MIRNA_BAM_TEMPLATE, DEDUP_PARAM = DEDUP_PARAM, RNA_TYPE = ['smallRNA'], TREATMENT = TREATMENTS),
        expand(VYRNA_BAM_TEMPLATE, DEDUP_PARAM = DEDUP_PARAM, RNA_TYPE = ['smallRNA'], TREATMENT = TREATMENTS, STRAND = STRANDS),
        expand(FILTERED_tRNA_FRAG_BAM_TEMPLATE, DEDUP_PARAM = DEDUP_PARAM, 
                RNA_TYPE = ['smallRNA'], FRAG_TYPE = ['five','three'], TREATMENT = TREATMENTS, STRAND = STRANDS),
        expand(tRNA_FRAG_BAM_TEMPLATE, DEDUP_PARAM = DEDUP_PARAM, 
                RNA_TYPE = ['smallRNA'], TREATMENT = TREATMENTS, STRAND = STRANDS),
        expand(BIGWIG_TEMPLATE, DEDUP_PARAM = DEDUP_PARAM, RNA_TYPE = RNA_TYPES, TREATMENT = TREATMENTS, STRAND = STRANDS),


rule bg_to_bw:
    input:
        BG = BG_TEMPLATE

    params:
        GENOME = RNA_FAI

    output:
        BW = BIGWIG_TEMPLATE

    shell:
        'bedGraphToBigWig {input.BG} {params.GENOME} {output.BW}'


rule small_rna_coverage:
    input:
        BED = BED_TEMPLATE

    params:
        GENOME = RNA_FAI,
        STRAND = lambda w: '"+"' if w.STRAND == "fwd" else '"-"'

    output:
        BG = BG_TEMPLATE

    shell:
        'zcat {input.BED} '\
        "| awk '$6=={params.STRAND}' "\
        "| awk {{'print $1, $2+50, $3+50'}} OFS='\\t' "\
        '| bedtools genomecov -i - -g {params.GENOME} -bga '\
        '| bedtools sort -i - '\
        '> {output.BG}'


rule small_rna_bed:
    input:
        BAM = SORT_BAM_TEMPLATE

    params:
        TEMP_DIR = OUT_PATH

    output:
        BED = BED_TEMPLATE,
        TABIX = BED_TEMPLATE + '.tbi'

    shell:
        'bam_to_bed.py -i {input.BAM} -o - '\
        '| sort -k1,1 -k2,2n -k3,3n -T {params.TEMP_DIR}'\
        '| bgzip '\
        '> {output.BED}'\
        '; tabix -p bed -f {output.BED}'

rule reverse_bam:
    input:
        BAM_FILE = OUT_BAM_TEMPLATE

    params:
        TMP_DIR = REV_OUT_BAM_TEMPLATE + '_tmp'

    output:
        BAM_FILE = REV_OUT_BAM_TEMPLATE

    threads: THREADS
    shell:
        'mkdir -p {params.TMP_DIR} '\
        '; bamtools filter -in {input} -script reverse_filter.json '\
        '| sambamba sort -t {threads} --show-progress '\
        '-o {output.BAM_FILE} --tmpdir={params.TMP_DIR} /dev/stdin'\
        '; rm -rf {params.TMP_DIR}'


rule sort_bam:
    input:
        BAM_FILE = SORT_BAM_TEMPLATE

    output:
        BAM_FILE = OUT_BAM_TEMPLATE

    params:
        THREADS = THREADS,
        TMP_DIR = OUT_BAM_TEMPLATE + '_tmp'

    shell:
        'mkdir -p {params.TMP_DIR} '\
        '; cat {input.BAM_FILE}' \
        '| python ~/ngs_qc_plot/bam_viz.py '\
        '| samtools view -b@ {params.THREADS}' \
        '| sambamba sort -t {params.THREADS} --show-progress '\
        '-o {output.BAM_FILE} --tmpdir={params.TMP_DIR} /dev/stdin'\
        '; rm -rf {params.TMP_DIR}'

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
        THREADS = THREADS,
        CAT_COMMAND = lambda w: cat_bam_command(w),
        TEMP_DIR = SORT_BAM_TEMPLATE + '_TMP'

    shell:
        'mkdir -p {params.TEMP_DIR} '
        '; {params.CAT_COMMAND} /dev/stdout {input.BAMS} '\
        '| samtools view -b@ {params.THREADS} -F4 '\
        '| sambamba sort -n -t {params.THREADS} --tmpdir {params.TEMP_DIR} -o {output.BAM} /dev/stdin'\
        '; rm -rf {params.TEMP_DIR} '



rule tsFrag:
    input:
        SORT_BAM_TEMPLATE

    params:
        FRAG_TYPE = lambda w: w.FRAG_TYPE,
        THREADS = THREADS,
        TEMP_DIR = FILTERED_tRNA_FRAG_BAM_TEMPLATE + '_tmp'

    output:
        FILTERED_tRNA_FRAG_BAM_TEMPLATE

    shell:
        'mkdir -p {params.TEMP_DIR} '\
        '; cat {input} '\
        '| filter_soft_clip.py -i - --pe '\
        '| python filter_end.py -i - -o - --type {params.FRAG_TYPE} '\
        '| samtools view -h '\
        '| python ~/ngs_qc_plot/bam_viz.py '\
        '| samtools view -b '\
        '| sambamba sort -t {params.THREADS} -o {output} --tmpdir {params.TEMP_DIR} /dev/stdin '\
        '; rm -rf {params.TEMP_DIR} '


rule tRNA_fragments:
    input:
        SORT_BAM_TEMPLATE
    
    params:
        THREADS = THREADS,
        REF = SMALL_RNA_BED,
        TEMP_DIR = tRNA_FRAG_BAM_TEMPLATE + '_tmp'

    output:
        tRNA_FRAG_BAM_TEMPLATE

    shell:
        'mkdir -p {params.TEMP_DIR} '\
        '; cat {params.REF} | egrep "^TR" --color=no '\
        '| bedtools pairtobed -abam {input} -b - -type both '\
        '| filter_soft_clip.py --pe -i - -o - '\
        '| bamtools filter -script filter.json '\
        '| python ~/ngs_qc_plot/bam_viz.py '\
        '| samtools view -b@ {params.THREADS} ' \
        '| sambamba sort -t {params.THREADS} --tmpdir {params.TEMP_DIR}  -o {output} /dev/stdin'\
        '; rm -rf {params.TEMP_DIR} '

rule vt_yRNA_bam:
    input:
        SORT_BAM_TEMPLATE
    params:
        THREADS = THREADS,
        REF = SMALL_RNA_BED,
        TEMP_DIR = VYRNA_BAM_TEMPLATE + '_tmp'

    output:
        VYRNA_BAM_TEMPLATE
    shell:
        'mkdir -p {params.TEMP_DIR} '\
        '; cat {params.REF} | egrep "vault|RN[Y7]" --color=no '\
        '| bedtools intersect -abam {input} -b - '\
        '| python ~/ngs_qc_plot/bam_viz.py '\
        '| samtools view -b@ {params.THREADS}' \
        '| sambamba sort -t {params.THREADS} -o {output} --tmpdir {params.TEMP_DIR} /dev/stdin'
        '; rm -rf {params.TEMP_DIR} '
        
        

rule miRNA_bam:
    input:
        SORT_BAM_TEMPLATE

    params:
        THREADS = THREADS,
        BED = SMALL_RNA_BED

    output:
        MIRNA_BAM_TEMPLATE
    
    shell:
        'cat {params.BED} | grep hsa --color=no '\
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
