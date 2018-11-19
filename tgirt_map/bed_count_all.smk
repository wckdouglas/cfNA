
import glob
import os
import sys
import re
from insert_size import bed_fragments, sample_fragments, add_dict
from collections import Counter
from functools import reduce
import pandas as pd

PROJECT_PATH='/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
COUNT_PATH = PROJECT_PATH + '/Counts/all_counts'
REF_BED_PATH = os.environ['REF'] + '/hg19/new_genes'
SAMPLE_FOLDERS = glob.glob(PROJECT_PATH + '/*001')
SAMPLE_FOLDERS = filter(lambda x: 'try' not in x, SAMPLE_FOLDERS)
SAMPLE_NAMES = list(map(os.path.basename, SAMPLE_FOLDERS))
RNA_TYPES = ['counts','sncRNA','small_RNA','rRNA_mt','reapeats']
DEDUP_TYPES = ['dedup','all']
STRANDS = ['sense', 'antisense']
COUNT_TEMPLATE = COUNT_PATH + '/{RNA_TYPE}/{SAMPLE_NAME}.{DEDUP}.{STRAND}.counts' 

SAMPLE_FOLDER_TEMPLATE = PROJECT_PATH + '/{SAMPLE_NAME}'
SAMPLE_COUNT_PATH = SAMPLE_FOLDER_TEMPLATE + '/count_temp'
INTERSECTED_TEMPLATE = SAMPLE_COUNT_PATH + '/{RNA_TYPE}.{DEDUP}.intersected.bed.gz'
BED_TEMPLATE = SAMPLE_COUNT_PATH + '/{RNA_TYPE}.bed.gz'
DEDUP_BED_TEMPLATE = SAMPLE_COUNT_PATH + '/{RNA_TYPE}.dedup.bed.gz'
PRIMARY_BAM = SAMPLE_FOLDER_TEMPLATE + '/Combined/primary_no_sncRNA_tRNA_rRNA.bam'
INSERT_SIZE_TABLE = PROJECT_PATH + '/fragment_sizes/{TREATMENT}.tsv'

TREATMENT_REGEX = ['Q[Cc][Ff][0-9]+|[ED][DE]|Exo|HS', 'Frag','[pP]hos', 
                  'L[1234]','All','N[aA][0-9]+',
                  'ED|DE','HS[123]','genome']
TREATMENTS = ['unfragmented','fragmented','phosphatase',
                'polyA','untreated', 'alkaline_hydrolysis',
                'exonuclease','high_salt','genome-sim'] 
treatment_regex_dict = {t:tr for t, tr in zip(TREATMENTS, TREATMENT_REGEX)}
def select_sample(wildcards, return_count = False):
    regex = treatment_regex_dict[wildcards.TREATMENT] 
    selected_samples = filter(lambda x: re.search(regex, x), SAMPLE_NAMES)
    selected_samples = list(selected_samples)
    if return_count:
        return len(selected_samples)
    else:
        return selected_samples

def select_treatment_bed(wildcards):
    selected_samples = select_sample(wildcards, return_count=False)
    selected_folders = expand(SAMPLE_FOLDER_TEMPLATE, SAMPLE_NAME = selected_samples)
    BEDS = []
    for selected_folder in selected_folders:
        BED = sample_fragments(selected_folder, return_beds = True)
        BEDS.extend(BED)
    return BEDS


rule all:
    input:
        expand(COUNT_TEMPLATE,
                SAMPLE_NAME = SAMPLE_NAMES,
                RNA_TYPE = RNA_TYPES,
                DEDUP = DEDUP_TYPES,
                STRAND = STRANDS),
        expand(INSERT_SIZE_TABLE,
                TREATMENT = TREATMENTS),


rule collect_insert_size:
    input:
        BEDS = lambda w:  select_treatment_bed(w)
    
    output:
        TABLE_NAME = INSERT_SIZE_TABLE
    
    run:
        size_dict = Counter()
        for bed in input.BEDS:
            size_dict += bed_fragments(bed)
        pd.DataFrame({'isize': list(size_dict.keys()),
                  'size_count': list(size_dict.values())})\
            .to_csv(output.TABLE_NAME, index=False, sep='\t')

rule count_bed:
    input:
        INTERSECTED_BED = INTERSECTED_TEMPLATE
    
    params:
        STRAND_FILTER = lambda w: strand_selection(w),
        TEMP = COUNT_TEMPLATE + '.TEMP',
        FIELDS = lambda w: field_selection(w)

    output:
        TABLE = COUNT_TEMPLATE
    
    shell:
        'mkdir -p {params.TEMP} '\
        '; zcat {input.INTERSECTED_BED} '\
        '| {params.STRAND_FILTER} '\
        '| {params.FIELDS} ' \
        '| sort --temporary-directory={params.TEMP} '\
        '| uniq -c '\
        "| awk {{'print $2,$3,$4,$5,$6,$7,$8,$9,$1'}} OFS='\\t' "\
        '> {output.TABLE}'
        '; rm -rf {params.TEMP}'


rule intersect_bed:
    input:
        BED = lambda w: BED_TEMPLATE if w.DEDUP=='all' else DEDUP_BED_TEMPLATE
    
    params:
        REF_BED = lambda w: get_parameter(w, return_format = 'ref')

    output:
        INTERSECTED_BED = INTERSECTED_TEMPLATE
    
    shell:
        'bedtools intersect -a {input.BED} -b {params.REF_BED} -wao '\
        '| bgzip > {output.INTERSECTED_BED}'


rule dedup_bed:
    input:
        BED = BED_TEMPLATE

    params:
        DEDUP = lambda w: deduplicate(w)

    output:
        BED = DEDUP_BED_TEMPLATE
    
    shell:
        'zcat {input.BED} '\
        '| {params.DEDUP} '\
        '| bgzip '\
        '> {output.BED} '



rule bam_to_bed:
    input:
        BAM = lambda w: get_parameter(w, return_format = 'bam')
    
    params:
        TEMP = BED_TEMPLATE + '.TEMP' 

    output:
        BED = BED_TEMPLATE

    shell:
        'mkdir -p {params.TEMP} '\
        '; cat {input.BAM} '\
        '| samtools view -bF 4 -F256 -F2048 '\
        '| bam_to_bed.py -i - --add_cigar --primary '\
        '| sort -k1,1 -k2,2n -k3,3n -k6,6 --temporary-directory={params.TEMP} '\
        '| bgzip '\
        '> {output.BED}'\
        ';rm -rf {params.TEMP} '\
    

rule filter_repeats:
    input:
        BAM = PRIMARY_BAM
    
    params:
        RMSK = os.environ['REF'] + '/hg19/genome/rmsk.bed.gz'

    output:
        BAM = PRIMARY_BAM.replace('.bam','_repeats.bam')
    
    shell:
        'bedtools pairtobed -type neither -abam {input.BAM} -b {params.RMSK}' \
        '> {output.BAM}'



def deduplicate(wildcards):
    '''
    genearte dedup command
    '''
    if not re.search('genome-sim|L[12E]|PEV', wildcards.SAMPLE_NAME):
        toleration = 0 
        DEDUP_COMMAND = "  deduplicate_bed.py -i - -d '_' -f 0 -t {TOLERATE} --ct 6" \
                        "| poisson_umi_adjustment.py -i - -o - --umi 6 "\
                        .format(TOLERATE=toleration)
    else:
        filter_command = "awk '$NF!~/N/ && ($3-$2) <= 300' |" if wildcards.RNA_TYPE == 'snc_all' else ''
        DEDUP_COMMAND = filter_command + ' sort -k1,1 -k2,2n -k3,3n -k6,6n -u ' \
                        '| cut -f1-6 '
    return DEDUP_COMMAND

def get_parameter(wildcards, return_format='bam'):

    SAMPLE_FOLDER = SAMPLE_FOLDER_TEMPLATE.format(SAMPLE_NAME = wildcards.SAMPLE_NAME) 
    if wildcards.RNA_TYPE == 'counts':
        REF_BED = REF_BED_PATH + '/genes.bed'
        BAM = SAMPLE_FOLDER  + '/Combined/primary_no_sncRNA_tRNA_rRNA_repeats.bam'
    elif wildcards.RNA_TYPE == 'sncRNA':
        REF_BED = REF_BED_PATH + '/sncRNA_x_protein.bed'
        BAM = SAMPLE_FOLDER + '/Combined/sncRNA.bam'
    elif wildcards.RNA_TYPE == 'small_RNA':
        REF_BED = REF_BED_PATH + '/smallRNA.bed'
        BAM = SAMPLE_FOLDER + '/smallRNA/aligned.bam'
    elif wildcards.RNA_TYPE == 'rRNA_mt':
        REF_BED = REF_BED_PATH + '/rRNA_mt.bed'
        BAM = SAMPLE_FOLDER + '/rRNA_mt/aligned.bam'
    elif wildcards.RNA_TYPE == 'repeats' or wildcards.RNA_TYPE == 'reapeats':
        REF_BED = os.environ['REF'] + '/hg19/genome/rmsk.bed.gz'
        BAM = SAMPLE_FOLDER + '/Combined/repeats.bam'
    
    if return_format == 'bam':
        return BAM

    elif return_format == 'ref':
        return REF_BED
    
    
def strand_selection(wildcards):
    if wildcards.DEDUP == "dedup":
        REF_STRAND = 12
    elif wildcards.DEDUP == 'all' :
        REF_STRAND = 13
    
    if wildcards.STRAND == 'sense':
        operator = '=='
    elif wildcards.STRAND == 'antisense':
        operator = '!='
    return "awk '$6 {operator} ${REF_STRAND} || ${REF_STRAND} == \".\" '".format(REF_STRAND = REF_STRAND, operator = operator)

def field_selection(wildcards):
    return ' cut -f7- ' if wildcards.DEDUP == 'dedup' else ' cut -f8-'
