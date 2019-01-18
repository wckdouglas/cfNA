import glob
import os
import sys
import re
from insert_size import bed_fragments, sample_fragments, add_dict
from collections import Counter, defaultdict
from functools import reduce
import pandas as pd

RNA_TYPES = ['counts','sncRNA','small_RNA','rRNA_mt','repeats']
DEDUP_TYPES = ['dedup','all']
STRANDS = ['sense', 'antisense']
wildcard_constraints:
    RNA_TYPE = '[A-Za-z_]+',
    DEDUP = 'dedup|all',
    STRAND = '[antisense]+',
    SAMPLE_NAME = '[A-Za-z0-9_\-]+',
    TREATMENT = '[A-Za-z_\-0-9]+'

PROJECT_PATH='/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
COUNT_PATH = PROJECT_PATH + '/Counts/all_counts'
REF_BED_PATH = os.environ['REF'] + '/hg19_ref/genes'
SAMPLE_FOLDERS = glob.glob(PROJECT_PATH + '/*001')
SAMPLE_FOLDERS = filter(lambda x: 'try' not in x, SAMPLE_FOLDERS)
SAMPLE_NAMES = list(map(os.path.basename, SAMPLE_FOLDERS))
COUNT_TEMPLATE = COUNT_PATH + '/{RNA_TYPE}/{SAMPLE_NAME}.{DEDUP}.{STRAND}.counts' 

SAMPLE_FOLDER_TEMPLATE = PROJECT_PATH + '/{SAMPLE_NAME}'
SAMPLE_COUNT_PATH = SAMPLE_FOLDER_TEMPLATE + '/count_temp'
INTERSECTED_TEMPLATE = SAMPLE_COUNT_PATH + '/intersected/{RNA_TYPE}.{DEDUP}.bed.gz'
BED_TEMPLATE = SAMPLE_COUNT_PATH + '/{RNA_TYPE}.bed.gz'
DEDUP_BED_TEMPLATE = SAMPLE_COUNT_PATH + '/{RNA_TYPE}.{DEDUP}.bed.gz'
PRIMARY_BAM = SAMPLE_FOLDER_TEMPLATE + '/Combined/primary_no_sncRNA_tRNA_rRNA.bam'
INSERT_SIZE_TABLE = PROJECT_PATH + '/fragment_sizes/{TREATMENT}.feather'
COUNT_TABLE = COUNT_PATH + '/spreaded_all_counts.feather'
LONG_COUNT_TABLE = COUNT_PATH + '/all_counts.feather'

TREATMENT_REGEX = ['Q[Cc][Ff][0-9]+|[ED][DE]|Exo|HS', 'Frag','[pP]hos', 
                  'L[1234]','All','N[aA][0-9]+',
                  'ED|DE','HS[123]','genome',
                    'MPF4','MPF10','MPCEV','^GC',
                    'PPF4','PPF10','PPCEV']
TREATMENTS = ['unfragmented','fragmented','phosphatase',
                'polyA','untreated', 'alkaline_hydrolysis',
                'exonuclease','high_salt','genome-sim',
                'EV','RNP','RNP-EV','HEK293',
                'MNase_EV','MNase_RNP','MNase_EV-RNP'] 


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
        LONG_COUNT_TABLE, COUNT_TABLE,
        expand(INSERT_SIZE_TABLE,
            TREATMENT = TREATMENTS),

rule make_table:
    input:
        expand(COUNT_TEMPLATE,
                SAMPLE_NAME = SAMPLE_NAMES,
                RNA_TYPE = RNA_TYPES,
                DEDUP = DEDUP_TYPES,
                STRAND = STRANDS),
        
    output:
        COUNT_TABLE,
        LONG_COUNT_TABLE
        
    shell:
        'python make_count_table.py'
        

rule collect_insert_size:
    input:
        BEDS = lambda w:  select_treatment_bed(w)
    
    output:
        TABLE_NAME = INSERT_SIZE_TABLE
    
    run:
        size_dict = defaultdict(Counter)
        for bed in input.BEDS:
            size_dict[bed] += bed_fragments(bed)

        dfs = []
        for bed, bed_dict in size_dict.items():
            df = pd.DataFrame({'isize': list(bed_dict.keys()),
                      'size_count': list(bed_dict.values())})\
                .assign(bed = bed)
            dfs.append(df)
        dfs = pd.concat(dfs, sort=False)\
            .reset_index(drop=True)
        dfs.to_feather(output.TABLE_NAME)

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
        '| python selective_count.py '\
        '| cut -f5,6,7,8,9,10,11,12'\
        '| sort --temporary-directory={params.TEMP} '\
        '| uniq -c '\
        "| awk {{'print $2,$3,$4,$5,$6,$7,$8,$9,$1'}} OFS='\\t' "\
        '> {output.TABLE}'
        '; rm -rf {params.TEMP}'

rule intersect_bed:
    input:
        BED = DEDUP_BED_TEMPLATE
    
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
    if wildcards.DEDUP == "dedup":
        if not re.search('genome-sim|L[12E]|PEV', wildcards.SAMPLE_NAME):
            toleration = 0 
            DEDUP_COMMAND = "  deduplicate_bed.py -i - -d '_' -f 0 -t {TOLERATE} --ct 6" \
                            "| poisson_umi_adjustment.py -i - -o - --umi 6 "\
                            .format(TOLERATE=toleration)
        else:
            filter_command = "awk '$NF!~/N/ && ($3-$2) <= 300' |" if wildcards.RNA_TYPE == RNA_TYPES[1] else ''
            DEDUP_COMMAND = filter_command + ' sort -k1,1 -k2,2n -k3,3n -k6,6n -u ' \
                            '| cut -f1-6 '
    else:
        DEDUP_COMMAND = ' cat '
    return DEDUP_COMMAND

def get_parameter(wildcards, return_format='bam'):


    SAMPLE_FOLDER = SAMPLE_FOLDER_TEMPLATE.format(SAMPLE_NAME = wildcards.SAMPLE_NAME) 

    if wildcards.RNA_TYPE == RNA_TYPES[0]:
        REF_BED = REF_BED_PATH + '/genes.bed'
        BAM = SAMPLE_FOLDER  + '/Combined/primary_no_sncRNA_tRNA_rRNA_repeats.bam'
    elif wildcards.RNA_TYPE == RNA_TYPES[1]:
        REF_BED = REF_BED_PATH + '/sncRNA_x_protein.bed'
        BAM = SAMPLE_FOLDER + '/Combined/sncRNA.bam'
    elif wildcards.RNA_TYPE == RNA_TYPES[2]:
        REF_BED = REF_BED_PATH + '/smallRNA.bed'
        BAM = SAMPLE_FOLDER + '/smallRNA/aligned.bam'
    elif wildcards.RNA_TYPE == RNA_TYPES[3]:
        REF_BED = REF_BED_PATH + '/rRNA_mt.bed'
        BAM = SAMPLE_FOLDER + '/rRNA_mt/aligned.bam'
    elif wildcards.RNA_TYPE == RNA_TYPES[4] or wildcards.RNA_TYPE:
        REF_BED = os.environ['REF'] + '/hg19/genome/rmsk.bed.gz'
        BAM = SAMPLE_FOLDER + '/Combined/repeats.bam'
#    else:
#        print(wildcards.SAMPLE_NAME, wildcards.RNA_TYPE)
#        BAM=''
#        REF_BED = ''
    
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
    return ' cut -f1,2,3,4,7- ' if wildcards.DEDUP == 'dedup' else ' cut -f1,2,3,4,8-'
