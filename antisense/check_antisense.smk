import glob
import os
import sys
from check_r1 import find_r2
from pandas import DataFrame


PROJECT_PATH = '/stor/work/Lambowitz/cdw2854/cfNA'
DATA_PATH = PROJECT_PATH + '/data'
FOLDER_PATH = PROJECT_PATH + '/tgirt_map'
SAMPLE_FOLDERS = glob.glob(FOLDER_PATH + '/Q*001')
SAMPLENAME_REGEX = '[Q][cC][fF].*001$'
SAMPLENAMES = map(os.path.basename, SAMPLE_FOLDERS)
SAMPLENAMES = filter(lambda x: re.search(SAMPLENAME_REGEX,x ), SAMPLENAMES)
SAMPLENAMES = list(SAMPLENAMES)

#SNAKEMAKE VARIABLE
SUMMARY_TABLE = 'recopy.csv'
RAW_FQ = DATA_PATH + '/{SAMPLENAME}.fastq.gz'
SAMPLE_FOLDER_TEMPLATE = FOLDER_PATH + '/{SAMPLENAME}'
SMALL_RNA_BED = SAMPLE_FOLDER_TEMPLATE + '/smallRNA/aligned.bed'
ANTISENSE_READ = SAMPLE_FOLDER_TEMPLATE + '/smallRNA/antisense_read.txt'
ANTISENSE_FQ = SAMPLE_FOLDER_TEMPLATE + '/smallRNA/antisense_read.fq.gz'
R2_ADAPTER_CONTAM_ANTISENSE_FQ = ANTISENSE_FQ.replace('.fq.gz','_contam.txt')


rule all:
    input:
        SUMMARY_TABLE

rule count_R2_contam:
    input:
        FQ = expand(ANTISENSE_FQ, SAMPLENAME = SAMPLENAMES)
    
    output:
        TXT = expand(R2_ADAPTER_CONTAM_ANTISENSE_FQ, SAMPLENAME = SAMPLENAMES),
        TABLE = SUMMARY_TABLE
    
    run: 
        rows = []
        for FQ, TXT, SAMPLENAME in  zip(input.FQ, output.TXT, SAMPLENAMES):
            seq_count, contam = find_r2(FQ, TXT, nt_cutoff=6)
            rows.append((seq_count, contam, SAMPLENAME))
        DataFrame(rows, columns = ['anti_seq_count','with_R2', 'samplename'])\
            .to_csv(output.TABLE, index=False)




rule extract_fq:
    input:
        READ_IDS = ANTISENSE_READ
    
    params:
        FQ = RAW_FQ,

    output:
        FQ = ANTISENSE_FQ
    
    shell:
        'cat  {params.FQ} '\
        '| seqkit grep -f {input.READ_IDS} -o {output.FQ}' 


rule find_antisense:
    input:
        BED = SMALL_RNA_BED

    output:
        TXT = ANTISENSE_READ

    shell:
        'cat {input.BED} '\
        "| awk '$6==\"-\"'" \
        '| cut -f4 ' \
        "| sed 's/^[ACTG]\\{{6\\}}_//g' "\
        "> {output.TXT}"



