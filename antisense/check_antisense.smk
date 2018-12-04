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
ANTISENSE_BED = SAMPLE_FOLDER_TEMPLATE + '/smallRNA/aligned.antisense.bed'
ANTISENSE_READ = SAMPLE_FOLDER_TEMPLATE + '/smallRNA/antisense_read.txt'
ANTISENSE_FQ = SAMPLE_FOLDER_TEMPLATE + '/smallRNA/antisense_read.fq.gz'
R2_ADAPTER_CONTAM_ANTISENSE_FQ = ANTISENSE_FQ.replace('.fq.gz','_contam.txt')
ANTISENSE_ANNOTATED_BED = SAMPLE_FOLDER_TEMPLATE + '/smallRNA/r2_annotated_antisense.bed'
NT_CUTOFF = 3


rule all:
    input:
        SUMMARY_TABLE,
        expand(ANTISENSE_ANNOTATED_BED, SAMPLENAME = SAMPLENAMES)

rule tag_antisense:
    input:
        BED = ANTISENSE_BED,
        TXT = R2_ADAPTER_CONTAM_ANTISENSE_FQ,

    output:
        BED = ANTISENSE_ANNOTATED_BED
    
    run:
        contam_info = {}
        with open(input.TXT) as contam:
            for line in contam:
                line = line.strip()
                if line.startswith('[') and line.endswith(')'):
                    seq_id = line.split(' ')[0].strip('[]')
                    umi = seq_id.split(':')[-1]
                    seq_id = ':'.join(seq_id.split(':')[:-1])
                    matched_nucleotide = line.split('(')[-1].split(',')[-2]
                    contam_info[umi + '_' + seq_id] = matched_nucleotide
        
        with open(input.BED) as bed, open(output.BED, 'w') as out:
            for line in bed:
                line = line.strip()
                seq_id = line.split('\t')[3]
                r2_num = contam_info.setdefault(seq_id, '0')
                print(line + '\t' + r2_num, file = out)








rule count_R2_contam:
    input:
        FQ = expand(ANTISENSE_FQ, SAMPLENAME = SAMPLENAMES)
    
    output:
        TXT = expand(R2_ADAPTER_CONTAM_ANTISENSE_FQ, SAMPLENAME = SAMPLENAMES),
        TABLE = SUMMARY_TABLE
    
    run: 
        rows = []
        for FQ, TXT, SAMPLENAME in  zip(input.FQ, output.TXT, SAMPLENAMES):
            seq_count, contam = find_r2(FQ, TXT, nt_cutoff=NT_CUTOFF)
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
        TXT = ANTISENSE_READ,
        BED = ANTISENSE_BED

    shell:
        'cat {input.BED} '\
        "| awk '$6==\"-\"'" \
        '| tee {output.BED} '\
        '| cut -f4 ' \
        "| sed 's/^[ACTG]\\{{6\\}}_//g' "\
        "> {output.TXT}"



