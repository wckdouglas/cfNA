
import glob
import os
import sys
from check_r1 import find_r2
from pandas import DataFrame

wildcard_constraints:
    READ_END = '[12]',
    SAMPLENAME = '[A-Za-z_0-9]+'

PROJECT_PATH = '/stor/work/Lambowitz/cdw2854/cfNA'
DATA_PATH = PROJECT_PATH + '/data'
FOLDER_PATH = PROJECT_PATH + '/tgirt_map'
SAMPLE_FOLDERS = glob.glob(FOLDER_PATH + '/Q*001')
SAMPLENAME_REGEX = '[Q][cC][fF].*001$'
SAMPLENAMES = map(os.path.basename, SAMPLE_FOLDERS)
SAMPLENAMES = filter(lambda x: re.search(SAMPLENAME_REGEX,x ), SAMPLENAMES)
SAMPLENAMES = list(SAMPLENAMES)
SMALL_RNA_INDEX = os.environ['REF'] + '/hg19_ref/genes/smallRNA'
TRANSCRIPTOME_INDEX = os.environ['REF'] + '/hg19_ref/genes/protein.fa'
GENOME_INDEX = os.environ['REF'] + '/hg19_ref/genome/hg19_genome.fa'

#SNAKEMAKE VARIABLE
SUMMARY_TABLE = 'recopy.csv'
SAMPLE_FOLDER_TEMPLATE = FOLDER_PATH + '/{SAMPLENAME}/smallRNA'
SMALL_RNA_BED = SAMPLE_FOLDER_TEMPLATE + '/aligned.bed'
ANTISENSE_BED = SAMPLE_FOLDER_TEMPLATE + '/aligned.antisense.bed'
ANTISENSE_READ = SAMPLE_FOLDER_TEMPLATE + '/antisense_read.txt'
ANTISENSE_FQ = SAMPLE_FOLDER_TEMPLATE + '/antisense_read.{READ_END}.fq.gz'
ANTISENSE_FILTER_CONTAM_FQ = SAMPLE_FOLDER_TEMPLATE + '/antisense_read.filtered.{READ_END}.fq.gz'
TRANSCRIPTOME_FILTERED_FQ = SAMPLE_FOLDER_TEMPLATE + '/antisense_read.filtered.transcriptome.{READ_END}.fq.gz'
GENOME_FILTERED_FQ = SAMPLE_FOLDER_TEMPLATE + '/antisense_read.filtered.genome.{READ_END}.fq.gz'
TRANSCRIPTOME_BAM = SAMPLE_FOLDER_TEMPLATE + '/transcriptome.bam'
GENOME_BAM = SAMPLE_FOLDER_TEMPLATE + '/genome.bam'
R2_ADAPTER_CONTAM_ANTISENSE_FQ = SAMPLE_FOLDER_TEMPLATE + '/antisense_read_R1_contam.txt'
ANTISENSE_ANNOTATED_BED = SAMPLE_FOLDER_TEMPLATE + '/r2_annotated_antisense.bed'
REVERSE_BAM = SAMPLE_FOLDER_TEMPLATE + '/reverse.bam'

COMBINED_BAM_PATH = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/merged_bam/small_rna'
COMBINED_BAM = COMBINED_BAM_PATH + '/{TREATMENT}.reverse.bam'
NT_CUTOFF = 3
THREADS = 3

TREATMENTS = ['unfragmented','phosphatase','fragmented', 'NaOH']
REGEXES = ['[Qq][cC][fF][0-9]+','[pP]hos[0-9]+','[fF]rag[0-9]+', '[Nn][aA][0-9]+']
REGEX_DICT = {TREATMENT:REGEX for TREATMENT, REGEX in zip(TREATMENTS, REGEXES)}
def select_sample(wildcards):
    REGEX = REGEX_DICT[wildcards.TREATMENT]
    return list(filter(lambda x: re.search(REGEX, x), SAMPLENAMES))


rule all:
    input:
        SUMMARY_TABLE,
        expand(ANTISENSE_ANNOTATED_BED, SAMPLENAME = SAMPLENAMES),
        expand(COMBINED_BAM, TREATMENT = TREATMENTS),


rule combine_bam:
    input:
        #BAMS = lambda w: expand(SORTED_REVERSE_BAM, SAMPLENAME = select_sample(w))
        BAMS = lambda w: expand(REVERSE_BAM, SAMPLENAME = select_sample(w))
    
    threads: THREADS
    params:
        TEMP_DIR = COMBINED_BAM_PATH

    output:
        BAM = COMBINED_BAM
    
    shell:
        'sambamba merge -t {threads} /dev/stdout {input.BAMS} ' \
        '| python ~/ngs_qc_plot/bam_viz.py '\
        '| samtools view -b@ {threads} '\
        '| sambamba sort -t {threads} -o {output.BAM} --tmpdir {params.TEMP_DIR} /dev/stdin'



rule map_antisense:
    input:
        FQ1 = TRANSCRIPTOME_FILTERED_FQ.replace('{READ_END}','1'),
        FQ2 = TRANSCRIPTOME_FILTERED_FQ.replace('{READ_END}','2')
    
    threads: THREADS
    params:
        INDEX = SMALL_RNA_INDEX,
        RG = '{SAMPLENAME}',
        TEMP_DIR = SAMPLE_FOLDER_TEMPLATE

    output:
        BAM = REVERSE_BAM
    
    shell:
        'bowtie2  ' \
        '--very-sensitive-local --mm ' \
        '-L 8  --mp 4,2 -N 1 '\
        '--no-mixed --no-discordant --dovetail '\
        '-p {threads} -x {params.INDEX} ' \
        '-1 {input.FQ1} -2 {input.FQ2}' \
        '| samtools view -b@ {threads} '\
        '| samtools addreplacerg -@ {threads} -r ID:{params.RG} -r SM:{params.RG} -o - - ' \
        '| samtools view -b@ {threads} '\
        '| samtools sort -O bam -@ {threads} -T {params.TEMP_DIR} -o {output.BAM} ' 


rule filter_genome:
    input:
        FQ1 = TRANSCRIPTOME_FILTERED_FQ.replace('{READ_END}','1'),
        FQ2 = TRANSCRIPTOME_FILTERED_FQ.replace('{READ_END}','2'),

    threads: THREADS
    params:
        INDEX = GENOME_INDEX,

    output:
        FQ1 = GENOME_FILTERED_FQ.replace('{READ_END}','1'),
        FQ2 = GENOME_FILTERED_FQ.replace('{READ_END}','2'),
        BAM = GENOME_BAM


    shell:
        'hisat2 -p {threads} --mm --no-mixed --no-discordant '\
        '-x {params.INDEX} -1 {input.FQ1} -2 {input.FQ2}' \
        '| samtools view -b@ {threads} '\
        '| tee {output.BAM}'\
        '| samtools view -bf4 -@ {threads}'\
        '| samtools fixmate -@ {threads} - -'\
        '| samtools fastq -f0x1 -f0x4 -@ {threads} -1 {output.FQ1} -2 {output.FQ2} - '


rule filter_transcriptome:
    input:
        FQ1 = ANTISENSE_FILTER_CONTAM_FQ.replace('{READ_END}','1'),
        FQ2 = ANTISENSE_FILTER_CONTAM_FQ.replace('{READ_END}','2')

    threads: THREADS
    params:
        INDEX = TRANSCRIPTOME_INDEX,
        GENOME_INDEX = GENOME_INDEX,

    output:
        FQ1 = TRANSCRIPTOME_FILTERED_FQ.replace('{READ_END}','1'),
        FQ2 = TRANSCRIPTOME_FILTERED_FQ.replace('{READ_END}','2'),
        BAM = TRANSCRIPTOME_BAM

    shell:
        'bowtie2 -p {threads} --mm --no-mixed --no-discordant --norc '\
        '-x {params.INDEX} -1 {input.FQ1} -2 {input.FQ2}' \
        '| samtools view -b@ {threads} '\
        '| python ~/ngs_qc_plot/exogenous_filter.py '\
        '-i - -o - --nm 3 -x {params.GENOME_INDEX}'
        '| tee {output.BAM}'\
        '| samtools view -bf4 -@ {threads}'\
        '| samtools fixmate -@ {threads} - -'\
        '| samtools fastq -f0x1 -f0x4 -@ {threads} -1 {output.FQ1} -2 {output.FQ2} - '


rule filter_r1_recopy:
    input:
        FQ1 = ANTISENSE_FQ.replace('{READ_END}','1'),
        FQ2 = ANTISENSE_FQ.replace('{READ_END}','2')

    threads: THREADS
    params:
        R2R_ADAPTER = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTT'

    output:
        FQ1 = ANTISENSE_FILTER_CONTAM_FQ.replace('{READ_END}','1'),
        FQ2 = ANTISENSE_FILTER_CONTAM_FQ.replace('{READ_END}','2')

    shell:
        'seqtk mergepe {input.FQ1} {input.FQ2} '\
        '| cutadapt  -g {params.R2R_ADAPTER} '\
        '-O 4 --discard-trimmed  '\
        '--quiet --interleaved - ' \
        '| cutadapt '\
        '-a AAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC '\
        '-A GATCGTCGGACTGTAGAACTCTGAACGTGTAGAT '\
        ' -O 4 --quiet --interleaved - '\
        '| deinterleave_fastq.py -i - -1 {output.FQ1}  -2 {output.FQ2}'

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
        FQ = expand(ANTISENSE_FQ, SAMPLENAME = SAMPLENAMES, READ_END = ['1'])
    
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
        DATA_PATH = DATA_PATH,
        FQ_NAME = lambda w: w.SAMPLENAME.replace('_R1_001','_R%s_001' %str(w.READ_END))

    output:
        FQ = ANTISENSE_FQ
    
    shell:
        'cat  {params.DATA_PATH}/{params.FQ_NAME}.fastq.gz'\
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


