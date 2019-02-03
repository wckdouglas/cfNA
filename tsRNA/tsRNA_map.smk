import glob
import os
import sys

wildcard_constraints:
    SAMPLENAME = '[Qq][cC][fF][0-9A-Za-z_]+'


PROJECT_PATH = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
SAMPLENAMES = glob.glob(PROJECT_PATH + '/*R1_001')
SAMPLENAMES = map(lambda x: os.path.basename(x).split('_R1')[0], SAMPLENAMES)
SAMPLENAMES = filter(lambda x: re.search('^[Qq][cC][fF]',x), SAMPLENAMES)
SAMPLENAMES = list(set(SAMPLENAMES))

tRF_BLAST_DB = '/stor/work/Lambowitz/ref/hg19_ref/tRNA_fragments/tRNA_frag'
tRNA_BLAST_DB = '/stor/work/Lambowitz/ref/TGSEQ/blast/tRNA_miRNA.fa'
DB_INDEX = '/stor/work/Lambowitz/ref/hg19_ref/tRNA_fragments/tRNA_frag'
THREADS = 4

#SNAKEMAKE VARIABLE
FQ1 = PROJECT_PATH +'/{SAMPLENAME}_R1_001/Trimmed/trim.1.fq.gz'
FQ2 = PROJECT_PATH +'/{SAMPLENAME}_R1_001/Trimmed/trim.2.fq.gz'
PROJECT_PATH = '/stor/work/Lambowitz/cdw2854/cfNA/tsRNA'
SAMPLE_FOLDER = PROJECT_PATH + '/{SAMPLENAME}'
PEAR_FQ = SAMPLE_FOLDER + '/insert.assembled.fastq'
INSERT_FQ = SAMPLE_FOLDER + '/insert.filtered.fastq'
INSERT_FA = SAMPLE_FOLDER + '/insert.filtered.fa'
BLAST_TAB = SAMPLE_FOLDER + '/blast.{REF_TYPE}.tsv'
BLAST_tRNA_TAB = SAMPLE_FOLDER + '/'
BAM = SAMPLE_FOLDER + '/aligned.bam'
REF_TYPES = ['tRF','smRNA']


rule all:
    input:
        expand(BAM, SAMPLENAME = SAMPLENAMES),
        expand(BLAST_TAB, SAMPLENAME = SAMPLENAMES, REF_TYPE = REF_TYPES),

rule bowtie:
    input:
        INSERT_FQ
    
    threads: THREADS
    params:
        INDEX = DB_INDEX

    output:
        BAM
    
    shell:
        'bowtie -k 5 -q --seedmms 3 '\
        '--tryhard  --phred33-quals --sam '\
        '--threads {threads} '
        ' --norc --seedlen 12 {params.INDEX} {input} '\
        '| samtools view -bF 4 -@ {threads} '
        '> {output} '

rule blast:
    input:
        INSERT_FA
    
    threads: THREADS
    params:
        BLAST_DB = lambda w: tRF_BLAST_DB if w.REF_TYPE == "tRF" else tRNA_BLAST_DB
    
    output:
        BLAST_TAB

    shell:
        'blastn  -strand plus -num_threads {threads} '\
        '-query {input}  -db {params.BLAST_DB} '\
        '-evalue 1e-6 -outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue" '\
        '-out {output}'


rule fq2fa:
    input:
        INSERT_FQ

    output:
        INSERT_FA

    shell:
        'cat {input} |seqtk seq -A > {output}'


rule filter_fq:
    input:
        PEAR_FQ

    output:
        INSERT_FQ
    
    shell:
        'cat {input} '\
        '| seqkit seq --max-len 40 --min-len 12 -o {output}'


rule merge_trimmed:
    input:
        FQ1 = FQ1,
        FQ2 = FQ2
    
    threads: THREADS
    params:
        OUT_PREFIX = INSERT_FQ.replace('.assembled.fastq.gz',''),

    output:
        PEAR_FQ
    
    shell:
        'pear -f {input.FQ1} -r {input.FQ2} -o {params.OUT_PREFIX} '\
        ' -n 12 -m 40 -j {threads} -p 0.01 '
