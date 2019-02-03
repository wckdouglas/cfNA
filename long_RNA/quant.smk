
import glob
import os
from collections import deque
import re

wildcard_constraints:
    TREATMENT = '[a-zA-Z_-]+',
    SAMPLE = '[A-Za-z0-9_-]+',
    RNA_TYPE = '[a-z]+'

project_path = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
SAMPLES = glob.glob(project_path + '/*R1_001')
SAMPLES = map(lambda x: os.path.basename(x.replace('_R1_001','')), SAMPLES)
SAMPLES = list(SAMPLES)


BAM_template = project_path + '/{SAMPLE}_R1_001/Combined/primary_no_sncRNA_tRNA_rRNA_repeats.bam'
FQ1_template = BAM_template.replace('.bam', '.1.fq.gz')
FQ2_template = BAM_template.replace('.bam', '.2.fq.gz')

#KALLISTO
REFFLAT = '/stor/work/Lambowitz/ref/hg19/new_genes/proteins.refflat'
KALLISTO_FOLDER = project_path + '/kallisto_{RNA_TYPE}_result'
KALLISTO_template = KALLISTO_FOLDER + '/{SAMPLE}'
KALLISTO_REF = os.environ['REF'] + '/hg19_ref/genes/genes_{RNA_TYPE}.kallisto_idx'
KALLISTO_TREATMENT_BAM = KALLISTO_FOLDER + '/bam_files/{TREATMENT}_kallisto.bam'
KALLISTO_TREATMENT_PICARD = KALLISTO_FOLDER + '/bam_files/picard/{TREATMENT}_kallisto.RNAseq_metrics'
KALLISTO_SAMPLE_BAM = KALLISTO_template + '/pseudoalignments.bam'
KALLISTO_READGROUP_SAMPLE_BAM = KALLISTO_template + '/pseudoalignments_rg.bam'

GENE_MAP_PROTEIN = os.environ['REF'] + '/hg19_ref/genes/genes.gtf',
GENOME = os.environ['REF'] + '/hg19_ref/genome/hg19_genome.genome'
num_threads = 4
def get_LIBTYPE(wildcards):
    lt = 'IU' if re.search('L[12]', wildcards.SAMPLE) else 'ISF'
    print(wildcards.SAMPLE, lt) 
    return lt

RNA_TYPES = ['all','protein']
TREATMENT_REGEX = ['Q[Cc][Ff][0-9]+|[ED][DE]|Exo', '[fF]rag[0-9]+','[pP]hos',
                  'L[1234]','MPF4','MPF10','MPCEV',
                    'PPF4','PPF10','PPCEV']
TREATMENTS = ['unfragmented','fragmented','phosphatase',
                'polyA','EV','RNP','RNP-EV',
                'MNase_EV','MNase_RNP','MNase_EV-RNP']

TREATMENTS_regex_dict = {t:tr for t, tr in zip(TREATMENTS, TREATMENT_REGEX)}
print(TREATMENTS_regex_dict)
def regex_samples(w):
    regex = TREATMENTS_regex_dict[w.TREATMENT]
    sample = list(filter(lambda x: re.search(regex, x), SAMPLES))
    return sample


rule all:
    input: 
        expand(KALLISTO_TREATMENT_PICARD,
                RNA_TYPE = RNA_TYPES, 
                TREATMENT = TREATMENTS),

rule run_picard:
    input:
        KALLISTO_TREATMENT_BAM
    
    params:
        REFFLAT = REFFLAT
    output:
        KALLISTO_TREATMENT_PICARD
    
    shell:
        'picard CollectRnaSeqMetrics ' \
        'I={input} '  \
        'O={output} '  \
        'REF_FLAT={params.REFFLAT} '  \
        'STRAND=FIRST_READ_TRANSCRIPTION_STRAND AS=false'


rule make_fq:
    input:
        bam = BAM_template

    params:
        THREADS = num_threads

    output:
        FQ1 = FQ1_template,
        FQ2 = FQ2_template
    shell:
        'samtools fastq -@ {params.THREADS} -1 {output.FQ1} -2 {output.FQ2} {input.bam}'


rule add_read_group:
    input:
        BAM = KALLISTO_SAMPLE_BAM

    params:
        SAMPLENAME = lambda w: w.SAMPLE,
        THREADS = num_threads

    output:
        BAM = KALLISTO_READGROUP_SAMPLE_BAM

    shell:
        'samtools addreplacerg -r ID:{params.SAMPLENAME} -r SM:{params.SAMPLENAME} '\
        '-o - {input.BAM}' \
        '| samtools sort -O bam -@ {params.THREADS} -o {output.BAM} -'

rule kallisto_quant:
    input:
        FQ1 = FQ1_template,
        FQ2 = FQ2_template,

    params:
        LT = lambda wildcards, output: '' if re.search('L[12]', wildcards.SAMPLE) else '--fr-stranded',
        THREADS = num_threads,
        GENE_MAP_PROTEIN = GENE_MAP_PROTEIN,
        INDEX = KALLISTO_REF,
        GENOME = GENOME,
        OUT = KALLISTO_template
    
    output: 
        BAM = KALLISTO_SAMPLE_BAM
    
    shell:
        'kallisto quant --bias '\
        '--output-dir {params.OUT} '\
        '--index {params.INDEX} ' \
        '--threads {params.THREADS} ' \
        '--gtf {params.GENE_MAP_PROTEIN} '\
        '--genomebam --chromosomes {params.GENOME} '\
        ' {params.LT} ' \
        ' {input.FQ1} {input.FQ2} '


rule kallisto_merge:
    input:
        BAMS = lambda w: expand(KALLISTO_READGROUP_SAMPLE_BAM,
                                SAMPLE = regex_samples(w),
                                RNA_TYPE = [w.RNA_TYPE])

    params:
        THREADS = num_threads

    output:
        BAM = KALLISTO_TREATMENT_BAM
        
    shell:
        'samtools merge -@ {params.THREADS} - {input.BAMS}'\
        '| samtools view -h@ {params.THREADS} '\
        "| egrep -v '[0-9]+M0N[0-9]+M' "\
        '| samtools view -bF4 > {output.BAM}'\
        '; samtools index -@ {params.THREADS} {output.BAM}'
