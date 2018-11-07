import glob
import os
from collections import deque
import re

    

project_path = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
SAMPLES = glob.glob(project_path + '/*R1_001')
SAMPLES = filter(lambda x: re.search('_L[0-9]+|Frag|Phos|[qQ][cC][fF][0-9]+|[ED][DE]|Exo', x), SAMPLES)
SAMPLES = map(lambda x: os.path.basename(x.replace('_R1_001','')), SAMPLES)
SAMPLES = list(SAMPLES)


BAM_template = project_path + '/{SAMPLE}_R1_001/Combined/primary_no_sncRNA_tRNA_rRNA.bam'
FQ1_template = BAM_template.replace('.bam', '.1.fq.gz')
FQ2_template = BAM_template.replace('.bam', '.2.fq.gz')

#SALMON
SALMON_template = project_path + '/salmon_result/{SAMPLE}'
SALMON_REF = os.environ['REF'] + '/hg19/new_genes/salmon_proteins_idx'

#KALLISTO
KALLISTO_template = project_path + '/kallisto_result/{SAMPLE}'
KALLISTO_REF = os.environ['REF'] + '/hg19/new_genes/kallisto_proteins_idx'
KALLISTO_TREATMENT_BAM = project_path + '/kallisto_result/bam_files/{TREATMENT}_kallisto.bam'
KALLISTO_SAMPLE_BAM = KALLISTO_template + '/pseudoalignments.bam'
KALLISTO_READGROUP_SAMPLE_BAM = KALLISTO_template + '/pseudoalignments_rg.bam'

GENE_MAP_PROTEIN = os.environ['REF'] + '/hg19/new_genes/proteins.gtf',
GENOME = os.environ['REF'] + '/hg19/genome/hg19_genome.genome'
num_threads = 4
def get_LIBTYPE(wildcards):
    lt = 'IU' if re.search('L[12]', wildcards.SAMPLE) else 'ISF'
    print(wildcards.SAMPLE, lt) 
    return lt

TREATMENTS = ['polyA','unfragmented','fragmented','phosphatase']
TREATMENTS_regexes = ['_L[0-9]+','[Qq][cC][fF][0-9]+|[ED][ED]|Exo', '[fF]rag[0-9]+', '_Phos']
TREATMENTS_regex_dict = {t:tr for t, tr in zip(TREATMENTS, TREATMENTS_regexes)}
def regex_samples(w):
    regex = TREATMENTS_regex_dict[w.TREATMENT]
    sample = filter(lambda x: re.search(regex, x), SAMPLES)
    return list(sample)


rule all:
    input: 
#        expand(SALMON_template, SAMPLE = SAMPLES),
        expand(KALLISTO_TREATMENT_BAM, TREATMENT = TREATMENTS),

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


rule salmon_quant:
    input:
        FQ1 = FQ1_template,
        FQ2 = FQ2_template,

    params:
        LT = lambda wildcards, output: 'IU' if re.search('L[12]', wildcards.SAMPLE) else 'ISF',
        THREADS = num_threads,
        GENE_MAP_PROTEIN = GENE_MAP_PROTEIN,
        INDEX = SALMON_REF
    
    output: 
        SALMON_OUT = SALMON_template,
    
    shell:
        'salmon quant --mates1 {input.FQ1} --mates2 {input.FQ2} '\
        '-o {output.SALMON_OUT} --threads {params.THREADS} ' \
        '--index {params.INDEX} '\
        '--gcBias --seqBias '\
        '--geneMap {params.GENE_MAP_PROTEIN} '\
        '--libType {params.LT} '


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
        BAM = KALLISTO_SAMPLE_BAM.replace('{SAMPLE}','{SAMPLE, [a-zA-Z0-9_]+}')
    
    shell:
        'kallisto quant --bias '\
        '--output-dir {params.OUT} '\
        '--index {params.INDEX} ' \
        '--threads {params.THREADS} ' \
        '--gtf {params.GENE_MAP_PROTEIN} '\
        '--genomebam --chromosomes {params.GENOME} '\
        ' {input.FQ1} {input.FQ2} '


rule kallisto_merge:
    input:
        BAMS = lambda w: expand(KALLISTO_READGROUP_SAMPLE_BAM, SAMPLE = regex_samples(w))

    params:
        THREADS = num_threads

    output:
        BAM = KALLISTO_TREATMENT_BAM
        
    shell:
        'samtools merge -@ {params.THREADS} {output.BAM} {input.BAMS}'\
        '; samtools index -@ {params.THREADS} {output.BAM}'
