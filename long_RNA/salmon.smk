import glob
import os
from collections import deque
import re


project_path = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
samples = glob.glob(project_path + '/*R1_001')
samples = filter(lambda x: re.search('L[12]|Frag|Phos|[qQ][cC][fF][0-9]+|[ED][DE]|Exo', x), samples)
samples = list(map(lambda x: os.path.basename(x.replace('_R1_001','')), samples))

BAM_template = project_path + '/{sample}_R1_001/Combined/primary_no_sncRNA_tRNA_rRNA.bam'
FQ1_template = BAM_template.replace('.bam', '.1.fq.gz')
FQ2_template = BAM_template.replace('.bam', '.2.fq.gz')
SALMON_template = project_path + '/salmon_result/{sample}'
SALMON_BAM = project_path + '/salmon_result/bam_files/{sample}.bam'
GENE_MAP_PROTEIN = os.environ['REF'] + '/hg19/new_genes/proteins.gtf',
SALMON_REF = os.environ['REF'] + '/hg19/new_genes/salmon_proteins_idx'
num_threads = 12
def get_LIBTYPE(wildcards):
    lt = 'IU' if re.search('L[12]', wildcards.sample) else 'ISF'
    print(wildcards.sample, lt) 
    return lt



rule all:
    input: 
        expand(SALMON_template, sample = samples) +\
        expand(SALMON_BAM, sample = samples) 

rule make_fq:
    input:
        bam = BAM_template

    output:
        FQ1 = FQ1_template,
        FQ2 = FQ2_template
    threads: num_threads
    shell:
        'samtools fastq -@ {threads} -1 {output.FQ1} -2 {output.FQ2} {input.bam}'


rule salmon_quant:
    input:
        FQ1 = FQ1_template,
        FQ2 = FQ2_template,
        GENE_MAP_PROTEIN = GENE_MAP_PROTEIN,
        INDEX = SALMON_REF

    params:
        LT = lambda wildcards, output: 'IU' if re.search('L[12]', wildcards.sample) else 'ISF'
    
    output: 
        SALMON_OUT = SALMON_template,
        BAM = SALMON_BAM
    threads: num_threads
    
    shell:
        'salmon quant --mates1 {input.FQ1} --mates2 {input.FQ2} '\
            '-o {output.SALMON_OUT} --threads {threads} ' \
            '--index {input.INDEX} '\
            '--gcBias --seqBias '\
            '--geneMap {input.GENE_MAP_PROTEIN} '\
            '--libType {params.LT} '\
            '--writeMappings '\
            '> {output.BAM} '
