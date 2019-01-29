
import os

wildcard_constraints:
    RNA_TYPE = '[a-z]+'

REF_PATH = os.environ['REF'] + '/hg19_ref' 
GENE_PATH = REF_PATH + '/genes'
GENOME_PATH = REF_PATH + '/genome'
GENOME_FA = GENOME_PATH + '/hg19_genome.fa'
GENES_GTF = GENE_PATH + '/genes.gtf'

#SNAKEMAKE VARIABLE
GENE_FA = GENE_PATH + '/genes_{RNA_TYPE}.fa'
KALLISTO_INDEX = GENE_FA.replace('.fa','.kallisto_idx')
GENE_ANNOTATION = GENE_PATH + '/genes.annot'
RNA_TYPES = ['all','protein']

def extract_gtf_info(l, field = 'gene_id'):
    info = l.split(field)[1].split(';')[0]
    return info.strip('" ')

rule all:
    input:
        expand(KALLISTO_INDEX, RNA_TYPE = RNA_TYPES),
        GENE_ANNOTATION
    
rule make_annot:
    input:

    params:
        GTF = GENES_GTF    

    output:
        TAB = GENE_ANNOTATION
    
    run:
        tdf = {}
        with open(params.GTF) as gtf:
            for line in gtf:
                if not line.startswith('#') and line.split('\t')[2] == 'transcript':
                    tid = extract_gtf_info(line, 'transcript_id')
                    gid = extract_gtf_info(line, 'gene_id')
                    gtype = extract_gtf_info(line, 'gene_type')
                    tdf[tid] = (gid, gtype)

        with open(output.TAB, 'w') as out:
            for tid, gene in tdf.items():
                print('%s\t%s\t%s' %(tid, gene[0], gene[1]), file = out)


rule make_index:
    input:
        FA = GENE_FA
    
    output:
        IDX = KALLISTO_INDEX
    
    shell:
        'kallisto index -i {output.IDX} {input.FA}'

rule make_fa:
    input:

    params:
        FA = GENOME_FA,
        GTF = GENES_GTF,
        FILTER = lambda w: '' if w.RNA_TYPE == "all" else "| grep --color=no 'protein_coding'"
    
    output:
        FA = GENE_FA
    
    shell:
        'cat {params.GTF}'\
        '{params.FILTER} '\
        '| gffread - -g {params.FA} -w {output.FA}'

