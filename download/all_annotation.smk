import os
REF_PATH = os.environ['REF'] + '/hg19_ref'
GENOME_PATH = REF_PATH + '/genome'
GENE_PATH = REF_PATH + '/genes'
DASHR_URL='http://dashr2.lisanwanglab.org/downloads/dashr.v2.sncRNA.annotation.hg19.bed'
DASHR_BED = GENE_PATH + '/dashr.bed'
RBP_BED = GENE_PATH + '/RBP.bed'
RBP_BED = GENE_PATH + '/RBP_IDR.reformated.bed.gz'
REFSEQ = GENE_PATH + '/hg19_refseq.bed'
REFSEQ_SMALL = GENE_PATH + '/hg19_refseq.sncRNA.bed'
ALL_ANNOT = GENE_PATH + '/all_annotation.bed.gz'

rule all:
    input:
        ALL_ANNOT + '.tbi',
        REFSEQ_SMALL


rule dashr:
    input:

    params:
        URL = DASHR_URL
    output:
        DASHR_BED

    shell:
        'curl {params.URL} ' \
        '| python dashr_db.py '\
        '> {output}'

rule refseq:
    input:

    output:
        REFSEQ

    shell:
        'python ucsc_Refseq.py {output}'

rule index:
    input:
        ALL_ANNOT

    output:
        ALL_ANNOT + '.tbi'

    shell:
        'tabix -p bed -f {input}'


rule all_annotation:
    input:
        RBP = RBP_BED,
        DASHR = DASHR_BED,
        GENES = GENE_PATH + '/genes.bed',
        REFSEQ = REFSEQ,
        RMSK = GENE_PATH + '/rmsk.bed',
        piRNA = GENE_PATH + '/piRNA.bed',
        

    output:
        ALL_ANNOT

    shell:
        #'cat {input.RBP} '\
        #"| awk '$8 > 2' "\
        #"| tr '_' '\\t' "\
        #"| awk '{{print $1,$2,$3,$4,$10,$8, \"RBP\",$4}}' OFS='\\t' "\
        "zcat {input.RBP} "
        "| cat - {input.GENES} "\
        ' {input.RMSK} {input.REFSEQ} {input.piRNA} {input.DASHR} '\
        '| sort -k1,1 -k2,2n -k3,3n ' \
        '| bgzip '\
        '> {output}'


rule refseq_small:
    input:
        REFSEQ

    output:
        REFSEQ_SMALL
        
    shell:
        'cat {input}'\
        "| egrep -v 'antisense|protein|lncRNA' "\
        "| awk '($3-$2) < 1000' "\
        "> {output} "
