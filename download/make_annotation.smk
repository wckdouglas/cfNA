import os
import sys
import re
from itertools import groupby
from ucsc_Refseq import main as download_refseq_function

REF = os.environ['REF']
REF_PATH = REF + '/hg19_ref'
GENOME_PATH = REF_PATH + '/genome'
GENE_PATH = REF_PATH + '/genes'
DASHR_URL = 'http://dashr2.lisanwanglab.org/downloads/dashr.v2.sncRNA.annotation.hg19.bed'
RBP_URL = 'https://www.encodeproject.org/metadata/type=Experiment&assay_title=eCLIP&limit=all/metadata.tsv'
ENCODE_FILE = 'encode_files.txt'
RMSK_URL = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz'
PIRNA_URL = 'http://www.regulatoryrna.org/database/piRNA/download/archive/v1.0/bed/piR_hg19_v1.0.bed.gz'
GTF = GENE_PATH + '/genes.gtf'

#out bed
RBP_TABLE = 'encode_url.txt'
DASHR_BED = GENE_PATH + '/dashr.bed.gz'
RBP_BED = GENE_PATH + '/RBP.bed.gz'
RBP_FILTERED_BED = GENE_PATH + '/RBP.filtered.bed.gz'
SORTED_RBP_BED = RBP_BED.replace('.bed.gz','sorted.bed.gz')
RMSK_BED = GENE_PATH + '/rmsk.bed.gz'
RMSK_SMRNA_BED = GENE_PATH + '/rmsk.smRNA.bed.gz'
REFSEQ_BED = GENE_PATH + '/hg19_refseq.bed.gz'
REFSEQ_SMRNA_BED = GENE_PATH + '/hg19_refseq.smRNA.bed.gz'
GENES_BED = GENE_PATH + '/genes.bed.gz'
PIRNA_BED = GENE_PATH + '/piRNA.bed.gz'
ALL_ANNOTATION = GENE_PATH + '/all_annotation.bed.gz'
ALL_EXONS = GENE_PATH + '/exons_all.bed'


rule all:
    input:
        ALL_ANNOTATION,
        RMSK_SMRNA_BED,
        REFSEQ_SMRNA_BED,
        ALL_EXONS,
        RBP_BED

rule merge_annotation:
    input:
        DASHR = DASHR_BED,
        REFSEQ = REFSEQ_BED,
        RBP = RBP_BED,
        GENE = GENES_BED,
        PIRNA = PIRNA_BED
    
    params:
        RMSK = GENE_PATH + '/rmsk.bed.gz',


    output:
        BED = ALL_ANNOTATION,
        RBP = RBP_FILTERED_BED
    
    shell:
        'zcat {input.RBP} '
        '| python clean_rbp.py ' \
        '| sort -k1,1 -k2,2n -k3,3n ' \
        '| bgzip ' \
        '| tee {output.RBP} ' \
        '| zcat - {input.GENE} {params.RMSK} ' \
        ' {input.DASHR} {input.REFSEQ} {input.PIRNA} ' \
        '| sort -k1,1 -k2,2n -k3,3n ' \
        '| bgzip '\
        '> {output.BED}'\
        '; tabix -f -p bed {output.BED}'\
        '; tabix -f -p bed {output.RBP}'

rule rmsk_smRNA:
    input:
        RMSK_BED
    
    output:
        RMSK_SMRNA_BED

    shell:
        'zcat {input} '\
        "| awk '$7~/^RNA$|tRNA|rRNA|scRNA|srpRNA|snRNA/'"\
        '|bgzip '\
        '> {output}'

rule refseq_smRNA:
    input:
        REFSEQ_BED
    
    output:
        REFSEQ_SMRNA_BED

    shell:
        'zcat {input} '\
        "| awk '$7~/guide_RNA|miRNA|scRNA|snoRNA|misc|snRNA|SRP_RNA|vault|Y_RNA|^rRNA$/'"\
        '| bgzip '\
        '> {output}'



rule gzip_genes:
    input:

    params:
        GENES = GENE_PATH + '/genes.bed'

    output:
        GENES_BED

    shell:
        'cat {params.GENES} '\
        '|sort -k1,1 -k2,2n -k3,3n '\
        '| bgzip > {output}'
        
rule make_piRNA:
    input:

    params:
        LINK = PIRNA_URL
    
    output:
        PIRNA_BED

    shell:
        'curl {params.LINK} '\
        '| zcat '\
        "| awk '$5 > 1000 {{print $0, \"piRNA\",\"piRNA\"}}' OFS='\\t' " \
        '| bgzip '\
        '> {output}'

rule download_refseq:
    input:
    
    output:
        BED = REFSEQ_BED

    run:
        download_refseq_function(output.BED.replace('.gz',''))
        command = 'cat {not_sorted} '\
                '| sort -k1,1 -k2,2n -k3,3n '\
                '| bgzip '\
                '> {sorted} '\
                .format(not_sorted = output.BED.replace('.gz','' ),
                        sorted = output.BED)
        shell(command)


rule make_dashr:
    input:
    
    params:
        LINK = DASHR_URL

    output:
        DASHR_BED

    shell:
        'curl {params.LINK} '\
        '| python dashr_db.py '\
        '| sort -k1,1 -k2,2n -k3,3n '\
        '| bgzip '
        '> {output}'

rule sort_rbp:
    input:
        RBP_BED
    
    params:
        TEMP = SORTED_RBP_BED + '_temp'

    output:
        SORTED_RBP_BED
    
    shell:
        'mkdir -p {params.TEMP} '\
        '; zcat {input} '\
        '| sort -k1,1 -k2,2n -k3,3n -k6,6 -T {params.TEMP} '\ 
        "| bedtools merge -i - -s -c 4 -o collapse  -delim ',' "\
        "| awk '{{print $1,$2,$3,\"RBP\",0,$4,$NF}}' OFS='\\t' "\
        '| bgzip '\
        '\> {output}'\
        '; rm -rf {params.TEMP}'

rule make_rbp:
    input:
        RBP_TABLE
    
    params:
        TEMP = RBP_BED + '_temp'

    output:
        RBP_BED

    shell:
        'rm -f {output}; touch {output};' \
        'for FILE_URL in $(cat {input}| grep -v meta );'\
        "do curl -L $FILE_URL | zcat | awk '$2~/^[0-9]+$/' >> {output};"\
        'done '\
        '; mv {output} {params.TEMP}'\
        '; cat {params.TEMP} | sort -k1,1 -k2,2 -k3,3 | bgzip > {output}'



rule download_rbp_table:
    input:

    params:
        LINK = RBP_URL,
        ENCODE_META = ENCODE_FILE

    output:
        RBP_TABLE

    shell:
        'cat {params.ENCODE_META}'\
        '| grep --color=no bed '\
	    '| grep -v bigBed '\
        '> {output}'

rule make_exons:
    input:

    params:
        GTF = GTF
    output:
        EXONS = ALL_EXONS

    run:
        def extract_gene_name(line, info = 'gene_id'):
            gn = line.split('%s '%info)[1]
            return gn.split(';')[0].strip('"')

        smallRNA = re.compile('miRNA|misc|rRNA|tRNA|snoRNA|snRNA|vault|scRNA')
        with open(params.GTF) as gtf:
            for hash_count, line in enumerate(gtf):
                if not line.startswith('##'):
                    break
        
        temp_file = output.EXONS + '_temp'
        with open(params.GTF) as gtf,\
                open(temp_file, 'w') as exon_out:
            for i in range(hash_count):
                line = next(gtf)

            for gname, lines in groupby(gtf, extract_gene_name):
                lines = [line for line in lines if line.split('\t')[2]=="exon"]
                exon_count = len(lines)
                for line in  lines:
                    if not smallRNA.search(line):
                        fields = line.split('\t')
                        try:
                            outline = '{chrom}\t{start}\t{end}\t{name}\t{exon_count}\t{strand}\t{gene_type}'\
                                .format(chrom = fields[0], start = int(fields[3])-1, end = int(fields[4]) + 1, strand = fields[6],
                                        name = extract_gene_name(line, info = 'gene_name'), exon_count = exon_count,
                                        gene_type = extract_gene_name(line, info ='gene_type'))
                            print(outline, file = exon_out)
                        except IndexError:
                            print(line)
                            sys.exit()

        shell('cat {temp} | sort -k1,1 -k2,2n -k3,3n -k6,6 -u > {EXONS}'.format(temp = temp_file, EXONS = output.EXONS))
