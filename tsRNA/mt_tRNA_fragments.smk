import glob
import os
import re

PROJECT_PATH = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
SAMPLE_FOLDERS = glob.glob(PROJECT_PATH + '/*')
SAMPLENAMES = map(os.path.basename, SAMPLE_FOLDERS)
SAMPLENAMES = filter(lambda x: re.search('001$',x), SAMPLENAMES)
SAMPLENAMES = list(SAMPLENAMES)

GENE_REF = os.environ['REF'] + '/hg19_ref/genes/rRNA_mt.bed'
GENE_FA = GENE_REF.replace('.bed','.fa')
tRNA_REF = os.environ['REF'] + '/hg19_ref/genes/smallRNA'
MT_tRNA_REF = os.environ['REF'] + '/hg19_ref/genes/mt_tRNA.fa'
MT_tRNA_REF_INDEX = MT_tRNA_REF + '.1.bt2'
MT_tRNA_STRUCTURE = os.environ['REF'] + '/hg19_ref/genes/mt_tRNA.tRNA_ss'
MT_tRNA_SCAN = os.environ['REF'] + '/hg19_ref/genes/mt_tRNA.tRNA_scan'
ANTICODON_TABLE = MT_tRNA_SCAN.replace('.tRNA_scan','.anticodon_annotations.tsv')
SAMPLE_FOLDERS = PROJECT_PATH + '/{SAMPLENAME}'
MT_FOLDER = SAMPLE_FOLDERS + '/rRNA_mt'
MT_BAM = MT_FOLDER + '/aligned.bam'
MT_tRNA_FQ1 = MT_FOLDER + '/mt_tRNA.1.fq'
MT_tRNA_FQ2 = MT_FOLDER + '/mt_tRNA.2.fq'
MT_tRNA_BAM = MT_FOLDER + '/mt_tRNA.bam'
MT_tRNA_BED = MT_FOLDER + '/mt_tRNA.bed.gz'
MT_tRNA_BED_INDEX = MT_tRNA_BED + '.tbi'
MERGED_BAM = PROJECT_PATH + '/merged_bam/mt_tRNA/{TREATMENT}.bam'
SUBSAMPLE_BAM = MERGED_BAM.replace('.bam','.subsample.bam')
BAM_INDEX = SUBSAMPLE_BAM + '.bai'
THREADS = 1

TREATMENT_REGEX = ['Q[Cc][Ff][0-9]+|[ED][DE]|Exo|HS', 'Frag','[pP]hos',
                    'MPF4','MPF10','MPCEV','^GC',
                    'PPF4','PPF10','PPCEV']
TREATMENTS = ['unfragmented','fragmented','phosphatase',
                'EV','RNP','RNP-EV','HEK293',
                'MNase_EV','MNase_RNP','MNase_EV-RNP']


treatment_regex_dict = {t:tr for t, tr in zip(TREATMENTS, TREATMENT_REGEX)}
def select_sample(wildcards, return_count = False):
    regex = treatment_regex_dict[wildcards.TREATMENT]
    selected_samples = filter(lambda x: re.search(regex, x), SAMPLENAMES)
    selected_samples = list(selected_samples)
    if return_count:
        return len(selected_samples)
    else:
        return selected_samples


rule all:
    input:
        expand(MT_tRNA_BED_INDEX, SAMPLENAME = SAMPLENAMES),
        expand(BAM_INDEX, TREATMENT = TREATMENTS),
        ANTICODON_TABLE


rule index_bam:
    input:
        SUBSAMPLE_BAM

    output:
        BAM_INDEX
    
    shell:
        'samtools index {input}'


rule subsample_bam:
    input:
        MERGED_BAM

    threads: 24
    output:
        SUBSAMPLE_BAM

    shell:
        'samtools view -@ {threads} -bs 1.01 {input} '
        '| filter_soft_clip.py -i - --pe '\
        '| python ~/ngs_qc_plot/bam_viz.py '\
        '| samtools sort -@ {threads} -O bam -o {output} '

rule combine_bam:
    input:
        BAMS = lambda w: expand(MT_tRNA_BAM, SAMPLENAME = select_sample(w)) 

    threads: 24
    output:
        MERGED_BAM

    shell:
        'samtools cat {input.BAMS}'\
        '> {output}'


rule tRNA_bed_index:
    input:
        MT_tRNA_BED
    
    output:
        MT_tRNA_BED_INDEX
    
    shell:
        'tabix -p bed {input}'

rule tRNA_bed:
    input:
        MT_tRNA_BAM
    
    output:
        MT_tRNA_BED
    
    shell:
        'cat {input} '\
        '| samtools view -bF 4 ' \
        '| bam_to_bed.py -i - -c '\
        '| sort -k1,1 -k2,2n -k3,3n '\
        '| bgzip '\
        '> {output}'

rule map_tRNA:
    input:
        FQ1 = MT_tRNA_FQ1,
        FQ2 = MT_tRNA_FQ2,
        INDEX = MT_tRNA_REF_INDEX

    threads: THREADS
    params:
        INDEX = MT_tRNA_REF

    output:
        MT_tRNA_BAM

    shell:
        'bowtie2 --very-sensitive-local '\
        '-L 8  --mp 4,2 -N 1 --mm '\
        '--no-mixed --no-discordant --dovetail '\
        '-p {threads} -x {params.INDEX} '\
        '-1 {input.FQ1} -2 {input.FQ2} '\
        '| samtools view -b@ {threads} '\
        '> {output}'


rule make_anticodon:
    input:
        MT_tRNA_REF,
        MT_tRNA_STRUCTURE,
        MT_tRNA_SCAN

    output:
        ANTICODON_TABLE

    shell:
        'python tRNA_anticodon_annot.py mt'
        

rule make_mt_ref_index:
    input:
        MT_tRNA_REF

    output:
        MT_tRNA_REF_INDEX

    shell:
        'bowtie2-build {input} {input}'


rule extract_fq:
    input:
        MT_BAM
    
    params:
        REF_MODEL = GENE_REF

    output:
        FQ1 = MT_tRNA_FQ1,
        FQ2 = MT_tRNA_FQ2
    
    shell:
        'cat {params.REF_MODEL} '\
        '| grep tRNA --color=no '\
        '| bedtools pairtobed -abam {input} -b - -type both '\
        '| samtools fastq - -1 {output.FQ1} -2 {output.FQ2}'

rule extract_mt_tRNA:
    input:
        FA = GENE_FA,
        BED = GENE_REF
    
    output:
        MT_tRNA_REF

    shell:
        'cat {input.BED} '\
        '| grep Mt_tRNA --color=no '\
        "| awk '{{print $1,$2,$3+1,$4,$5,$6,$7,$8 }}' OFS='\\t'"\
        '| bedtools getfasta -name -fi {input.FA} -bed - -s ' \
        '| seqtk seq -U '\
        '| python add_cca.py ' \
        '> {output} '


rule find_anticodon:
    input:
        MT_tRNA_REF

    output:
        SCAN = MT_tRNA_SCAN,
        SS = MT_tRNA_STRUCTURE
        

    shell:
        'tRNAscan-SE -M mammal -o {output.SCAN} -f {output.SS} {input}'



