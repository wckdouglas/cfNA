import re
import glob
import os
from collections import Counter


wildcard_constraints:
    TREATMENT = '[a-zA-Z-_0-9]+',
    STRAND = 'antisense|sense',
    PMSTRAND = '[plusmin]+'

PROJECT_PATH = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
SAMPLE_FOLDERS = glob.glob(PROJECT_PATH + '/*001') 
#SAMPLE_FOLDERS = list(filter(lambda x: not re.search('genome-sim', x), SAMPLE_FOLDERS))
SAMPLE_NAMES = list(map(os.path.basename, SAMPLE_FOLDERS))
COMBINED_BAM_PATH = PROJECT_PATH + '/merged_bam'
FILTER_BAM_PATH = COMBINED_BAM_PATH + '/filtered_bam'
REFFLAT = '/stor/work/Lambowitz/ref/hg19/new_genes/proteins.refflat'
snc_annotation = os.environ['REF'] + '/hg19_ref/genes/sncRNA_rRNA_for_bam_filter.bed'
rmsk_annotation = os.environ['REF'] + '/hg19/genome/rmsk.bed'
protein_bed = os.environ['REF'] + '/hg19_ref/genes/protein.bed'
stranded_bed = os.environ['REF'] + '/hg19/new_genes/protein_{PMSTRAND}.bed'
chrM = os.environ['REF'] + '/hg19/genome/chrM.fa'
THREADS = 1

# set up templates
#STRAND: sense, antisense
#PMSTRAND: plus, minus

STRANDED_METRICS_TEMPLATE = FILTER_BAM_PATH + '/{TREATMENT}.{STRAND}.RNA_Metrics'
COMBINED_METRICS_TEMPLATE = FILTER_BAM_PATH + '/{TREATMENT}.RNA_Metrics'
COMBINED_BAM_TEMPLATE = COMBINED_BAM_PATH + '/{TREATMENT}.bam'
COMBINED_NAME_SORT_BAM_TEMPLATE = COMBINED_BAM_PATH + '/{TREATMENT}.nameSorted.bam'
COMBINED_FILTERED_BAM_TEMPLATE = FILTER_BAM_PATH + '/{TREATMENT}.protein.bam'
FILTERED_STRAND_BAM_TEMPLATE = FILTER_BAM_PATH + '/{TREATMENT}.protein.{STRAND}.bam'
FILTERED_PMSTRAND_BAM_TEMPLATE = FILTER_BAM_PATH + '/{TREATMENT}.{PMSTRAND}_{STRAND}.bam'
PMSTRAND_BAM_TEMPLATE = FILTER_BAM_PATH + '/{TREATMENT}.{PMSTRAND}.bam'
chrM_FILTERED_BAM = COMBINED_BAM_PATH + '/dedup/unfragmented.chrM_filter.bam'
PAIRED_DEDUP_BAM = COMBINED_BAM_PATH + '/dedup/unfragmented.chrM_filter.dedup.bam'

SAMPLE_FOLDER = PROJECT_PATH + '/{SAMPLE}'
PICARD_FOLDER = SAMPLE_FOLDER + '/picard'
SAMPLE_PRIMARY_BAM = SAMPLE_FOLDER + '/Combined/primary.bam'
SAMPLE_DEDUP_BAM = SAMPLE_FOLDER + '/Combined/primary.deduplicated.bam'
SAMPLE_MARKDUP_BAM = SAMPLE_FOLDER + '/Combined/primary.mark_duplicated.bam'
SAMPLE_SORTED_BAM = SAMPLE_FOLDER + '/Combined/primary.sorted.bam'
SAMPLE_NAME_SORT_BAM = SAMPLE_FOLDER + '/Combined/primary.mark_duplicated.name_sorted.bam'
SAMPLE_FILTERED_STRAND_BAM_TEMPLATE = PICARD_FOLDER + '/protein.{STRAND}.bam'
SAMPLE_FILTERED_PMSTRAND_BAM_TEMPLATE = PICARD_FOLDER + '/protein.{PMSTRAND}_{STRAND}.bam'
SAMPLE_PMSTRAND_BAM_TEMPLATE = PICARD_FOLDER + '/protein.{PMSTRAND}.bam'
SAMPLE_METRIC_TEMPLATE = PICARD_FOLDER + '/protein.{STRAND}.RNA_Metrics'


# for combining samples
TREATMENT_REGEX = ['Q[Cc][Ff][0-9]+|[ED][DE]|Exo|HS', 'Frag','[pP]hos', 
                  'L[1234]','All','N[aA][0-9]+',
                  'ED|DE','HS[123]','genome',
                    'MPF4','MPF10','MPCEV','^GC',
                    'PPF4','PPF10','PPCEV']
TREATMENTS = ['unfragmented','fragmented','phosphatase',
                'polyA','untreated', 'alkaline_hydrolysis',
                'exonuclease','high_salt','genome-sim',
                'EV','RNP','RNP-EV','HEK293',
                'MNase_EV','MNase_RNP','MNase_EV-RNP'] 
treatment_regex_dict = {t:tr for t, tr in zip(TREATMENTS, TREATMENT_REGEX)}
def select_sample(wildcards, return_count = False):
    regex = treatment_regex_dict[wildcards.TREATMENT] 
    selected_samples = filter(lambda x: re.search(regex, x), SAMPLE_NAMES)
    selected_samples = list(selected_samples)
    if return_count:
        return len(selected_samples)
    else:
        return selected_samples


def get_bams(wildcards):
    samples = select_sample(wildcards)
    #bams = [SAMPLE_SORTED_BAM.format(SAMPLE = s) for s in samples]
    bams = [SAMPLE_DEDUP_BAM.format(SAMPLE = s) for s in samples]
    return bams

def get_preprocessing(wildcards):
    if re.search('genome-sim|[pP]olyA|PEV', wildcards.SAMPLE):
        command = ''
    else:
        command =  '| bam_umi_tag.py -i - -t RX '
    return command

# for strand filtering
def get_filter_command(wildcards):
    if wildcards.PMSTRAND == "plus":
        filtering = ' awk \'{if (($1~/^@/)|| ($2==99)||($2==147)|| ($2==355)||($2==403)|| '\
                        '($2==1123)||($2==1171)) print}\''
    elif wildcards.PMSTRAND == "minus":
        filtering = 'awk \'{if (($1~/^@/)|| ($2==83)||($2==163)|| ($2== 339)||($2==419)|| '\
                            '($2==1187)||($2==1107)) print}\''
    return  filtering


# for deduplication
def get_dedup_command(wildcards):
    if not re.search('genome|[pP]olyA|L[0-9E]+|PEV_|^GC', wildcards.SAMPLE):
        UMI_METRIC = SAMPLE_DEDUP_BAM.replace('.bam','.umi_metrics').format(SAMPLE=wildcards.SAMPLE)
        md = 'picard UmiAwareMarkDuplicatesWithMateCigar '\
            ' UMI_METRICS_FILE={UMI_METRIC} '\
            ' MAX_EDIT_DISTANCE_TO_JOIN=1 '\
            ' TAG_DUPLICATE_SET_MEMBERS=true ' \
            ' UMI_TAG_NAME=RX '\
            'ASSUME_SORT_ORDER=coordinate ' \
            .format(UMI_METRIC = UMI_METRIC)
    
    else:
        md = 'picard MarkDuplicates ASSUME_SORT_ORDER=coordinate ' 

    return md 

# for strand selections
def get_strand_bed(wildcards):
    if wildcards.STRAND == 'sense' and wildcards.PMSTRAND == "plus":
        STRANDED_BED = stranded_bed.format(PMSTRAND = 'plus')
    
    elif wildcards.STRAND == 'antisense' and wildcards.PMSTRAND == 'plus':
        STRANDED_BED = stranded_bed.format(PMSTRAND = 'minus')

    elif wildcards.STRAND == 'sense' and wildcards.PMSTRAND == 'minus':
        STRANDED_BED = stranded_bed.format(PMSTRAND = 'minus')

    elif wildcards.STRAND == 'antisense' and wildcards.PMSTRAND == 'minus':
        STRANDED_BED = stranded_bed.format(PMSTRAND = 'plus')
    return STRANDED_BED

    


#shared command
run_RNASeqMetrics = 'picard CollectRnaSeqMetrics ' \
            'I={input.BAM} '  \
            'O={output.METRIC} '  \
            'REF_FLAT={params.REFFLAT} '  \
            'STRAND=FIRST_READ_TRANSCRIPTION_STRAND AS=false'

run_StrandCombination = 'samtools cat {input.BAMS} '\
        '| samtools view -bF 1024 -@ {params.THREADS} '\
        '| sambamba sort -t {params.THREADS} --tmpdir={params.TMPDIR} '\
        ' -o {output.BAM}  /dev/stdin'

run_ProteinFilter =  'bedtools pairtobed -abam {input.BAM} -b {params.SNC_ANNOTATION} -type neither  '\
        '| bedtools pairtobed -abam - -b {params.RMSK_ANNOTATION} -type neither '\
        '| bedtools pairtobed -abam - -b {params.STRANDED_BED} > {output.BAM} '


run_SplitStrand = 'samtools view -h@ {params.THREADS} {input.BAM} | {params.FILTER_COMMAND}' \
        '| samtools view -b@ {params.THREADS} - > {output.PMSTRAND_BAM} '

run_NameSort = 'sambamba sort -t {params.THREADS} -n --tmpdir={params.TMPDIR} '\
        ' -o /dev/stdout {input.BAM}'\
        '| samtools fixmate -@ {params.THREADS} - {output.NAME_SORT_BAM} '


rule all:
    input:
        expand(STRANDED_METRICS_TEMPLATE, 
                TREATMENT = TREATMENTS,
                STRAND = ['sense', 'antisense']),
        expand(SAMPLE_METRIC_TEMPLATE,
            SAMPLE = SAMPLE_NAMES,
            STRAND = ['sense', 'antisense']),
        expand(COMBINED_METRICS_TEMPLATE,
                TREATMENT = TREATMENTS),
        PAIRED_DEDUP_BAM


rule dedup_umi:
    input:
        BAM = chrM_FILTERED_BAM

    output:
        PAIRED_DEDUP_BAM

    shell:
        'umi_tools dedup -I {input.BAM} '\
        '-S {output} --paired --extract-umi-method tag '\
        '--umi-tag RX --cell-tag RG'


rule chrM_filter_bam:
    input:
        BAM = expand(COMBINED_BAM_TEMPLATE, TREATMENT = ['unfragmented'])

    params:
        INDEX = chrM

    output:
        BAM = chrM_FILTERED_BAM

    shell:
        'python ~/cfNA/peak_callings/chrM_filter.py '\
        '-i {input.BAM} -o {output.BAM} '\
        '-x {params.INDEX} '


rule RNAseqPICARD:
    input:
        BAM = COMBINED_FILTERED_BAM_TEMPLATE

    params:
        REFFLAT = REFFLAT

    output:
        METRIC = COMBINED_METRICS_TEMPLATE

    shell:
        run_RNASeqMetrics
    


rule make_protein_bam:
    input:
        BAMS = expand(FILTERED_STRAND_BAM_TEMPLATE\
                        .replace('{TREATMENT}','{{TREATMENT}}'),
                    STRAND = ['sense','antisense'])
    
    params:
        THREADS = THREADS,
        TMPDIR = COMBINED_FILTERED_BAM_TEMPLATE.replace('.bam','')

    output:
        BAM = COMBINED_FILTERED_BAM_TEMPLATE

    shell:
        run_StrandCombination
        

rule Stranded_RNAseqPICARD:
    input:
        BAM = FILTERED_STRAND_BAM_TEMPLATE
    
    params:
        REFFLAT = REFFLAT
    
    output:
        METRIC = STRANDED_METRICS_TEMPLATE

    log:
        STRANDED_METRICS_TEMPLATE.replace('.RNA_Metrics','.log')

    shell:
        run_RNASeqMetrics

rule Combine_strand:
    # combingin plus_sense and minus_sense // plus_antisense and minus_antisense
    input:
        BAMS = expand(FILTERED_PMSTRAND_BAM_TEMPLATE \
                .replace('{STRAND}', '{{STRAND}}') \
                .replace('{TREATMENT}', '{{TREATMENT}}'),
            PMSTRAND = ['plus','minus'])

    params:
        THREADS = THREADS,
        TMPDIR = FILTERED_STRAND_BAM_TEMPLATE.replace('.bam','.log')
    
    output:
        BAM = FILTERED_STRAND_BAM_TEMPLATE
    
    log:
        FILTERED_STRAND_BAM_TEMPLATE.replace('.bam','.log')

    shell:
        run_StrandCombination


rule Filter_protein:
    input:
        BAM = PMSTRAND_BAM_TEMPLATE

    params:
        SNC_ANNOTATION = snc_annotation,
        RMSK_ANNOTATION = rmsk_annotation,
        STRANDED_BED = lambda w: get_strand_bed(w),
        THREADS = THREADS

    output:
        BAM = FILTERED_PMSTRAND_BAM_TEMPLATE 
    
    log:
        FILTERED_PMSTRAND_BAM_TEMPLATE.replace('.bam','.log')

    shell:
        run_ProteinFilter


rule Split_strand:
    input:
        BAM = COMBINED_NAME_SORT_BAM_TEMPLATE

    params:
        THREADS = THREADS,
        FILTER_COMMAND = lambda w: get_filter_command(w),

    output:
        PMSTRAND_BAM = PMSTRAND_BAM_TEMPLATE
    
    log:
        PMSTRAND_BAM_TEMPLATE.replace('.bam','.log')

    shell:
        run_SplitStrand


rule Name_sort:
    input:
        BAM = COMBINED_BAM_TEMPLATE

    params:
        THREADS = THREADS,
        TMPDIR = COMBINED_NAME_SORT_BAM_TEMPLATE.replace('.bam','')

    output:
        NAME_SORT_BAM = COMBINED_NAME_SORT_BAM_TEMPLATE

    log:
        COMBINED_NAME_SORT_BAM_TEMPLATE.replace('.bam','.log')

    shell:
        run_NameSort


rule Combined_bam:
    input:
        BAM_LIST = lambda w: get_bams(w)

    params:
        THREADS = THREADS,
        TMPDIR = COMBINED_BAM_TEMPLATE.replace('.bam','')

    output:
        BAM = COMBINED_BAM_TEMPLATE
    
    log:
        COMBINED_BAM_TEMPLATE.replace('.bam','.log')

    shell:
        'samtools merge -@ {params.THREADS} - {input.BAM_LIST} '\
        '| sambamba sort --tmpdir={params.TMPDIR} -t {params.THREADS} -o {output.BAM} /dev/stdin'
    

rule remove_dup_sample:
    input:
        BAM = SAMPLE_MARKDUP_BAM
    
    output:
        BAM = SAMPLE_DEDUP_BAM,
        BAI = SAMPLE_DEDUP_BAM.replace('.bam','.bam.bai')
    
    shell:
        'samtools view -bF 1024 {input.BAM}' \
        '> {output.BAM} ' \
        '; samtools index {output.BAM}'


rule markdup_bam_sample:
    input:
        BAM = SAMPLE_SORTED_BAM

    params:
        DEDUP_COMMAND = lambda w: get_dedup_command(w)

    output:
        MARKDUP_BAM = SAMPLE_MARKDUP_BAM,
        DEDUP_METRIC = SAMPLE_DEDUP_BAM.replace('.bam','.dedup_metrics')

    log:
        SAMPLE_DEDUP_BAM.replace('.bam','.log')

    shell:
        '{params.DEDUP_COMMAND} INPUT={input.BAM} '\
        'REMOVE_SEQUENCING_DUPLICATES=true '\
        'OUTPUT={output.MARKDUP_BAM} '\
        'METRICS_FILE={output.DEDUP_METRIC} '\
        'REMOVE_DUPLICATES=false '


rule sort_bam_sample:
    input:
        BAM = SAMPLE_PRIMARY_BAM

    params:
        ID = lambda w: w.SAMPLE,
        PREPROCESS = lambda w: get_preprocessing(w),
        TMPDIR = SAMPLE_SORTED_BAM.replace('.bam','')

    output:
        BAM = SAMPLE_SORTED_BAM

    log:
        SAMPLE_SORTED_BAM.replace('.bam','.log')

    shell:
        'cat {input.BAM} '\
        '{params.PREPROCESS} '\
        '| samtools view -bF 256 -F4 -F2048 '\
        '| samtools addreplacerg -r ID:{params.ID} -r SM:{params.ID}  - '\
        '| samtools view -b '\
        '| sambamba sort -n --tmpdir={params.TMPDIR}_1 -o /dev/stdout /dev/stdin'\
        '| picard FixMateInformation ADD_MATE_CIGAR=true '\
        ' ASSUME_SORTED=true INPUT=/dev/stdin OUTPUT=/dev/stdout ' \
        '| sambamba sort --tmpdir={params.TMPDIR}_2 -o {output.BAM} /dev/stdin'

rule RNAseqPICARD_sample:
    input:
        BAM = SAMPLE_FILTERED_STRAND_BAM_TEMPLATE 
    
    params:
        REFFLAT = REFFLAT
    
    output:
        METRIC = SAMPLE_METRIC_TEMPLATE

    log:
        SAMPLE_METRIC_TEMPLATE.replace('.RNA_Metrics','.log')
    
    shell:
        run_RNASeqMetrics


rule Combine_strand_sample:
    input:
        BAMS = expand(SAMPLE_FILTERED_PMSTRAND_BAM_TEMPLATE\
                .replace('{STRAND}', '{{STRAND}}')\
                .replace('{SAMPLE}', '{{SAMPLE}}'),
            PMSTRAND = ['plus','minus'])
    
    params:
        THREADS = THREADS,
        TMPDIR = SAMPLE_FILTERED_STRAND_BAM_TEMPLATE.replace('.bam','')

    output:
        BAM = SAMPLE_FILTERED_STRAND_BAM_TEMPLATE 
    
    log:
        SAMPLE_FILTERED_STRAND_BAM_TEMPLATE.replace('.bam','.log')
    shell:
        run_StrandCombination



rule Filter_protein_sample:
    input:
        BAM = SAMPLE_PMSTRAND_BAM_TEMPLATE
    
    params:
        SNC_ANNOTATION = snc_annotation,
        RMSK_ANNOTATION = rmsk_annotation,
        STRANDED_BED = lambda w: get_strand_bed(w),
        THREADS = THREADS

    output:
        BAM = SAMPLE_FILTERED_PMSTRAND_BAM_TEMPLATE
    
    log:
        SAMPLE_FILTERED_PMSTRAND_BAM_TEMPLATE.replace('.bam','.log')
    shell:
        run_ProteinFilter

rule Split_strand_sample:
    input:
        BAM = SAMPLE_NAME_SORT_BAM
    
    params:
        THREADS = THREADS,
        FILTER_COMMAND = lambda w: get_filter_command(w),

    output:
        PMSTRAND_BAM = SAMPLE_PMSTRAND_BAM_TEMPLATE
    
    log:
        SAMPLE_PMSTRAND_BAM_TEMPLATE.replace('.bam','.log')
    shell:
        run_SplitStrand


rule Name_sort_sample:
    input:
        BAM = SAMPLE_MARKDUP_BAM
    
    params:
        THREADS = THREADS,
        TMPDIR = SAMPLE_NAME_SORT_BAM.replace('.bam','')
    
    output:
        NAME_SORT_BAM = SAMPLE_NAME_SORT_BAM
    
    log:
        SAMPLE_NAME_SORT_BAM.replace('.bam','.log')

    shell:
        run_NameSort
