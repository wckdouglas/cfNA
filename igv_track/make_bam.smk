import re
import glob
import os

PROJECT_PATH = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
SAMPLE_FOLDERS = glob.glob(PROJECT_PATH + '/*001') 
SAMPLE_NAMES = list(map(os.path.basename, SAMPLE_FOLDERS))
COMBINED_BAM_PATH = PROJECT_PATH + '/merged_bam'
FILTER_BAM_PATH = COMBINED_BAM_PATH + '/filtered_bam'
REFFLAT = '/stor/work/Lambowitz/ref/hg19/new_genes/proteins.refflat'
snc_annotation = os.environ['REF'] + '/hg19/new_genes/sncRNA_rRNA_for_bam_filter.bed'
rmsk_annotation = os.environ['REF'] + '/hg19/genome/rmsk.bed'
protein_bed = os.environ['REF'] + '/hg19/new_genes/protein.bed'
stranded_bed = os.environ['REF'] + '/hg19/new_genes/protein_{PMSTRAND}.bed'
THREADS = 6



# set up templates
#STRAND: sense, antisense
#PMSTRAND: plus, minus

COMBINED_METRICS_TEMPLATE = FILTER_BAM_PATH + '/{TREATMENT}.{STRAND}.RNA_Metrics'
COMBINED_BAM_TEMPLATE = COMBINED_BAM_PATH + '/{TREATMENT}.bam'
COMBINED_NAME_SORT_BAM_TEMPLATE = COMBINED_BAM_PATH + '/{TREATMENT}.nameSorted.bam'

FILTERED_STRAND_BAM_TEMPLATE = FILTER_BAM_PATH + '/{TREATMENT}.protein.{STRAND,[a-z]+}.bam'
FILTERED_PMSTRAND_BAM_TEMPLATE = FILTER_BAM_PATH + '/{TREATMENT}.{PMSTRAND,[a-z]+}_{STRAND,[a-z]+}.bam'
PMSTRAND_BAM_TEMPLATE = FILTER_BAM_PATH + '/{TREATMENT}.{PMSTRAND,[a-z]+}.bam'

SAMPLE_FOLDER = PROJECT_PATH + '/{SAMPLE}'
SAMPLE_PRIMARY_BAM = SAMPLE_FOLDER + '/Combined/primary.bam'
SAMPLE_DEDUP_BAM = SAMPLE_FOLDER + '/Combined/primary.deduplicated.bam'
SAMPLE_MARKDUP_BAM = SAMPLE_FOLDER + '/Combined/primary.mark_duplicated.bam'
SAMPLE_SORTED_BAM = SAMPLE_FOLDER + '/Combined/primary.sorted.bam'
SAMPLE_NAME_SORT_BAM = SAMPLE_FOLDER + '/Combined/primary.mark_duplicated.name_sorted.bam'
SAMPLE_FILTERED_STRAND_BAM_TEMPLATE = SAMPLE_FOLDER + '/picard/protein.{STRAND,[a-z]+}.bam'
SAMPLE_FILTERED_PMSTRAND_BAM_TEMPLATE = SAMPLE_FOLDER + '/protein.{PMSTRAND,[a-z]+}_{STRAND,[a-z]+}.bam'
SAMPLE_PMSTRAND_BAM_TEMPLATE = SAMPLE_FOLDER + '/protein.{PMSTRAND,[a-z]+}.bam'
SAMPLE_METRIC_TEMPLATE = SAMPLE_FOLDER + '/picard/protein.{STRAND}.RNA_Metrics'


# for combining samples
TREATMENT_REGEX = ['Q[Cc][Ff][0-9]+|[ED][DE]|Exo|HS', 'Frag', 
                  'L[12]','All','N[aA][0-9]+',
                  'ED|DE','HS[123]','genome']
TREATMENTS = ['unfragmented','fragmented',
                'polyA','untreated', 'alkaline_hydrolysis',
                'exonuclease','high_salt','genome-sim'] 
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
    if re.search('poly|genome', wildcards.TREATMENT):
        bams = [SAMPLE_PRIMARY_BAM.format(SAMPLE = s) for s in samples]
    else:
        bams = [SAMPLE_MARKDUP_BAM.format(SAMPLE = s) for s in samples]

    return bams

    

# for strand filtering
def get_filter_command(wildcards):
    if wildcards.PMSTRAND == "plus":
        filtering = ' awk \'{{if (($1~/^@/)|| ($2==99)||($2==147)|| ($2==355)||($2==403)|| '\
                        '($2==1123)||($2==1171)) print}}\''
    elif wildcards.PMSTRAND == "minus":
        filtering = 'awk \'{{if (($1~/^@/)|| ($2==83)||($2==163)|| ($2== 339)||($2==419)|| '\
                            '($2==1187)||($2==1107)) print}}\''
    return  filtering


# for deduplication
def get_dedup_command(input, output, params, UMI=True):
    if UMI:
        md = 'picard UmiAwareMarkDuplicatesWithMateCigar UMI_METRICS_FILE={UMI_METRIC} '\
            ' MAX_EDIT_DISTANCE_TO_JOIN=1 TAG_DUPLICATE_SET_MEMBERS=true ' \
            ' UMI_TAG_NAME=RX ' \
            .format(UMI_METRIC = params.UMI_METRIC)
    
    else:
        md = 'picard MarkDuplicates ' 

    return md + 'INPUT={BAM} REMOVE_SEQUENCING_DUPLICATES=true '\
            'OUTPUT=/dev/stdout '\
            'METRICS_FILE={DEDUP_METRIC} REMOVE_DUPLICATES=false ASSUME_SORT_ORDER=coordinate '\
            '| tee {MARKDUP_BAM} '\
            '| samtools view -bF 1024 ' \
            '> {DEDUP_BAM} ' \
            '; samtools index {DEDUP_BAM}'\
            .format(BAM = input.BAM, 
                DEDUP_METRIC = output.DEDUP_METRIC,
                MARKDUP_BAM = output.MARKDUP_BAM,
                DEDUP_BAM = output.DEDUP_BAM)
    


#shared command
run_RNASeqMetrics = 'picard CollectRnaSeqMetrics ' \
            'I={input.BAM} '  \
            'O={output.METRIC} '  \
            'REF_FLAT={params.REFFLAT} '  \
            'STRAND=FIRST_READ_TRANSCRIPTION_STRAND AS=false'

run_StrandCombination = 'samtools cat {input.BAMS} '\
        '| samtools view -bF 1024 -@ {params.THREADS} '\
        '| sambamba sort -t {params.THREADS} '\
        '--show-progress -o {output.BAM}  /dev/stdin'

run_ProteinFilter =  '; bedtools pairtobed -abam {input.BAM} -b {params.SNC_ANNOTATION} -type neither  '\
        '| bedtools pairtobed -abam - -b {params.RMSK_ANNOTATION} -type neither '\
        '| bedtools pairtobed -abam - -b {params.STRANDED_BED} > {output.BAM} '


run_SplitStrand = 'samtools view -h@ {params.THREADS} {input.BAM} | {params.FILTER_COMMAND}' \
        '| samtools view -b@ {params.THREADS} - > {output.PMSTRAND_BAM} '

run_NameSort = 'sambamba sort -t {params.THREADS} -n '\
        '--show-progress -o /dev/stdout {input.BAM}'\
        '| samtools fixmate -@ {params.THREADS} - {output.NAME_SORT_BAM} '


rule all:
    input:
        expand(COMBINED_METRICS_TEMPLATE, 
                TREATMENT = TREATMENTS,
                STRAND = ['sense', 'antisense']),
        expand(SAMPLE_METRIC_TEMPLATE,
            SAMPLE = SAMPLE_NAMES,
            STRAND = ['sense', 'antisense'])


rule RNAseqPICARD:
    input:
        BAM = FILTERED_STRAND_BAM_TEMPLATE\
            .replace('{TREATMENT}', '{TREATMENT,[a-zA-Z-_]+}')
    
    params:
        REFFLAT = REFFLAT
    
    output:
        METRIC = COMBINED_METRICS_TEMPLATE\
            .replace('{TREATMENT}', '{TREATMENT,[a-zA-Z-_]+}')

    shell:
        run_RNASeqMetrics

rule Combine_strand:
    # combingin plus_sense and minus_sense // plus_antisense and minus_antisense
    input:
        BAMS = expand(FILTERED_PMSTRAND_BAM_TEMPLATE\
                .replace('{STRAND,[a-z]+}', '{{STRAND}}')\
                .replace('{TREATMENT}', '{{TREATMENT}}')\
                .replace('{PMSTRAND,[a-z]+}', '{PMSTRAND}'),
            PMSTRAND = ['plus','minus'])

    params:
        THREADS = THREADS

    output:
        BAM = FILTERED_STRAND_BAM_TEMPLATE\
            .replace('{TREATMENT}', '{TREATMENT,[a-zA-Z-_]+}')
    
    shell:
        run_StrandCombination


rule Filter_protein:
    input:
        BAM = PMSTRAND_BAM_TEMPLATE\
            .replace('{TREATMENT}', '{TREATMENT,[a-zA-Z-_]+}')

    params:
        SNC_ANNOTATION = snc_annotation,
        RMSK_ANNOTATION = rmsk_annotation,
        STRANDED_BED = lambda w: stranded_bed.format(PMSTRAND =  w.PMSTRAND),
        THREADS = THREADS

    output:
        BAM = FILTERED_PMSTRAND_BAM_TEMPLATE 
    
    shell:
        run_ProteinFilter


rule Split_strand:
    input:
        BAM = COMBINED_NAME_SORT_BAM_TEMPLATE\
            .replace('{TREATMENT}', '{TREATMENT,[a-zA-Z-_]+}')

    params:
        THREADS = THREADS,
        FILTER_COMMAND = lambda w: get_filter_command(w),

    output:
        PMSTRAND_BAM = PMSTRAND_BAM_TEMPLATE\
            .replace('{TREATMENT}', '{TREATMENT,[a-zA-Z-_]+}')
    
    shell:
        run_SplitStrand


rule Name_sort:
    input:
        BAM = COMBINED_BAM_TEMPLATE

    params:
        THREADS = THREADS,

    output:
        NAME_SORT_BAM = COMBINED_NAME_SORT_BAM_TEMPLATE

    shell:
        run_NameSort


rule Combined_bam:
    input:
        BAMS = lambda w: get_bams(w)

    params:
        THREADS = THREADS,

    output:
        BAM = COMBINED_BAM_TEMPLATE
    
    run:
        if select_sample(wildcards, return_count = True) > 1:
            COMBINED_COMMAND = 'sambamba merge -t {params.THREADS} --show-progress  /dev/stdout ' 
        else:
            COMBINED_COMMAND = 'cat '
        
        bam_out = output.BAM
        bams_in = ' '.join(input.BAMS)

        shell(COMBINED_COMMAND + ' %s > %s' %(bams_in, bam_out))
    

rule dedup_bam_sample:
    input:
        BAM = SAMPLE_SORTED_BAM

    params:
        UMI_METRIC = SAMPLE_DEDUP_BAM.replace('.bam','.umi_metrics'),

    output:
        DEDUP_BAM = SAMPLE_DEDUP_BAM,
        MARKDUP_BAM = SAMPLE_MARKDUP_BAM,
        DEDUP_METRIC = SAMPLE_DEDUP_BAM.replace('.bam','.dedup_metrics')

    run:
        DEDUP_COMMAND = get_dedup_command(input, output, params, 
                        UMI = re.search('genome-sim|[pP]olyA', wildcards.SAMPLE))
        shell(DEDUP_COMMAND)


rule sort_bam_sample:
    input:
        BAM = SAMPLE_PRIMARY_BAM

    params:
        ID = lambda w: w.SAMPLE,
        PREPROCESS = lambda w: '' if re.search('genome-sim|[pP]olyA', w.SAMPLE)  else '| bam_umi_tag.py -i - -t RX ' 

    output:
        BAM = SAMPLE_SORTED_BAM

    shell:
        'cat {input.BAM} '\
        '{params.PREPROCESS} '\
        '| samtools view -bF 256 -F4 -F2048 '\
        '| samtools addreplacerg -r ID:{params.ID} -r SM:{params.ID}  - '\
        '| samtools view -b '\
        '| sambamba sort -n -o /dev/stdout /dev/stdin '\
        '| picard FixMateInformation ADD_MATE_CIGAR=true '\
        ' ASSUME_SORTED=true INPUT=/dev/stdin OUTPUT=/dev/stdout ' \
        '| sambamba sort -o {output.BAM} /dev/stdin'

rule RNAseqPICARD_sample:
    input:
        BAM = SAMPLE_FILTERED_STRAND_BAM_TEMPLATE 
    
    params:
        REFFLAT = REFFLAT
    
    output:
        METRIC = SAMPLE_METRIC_TEMPLATE
    
    shell:
        run_RNASeqMetrics


rule Combine_strand_sample:
    input:
        BAMS = expand(SAMPLE_FILTERED_PMSTRAND_BAM_TEMPLATE\
                .replace('{STRAND,[a-z]+}', '{{STRAND}}')\
                .replace('{SAMPLE}', '{{SAMPLE}}')\
                .replace('{PMSTRAND,[a-z]+}', '{PMSTRAND}'),
            PMSTRAND = ['plus','minus'])
    
    params:
        THREADS = THREADS

    output:
        BAM = SAMPLE_FILTERED_STRAND_BAM_TEMPLATE 
    
    shell:
        run_StrandCombination



rule Filter_protein_sample:
    input:
        BAM = SAMPLE_PMSTRAND_BAM_TEMPLATE
    
    params:
        SNC_ANNOTATION = snc_annotation,
        RMSK_ANNOTATION = rmsk_annotation,
        STRANDED_BED = lambda w: stranded_bed.format(PMSTRAND =  w.PMSTRAND),
        THREADS = THREADS

    output:
        BAM = SAMPLE_FILTERED_PMSTRAND_BAM_TEMPLATE
    
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
    
    shell:
        run_SplitStrand


rule Name_sort_sample:
    input:
        BAM = SAMPLE_MARKDUP_BAM
    
    params:
        THREADS = THREADS,
    
    output:
        NAME_SORT_BAM = SAMPLE_NAME_SORT_BAM
    
    run:
        run_NameSort