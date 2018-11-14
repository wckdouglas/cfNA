import glob
import re
import os


PROJECT_PATH = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
SAMPLE_FOLDERS = glob.glob(PROJECT_PATH + '/*001')
SAMPLE_NAMES = map(os.path.basename, SAMPLE_FOLDERS)
SAMPLE_NAMES = filter(lambda x: re.search('Q[cC][fF]',x), SAMPLE_NAMES)
SAMPLE_NAMES = list(SAMPLE_NAMES)
SAMPLE_FOLDER = PROJECT_PATH + '/{SAMPLE_NAME}'
REPEAT_FOLDER = SAMPLE_FOLDER + '/repeats'
FQ1 = REPEAT_FOLDER + '/repeats.1.fq.gz'
FQ2 = REPEAT_FOLDER + '/repeats.2.fq.gz'
STR_BAM = REPEAT_FOLDER + '/STR.bam'
STR_BED = REPEAT_FOLDER + '/STR.bed'
DEDUP_STR_BED = REPEAT_FOLDER + '/STR.dedup.bed'
STR_COUNT_TEMPLATE = PROJECT_PATH + '/STR/{TREATMENT}.tsv'
REF_PATH = '/stor/work/Lambowitz/ref/hg19_ref/STR'
REF_INDEX = REF_PATH + '/STRdecoys.fasta'
THREADS = 24

def dedup_command(w):
    if re.search('L[0-9E]+|', w.SAMPLE_NAME):
        command = ' sort -k1,1 -k2,2 -k3,3 -k6,6 -u '
    else:
        command = ' deduplicate_bed.py -i - -t 1 -o - -d "_" -f 0 '
    return command


# for combining samples
TREATMENT_REGEX = ['Q[Cc][Ff][0-9]+|[ED][DE]|Exo|HS', 'Frag','[pP]hos', 
                  'All','N[aA][0-9]+',
                  'ED|DE','HS[123]']
TREATMENTS = ['unfragmented','fragmented','phosphatase',
                'untreated', 'alkaline_hydrolysis',
                'exonuclease','high_salt'] 
treatment_regex_dict = {t:tr for t, tr in zip(TREATMENTS, TREATMENT_REGEX)}
def select_sample(wildcards, return_count = False):
    regex = treatment_regex_dict[wildcards.TREATMENT] 
    selected_samples = filter(lambda x: re.search(regex, x), SAMPLE_NAMES)
    selected_samples = list(selected_samples)
    if return_count:
        return len(selected_samples)
    else:
        return selected_samples

def dedup_command(w):
    if re.search('L[0-9E]+|', w.SAMPLE_NAME):
        command = ' sort -k1,1 -k2,2 -k3,3 -k6,6 -u '
    else:
        command = ' deduplicate_bed.py -i - -t 1 -o - -d "_" -f 0 '
    return command



rule all:
    input:
        expand(STR_COUNT_TEMPLATE, TREATMENT = TREATMENTS)
    
rule combine_count:
    input:
        BED_LIST = lambda w: expand(DEDUP_STR_BED, SAMPLE_NAME=select_sample(w))
    
    output:
        STR_COUNT_TEMPLATE

    shell:
        'cat {input.BED_LIST} '\
        '| cut -f1,6 '\
        '| sort -k1,1 -k6,6 '\
        '| uniq -c '\
        "| awk '{{ print $2,$3,$1 }}' OFS='\\t' "\
        '| sort -k3,3nr '\
        '> {output}'


rule dedup_bed:
    input:
        BED = STR_BED

    params:
        DEDUP_COMMAND = lambda w: dedup_command(w)

    output:
        BED = DEDUP_STR_BED
    
    shell:
        'cat {input.BED} | {params.DEDUP_COMMAND} > {output.BED}'

rule make_bed:
    input:
        BAM = STR_BAM
    
    output:
        BED = STR_BED
    
    shell:
        'samtools view -F4 -F2048 -bF256 {input.BAM} '\
        '| bam_to_bed.py -i -  -o - '\
        '| sort -k1,1 -k2,2n -k3,3n -k6,6 ' \
        '> {output.BED} '


rule mapping:
    input:
        FQ1 = FQ1, 
        FQ2 = FQ2,

    params:
        THREADS = THREADS,
        INDEX = REF_INDEX

    output:
        BAM = STR_BAM
    
    shell:
        'bowtie2 --local -p {params.THREADS} -x {params.INDEX} '\
        '--no-discordant --no-mixed --dovetail --mm '\
        ' -1 {input.FQ1} -2 {input.FQ2} '\
        '| samtools view -b@ {params.THREADS} - '\
        '> {output.BAM}'
   
