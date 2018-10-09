#!/usr/bin/env python

import pandas as pd 
import glob
import os
import sys
import re
from functools import partial
from multiprocessing import Pool
from collections import deque
import time

PROJECT_PATH='/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
COUNT_PATH = PROJECT_PATH + '/Counts/all_counts'
REF_BED_PATH = os.environ['REF'] + '/hg19/new_genes'

def deduplicate(bam_tag, bam):
    '''
    genearte dedup command
    '''
    if not re.search('genome-sim|L[12E]', bam):
        toleration = 1 if  bam_tag in ['all','repeats'] else 0
        DEDUP_COMMAND = "  deduplicate_bed.py -i - -d '_' -f 0 -t {TOLERATE} --ct 6" \
                        "| poisson_umi_adjustment.py -i - -o - --umi 6 "\
                        .format(TOLERATE=toleration)
    else:
        filter_command = "awk '$NF!~/N/ && ($3-$2) <= 300' |" if bam_tag == 'snc_all' else ''
        DEDUP_COMMAND = filter_command + ' sort -k1,1 -k2,2n -k3,3n -k6,6n -u ' \
                        '| cut -f1-6 '
    return DEDUP_COMMAND

def get_parameter(bam_tag):
    REF_BED = ''
    OUT_PATH = ''

    if bam_tag == 'all':
        REF_BED = REF_BED_PATH + '/genes.bed'
        OUT_PATH = COUNT_PATH + '/counts'
    elif bam_tag == 'snc_all':
        REF_BED = REF_BED_PATH + '/sncRNA_x_protein.bed'
        OUT_PATH = COUNT_PATH + '/sncRNA'
    elif bam_tag == 'sncRNA':
        REF_BED = REF_BED_PATH + '/smallRNA.bed'
        OUT_PATH = COUNT_PATH + '/smallRNA'
    elif bam_tag == 'rRNA_mt':
        REF_BED = REF_BED_PATH + '/rRNA_mt.bed'
        OUT_PATH = COUNT_PATH + '/rRNA_mt'
    elif bam_tag == 'repeats':
        REF_BED = os.environ['REF'] + '/hg19/genome/rmsk.bed.gz'
        OUT_PATH = COUNT_PATH + '/repeats'
    
    assert OUT_PATH != "" and REF_BED != "", bam_tag
    if not os.path.isdir(OUT_PATH):
        os.makedirs(OUT_PATH)
    return REF_BED, OUT_PATH


def get_date_created(f):
    s = time.time() - os.path.getmtime(f)
    hr = s/3600
    return hr


def strand_count(TEMP_BED, REF_BED, TEMP_FOLDER, OUT_PATH, samplename, dedup_label = 'all'):

    strand_count_command = ''

    field_chooser = '| cut -f7- ' if dedup_label == 'dedup' else '| cut -f8-'

    for strand in ['sense','antisense']:
        STRAND = '-s' if strand == 'sense' else '-S'
        TEMP_INTERSECTED = TEMP_BED.replace('.bed.gz', '.%s.intersected.bed.gz' %strand ) 
        OUT_COUNT_FILE = OUT_PATH + '/' + samplename + '.' +  dedup_label + '.' + strand + '.counts'

        if not os.path.isfile(OUT_COUNT_FILE) or get_date_created(OUT_COUNT_FILE) > 24:
            strand_count_command += "; bedtools intersect -a {TEMP_BED} -b {REF_BED} -wao {STRAND} "\
                                    "| bgzip | tee  {TEMP_INTERSECTED} | zcat " \
                                    " {field_chooser} " \
                                    "| sort --temporary-directory={TEMP_FOLDER} "\
                                    "| uniq -c " \
                                    .format(TEMP_BED = TEMP_BED,
                                            REF_BED = REF_BED,
                                            TEMP_INTERSECTED = TEMP_INTERSECTED,
                                            TEMP_FOLDER = TEMP_FOLDER,
                                            field_chooser = field_chooser,
                                            STRAND = STRAND) 
            strand_count_command += "| awk {'print $2,$3,$4,$5,$6,$7,$8,$9, $1 '} OFS='\\t'"\
                                    "> " + OUT_COUNT_FILE
        else:
            print('Using old file: %s' %OUT_COUNT_FILE)
    return strand_count_command



def process_bam(bam_tag, bam, samplename):
    TEMP_FOLDER = PROJECT_PATH + \
            '/' + samplename + \
            '.' + os.path.basename(os.path.dirname(bam)) + \
            '.' + os.path.basename(bam).strip(' ') + \
            '_TEMP'
    TEMP_FOLDER = TEMP_FOLDER.replace(' ','.')

    REF_BED, OUT_PATH = get_parameter(bam_tag)
    TEMP_BED = bam.replace('.bam','.bed.gz')
    DEDUP_BED = TEMP_BED.replace('.bed.gz','.dedup.bed.gz') 
    DEDUP_COMMAND = deduplicate(bam_tag, bam)

    STRAND_COUNT_COMMAND = strand_count(TEMP_BED, REF_BED, TEMP_FOLDER, 
                                        OUT_PATH, samplename, 
                                        dedup_label = 'all')
    STRAND_COUNT_COMMAND += strand_count(DEDUP_BED, REF_BED, TEMP_FOLDER, 
                                        OUT_PATH, samplename, 
                                        dedup_label = 'dedup')
    
    assert os.stat(bam).st_size>0 , bam
    BASE_COMMAND = " mkdir -p {TEMP_FOLDER} ;"\
            "cat {BAM} "\
            "| samtools view -bF 4 -F 256 -F 2048 "\
            "| bam_to_bed.py -i - --add_cigar --primary " \
            "| sort -k1,1 -k2,2n -k3,3n -k6,6 --temporary-directory={TEMP_FOLDER}"\
            "| bgzip "\
            "| tee {TEMP_BED} "\
            "| zcat | {DEDUP_COMMAND} "\
            "| bgzip "\
            "> {DEDUP_BED} "\
            "{STRAND_COUNT_COMMAND}; rm -rf {TEMP_FOLDER}"\
            .format(TEMP_FOLDER = TEMP_FOLDER,
                    BAM = bam,
                    TEMP_BED = TEMP_BED,
                    DEDUP_BED = DEDUP_BED,
                    DEDUP_COMMAND = DEDUP_COMMAND,
                    STRAND_COUNT_COMMAND = STRAND_COUNT_COMMAND)
    #print(BASE_COMMAND)
    os.system(BASE_COMMAND)
    print('Finished %s: %s' %(samplename, bam_tag))


def process_sample(count_path, sample_folder):
    all_aligned_path = sample_folder + '/Combined'
    small_RNA_path = sample_folder + '/smallRNA'
    rRNA_path = sample_folder + '/rRNA_mt'
    samplename = os.path.basename(sample_folder)

    BAMs = {'all':all_aligned_path + '/primary_no_sncRNA_tRNA_rRNA_repeats.bam',
            'snc_all':all_aligned_path + '/sncRNA.bam',
            'sncRNA':small_RNA_path + '/aligned.bam',
            'rRNA_mt': rRNA_path + '/aligned.bam',
            'repeats': all_aligned_path + '/repeats.bam'}
    
    for bam_tag, bam in BAMs.items():
        process_bam(bam_tag, bam, samplename)


def main():
    if not os.path.isdir(COUNT_PATH):
        os.makedirs(COUNT_PATH)
    
    sample_folders = glob.glob(PROJECT_PATH + '/*001')
    sample_folders = filter(lambda x: 'try' not in x, sample_folders)
    COUNT_FUNCTION = partial(process_sample, COUNT_PATH)

    debug = True
    debug= False
    if debug:
        deque(map(COUNT_FUNCTION, sample_folders))
    else:
        p = Pool(24)
        p.map(COUNT_FUNCTION, sample_folders)
        p.close()
        p.join()



if __name__ == '__main__':
    main()
