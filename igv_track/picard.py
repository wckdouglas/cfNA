#!/usr/bin/env python

import os
import glob
from multiprocessing import Pool
import re
from functools import partial

def run_dedup(name_sorted_bam_file, sorted_bam, dedup_bam):
    sort_command = 'cat %s ' %name_sorted_bam_file +\
        '| bam_umi_tag.py -i - -t RX ' \
        '| samtools view -bF 256 -F4 -F2048 '\
        '| sambamba sort -n -o /dev/stdout /dev/stdin '\
        '| picard FixMateInformation ADD_MATE_CIGAR=true '\
        ' ASSUME_SORTED=true INPUT=/dev/stdin OUTPUT=/dev/stdout ' \
        '| sambamba sort -o %s /dev/stdin' %sorted_bam 
    
    metric_file = dedup_bam.replace('.bam','.dedup_metric')
    umi_metric_file = dedup_bam.replace('.bam','.umi_metric')
    mark_dup_bam = dedup_bam.replace('.deduplicated.','.marked_duplicate.')
    dedup_command = 'picard UmiAwareMarkDuplicatesWithMateCigar UMI_METRICS_FILE={UMI_METRIC} '\
                        .format(UMI_METRIC=umi_metric_file)+\
                    ' MAX_EDIT_DISTANCE_TO_JOIN=1 TAG_DUPLICATE_SET_MEMBERS=true '+\
                    ' UMI_TAG_NAME=RX INPUT={SORTED_BAM} OUTPUT=/dev/stdout REMOVE_SEQUENCING_DUPLICATES=true'\
                        .format(SORTED_BAM = sorted_bam) +\
		    ' METRICS_FILE={METRIC} REMOVE_DUPLICATES=false ASSUME_SORT_ORDER=coordinate '\
                        .format(METRIC=metric_file) +\
                    '| tee {MARKDUP_BAM}'.format(MARKDUP_BAM = mark_dup_bam) +\
                    '| samtools view -bF 1024 ' +\
                    '> {DEDUP_BAM}'.format(DEDUP_BAM = dedup_bam) +\
                    '; samtools index {DEDUP_BAM}'.format(DEDUP_BAM=dedup_bam)
    
    os.system(sort_command)
    os.system(dedup_command)
    #print(sort_command)
    #print(dedup_command)


def run_rna_seq_picard(refflat, outpath, bam_file, samplename, return_command = False):
    command = 'picard CollectRnaSeqMetrics ' +\
            'I=%s ' %(bam_file) +\
            'O=%s/%s.RNA_Metrics ' %(outpath, samplename) +\
            'REF_FLAT=%s ' %(refflat) +\
            'STRAND=FIRST_READ_TRANSCRIPTION_STRAND AS=false'
    if return_command:
        return command

    else:
        os.system(command)


def run_qc(refflat, outpath, sample_path):
    bam_file = sample_path + '/Combined/primary.bam'
    #bam_file = sample_path + '/Combined/primary_no_sncRNA_tRNA_rRNA_repeats.bam'

    samplename = os.path.basename(sample_path)
    sorted_bam_file = bam_file.replace('.bam','.sorted.bam')
    dedup_bam = bam_file.replace('.bam','.deduplicated.bam')

    run_dedup(bam_file, sorted_bam_file, dedup_bam)
    run_rna_seq_picard(refflat, outpath, sorted_bam_file, samplename)
    run_rna_seq_picard(refflat, outpath, dedup_bam, samplename + '_deduplicated')


def main():
    refflat = '/stor/work/Lambowitz/ref/hg19/genome/refFlat.txt'
    project_path = '/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map'
    outpath = project_path + '/picard_qc'
    if not os.path.isdir(outpath):
        os.mkdir(outpath)

    samples = glob.glob(project_path + '/Q*001')
    samples = filter(lambda x: not re.search('L[12]',x), samples)
    picard_func = partial(run_qc, refflat, outpath)
    p = Pool(24)
    p.map(picard_func, samples)
    p.close()
    p.join()


if __name__ == '__main__':
    main()



