#!/usr/bin/env python

import os
import glob
from multiprocessing import Pool
import re
from functools import partial



def filter_protein_bam(in_bam, filter_path, label, refflat, threads = 6, return_command=True):

    name_sorted_bam = in_bam.replace('.bam','.name_sort.bam')
    filtered_bam = filter_path + '/' + label + '.bam'
    plus_bam = filter_path + '/' + label + '.plus.bam' 
    plus_anti_bam = filter_path + '/' + label + '.plus_anti.bam' 
    plus_sense_bam = filter_path + '/' + label + '.plus_sense.bam' 
    minus_bam = filter_path + '/' + label + '.minus.bam'
    minus_anti_bam = filter_path + '/' + label + '.minus_anti.bam'
    minus_sense_bam = filter_path + '/' + label + '.minus_sense.bam'
    protein_sense_bam = filter_path + '/' + label + '.protein_sense.bam'
    protein_anti_bam = filter_path + '/' + label + '.protein_anti.bam'


    snc_annotation = os.environ['REF'] + '/hg19/new_genes/sncRNA_rRNA_for_bam_filter.bed'
    rmsk_annotation = os.environ['REF'] + '/hg19/genome/rmsk.bed'
    protein_bed = os.environ['REF'] + '/hg19/new_genes/protein.bed'
    plus_bed = os.environ['REF'] + '/hg19/new_genes/protein_plus.bed'
    minus_bed = os.environ['REF'] + '/hg19/new_genes/protein_minus.bed'


    filtering_plus = ' awk \'{if (($1~/^@/)|| ($2==99)||($2==147)|| ($2==355)||($2==403)|| ($2==1123)||($2==1171)) print}\''
    filtering_minus = 'awk \'{if (($1~/^@/)|| ($2==83)||($2==163)|| ($2== 339)||($2==419)|| ($2==1187)||($2==1107)) print}\''

    command = ''

    command += '; sambamba sort -t {threads} -n -p -o /dev/stdout {in_bam}'\
            '| samtools fixmate -@ {threads} - {name_sorted}'\
            '; samtools view -h@ {threads} {name_sorted} | {filtering_plus} ' \
            '| samtools view -b@ {threads} - > {plus_bam} '\
            '; samtools view -h@ {threads} {name_sorted} | {filtering_minus} '\
            '| samtools view -b@ {threads} - > {minus_bam} '\
            '; bedtools pairtobed -abam {plus_bam} -b {sncRNA} -type neither  '\
                '| bedtools pairtobed -abam - -b {rmsk} -type neither '\
                '| bedtools pairtobed -abam - -b {plus_bed} > {plus_sense_bam} '\
            '; bedtools pairtobed -abam {minus_bam} -b {sncRNA} -type neither '\
                '| bedtools pairtobed -abam - -b {rmsk} -type neither '\
                '| bedtools pairtobed -abam - -b {minus_bed} > {minus_sense_bam} '\
            '; bedtools pairtobed -type neither -abam {minus_bam} -b {minus_bed} '\
            '| bedtools pairtobed -abam - -b {plus_bed} > {minus_anti_bam} '\
            '; bedtools pairtobed -type neither -abam {plus_bam} -b {plus_bed} '\
            '| bedtools pairtobed -abam - -b {minus_bed} > {plus_anti_bam} '\
            '; sambamba merge -t {threads} /dev/stdout {plus_sense_bam} {minus_sense_bam} '\
            '| samtools view -bF 1024 -@ {threads} '\
            '| sambamba sort -t {threads} -o {protein_sense_bam} /dev/stdin'\
            '; sambamba merge -t {threads} /dev/stdout  {plus_anti_bam} {minus_anti_bam} '\
            '| samtools view -bF 1024 -@ {threads} '\
            '| sambamba sort -t {threads} -o {protein_anti_bam} /dev/stdin '\
            .format(threads = threads,
                    in_bam = in_bam,
                    name_sorted = name_sorted_bam,
                    filtering_plus = filtering_plus,
                    filtering_minus = filtering_minus,
                    plus_bam = plus_bam,
                    minus_bam = minus_bam,
                    sncRNA = snc_annotation, 
                    rmsk = rmsk_annotation,
                    plus_bed = plus_bed,
                    minus_bed = minus_bed,
                    plus_sense_bam = plus_sense_bam,
                    minus_sense_bam = minus_sense_bam,
                    minus_anti_bam = minus_anti_bam,
                    plus_anti_bam = plus_anti_bam,
                    protein_sense_bam = protein_sense_bam,
                    protein_anti_bam = protein_anti_bam)

    command += '; sambamba sort -p -n -o /dev/stdout -t {threads} {in_bam} '\
            '| samtools fixmate -@ {threads} - - '\
            '| bedtools pairtobed -abam - -b {protein} -type both '\
            ' | bedtools pairtobed -abam - -b {sncRNA} -type neither '\
            '| samtools view -bF 1024 -F 256 -F 2048 ' \
            '| sambamba sort -p -t {threads} -o {filtered_bam}  /dev/stdin '\
            .format(in_bam = in_bam,
                    threads = threads,
                    sncRNA = snc_annotation,
                    protein = protein_bed,
                    filtered_bam = filtered_bam)

    command += '; '
    command += run_rna_seq_picard(refflat, filter_path, protein_sense_bam, label + '.sense', return_command=True)
    command += '; '
    command += run_rna_seq_picard(refflat, filter_path, protein_anti_bam, label + '.anti', return_command=True)
    command += '; '
    command += run_rna_seq_picard(refflat, filter_path, filtered_bam, label + '.filtered', return_command=True)


    if return_command:
        return command
    else:
        os.system(command)



def run_dedup(name_sorted_bam_file, sorted_bam, dedup_bam, samplename, dry=False):
    sort_command = 'cat %s ' %name_sorted_bam_file +\
        '| bam_umi_tag.py -i - -t RX ' \
        '| samtools view -bF 256 -F4 -F2048 '\
        '| samtools addreplacerg -r ID:{id_name} -r SM:{id_name}  - '\
        '| samtools view -b '\
        '| sambamba sort -n -o /dev/stdout /dev/stdin '\
        '| picard FixMateInformation ADD_MATE_CIGAR=true '\
        ' ASSUME_SORTED=true INPUT=/dev/stdin OUTPUT=/dev/stdout ' \
        '| sambamba sort -o {sorted_bam} /dev/stdin' \
        .format(sorted_bam = sorted_bam ,
                id_name = samplename.replace('_R1_001',''))
    
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

    command = sort_command + '; ' + dedup_command
    
    if dry:
        print(command)
    else:
        os.system(command)
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


def run_qc(refflat, outpath, outpath_dedup, sample_path):
    bam_file = sample_path + '/Combined/primary.bam'
    #bam_file = sample_path + '/Combined/primary_no_sncRNA_tRNA_rRNA_repeats.bam'

    samplename = os.path.basename(sample_path)
    sorted_bam_file = bam_file.replace('.bam','.sorted.bam')
    dedup_bam = bam_file.replace('.bam','.deduplicated.bam')

    run_dedup(bam_file, sorted_bam_file, dedup_bam, samplename, dry=False)
    filter_protein_bam(sorted_bam_file, outpath, samplename, refflat)
    filter_protein_bam(dedup_bam, outpath_dedup, samplename, refflat, return_command=True)
    run_rna_seq_picard(refflat, outpath, sorted_bam_file, samplename)
    run_rna_seq_picard(refflat, outpath, dedup_bam, samplename + '_deduplicated')


def main():
    refflat = '/stor/work/Lambowitz/ref/hg19/new_genes/proteins.refflat'
    project_path = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
    outpath_dedup = project_path + '/picard_qc'
    outpath = project_path + '/picard_qc_all'
    if not os.path.isdir(outpath):
        os.mkdir(outpath)
    if not os.path.isdir(outpath_dedup):
        os.mkdir(outpath_dedup)

    samples = glob.glob(project_path + '/Q*001')
    samples = filter(lambda x: not re.search('L[12E]',x), samples)
    picard_func = partial(run_qc, refflat, outpath, outpath_dedup)
    p = Pool(24)
    p.map(picard_func, samples)
    p.close()
    p.join()


if __name__ == '__main__':
    main()



