#!/bin/bash

PEAK_PATH=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/merged_bed/genome_WPS/genome_peaks
ANNOTATED_PATH=$PEAK_PATH/annotated
ANNOTATION_FILE=$REF/hg19/genome/all_annotation.bed.gz
mkdir -p $ANNOTATED_PATH

for PEAK_FILE in `ls $PEAK_PATH | grep -v 'no' | cut -d '.' -f1  | sort | uniq`
do
    echo cat $PEAK_PATH/${PEAK_FILE}.rvs.bed \
        $PEAK_PATH/${PEAK_FILE}.fwd.bed  \
    \| awk \''$5 != 0 && $8 > 3'\' \
    \| sort -k1,1 -k2,2n -k3,3n \
    \| bedtools intersect -wao -a - \
        -b $ANNOTATION_FILE \
    \| datamash -g 1,2,3,4,5,6,7,8 collapse 12,14,15 \
    \| sort -k5,5nr \
    \> $ANNOTATED_PATH/${PEAK_FILE}.annotated_peaks.tsv
done
