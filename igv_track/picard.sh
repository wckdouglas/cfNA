#!/bin/bash

MERGED_FILTER_BAM_PATH=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/merged_bam/filtered_bam
REFFLAT=/stor/work/Lambowitz/ref/hg19/genome/refFlat.txt

for BAM in $MERGED_FILTER_BAM_PATH/*bam
do
    echo picard CollectRnaSeqMetrics \
        I=$BAM \
        O=${BAM%.bam}.RNA_Metrics \
        REF_FLAT=$REFFLAT \
        STRAND=FIRST_READ_TRANSCRIPTION_STRAND AS=false 
done

