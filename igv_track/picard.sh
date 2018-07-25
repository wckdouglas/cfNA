#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map
REFFLAT=/stor/work/Lambowitz/ref/hg19/new_genes/proteins.refflat

for SAMPLE_FOLDER in $PROJECT_PATH/Q*L*001
do
    BAM=$SAMPLE_FOLDER/Combined/primary_no_sncRNA_tRNA_rRNA.bam
    SORTED_BAM=${BAM/.bam/.sorted.bam}
    echo sambamba sort -p  -o $SORTED_BAM $BAM \
        \; picard CollectRnaSeqMetrics \
        I=$SORTED_BAM \
        O=${BAM%.bam}.RNA_Metrics \
        REF_FLAT=$REFFLAT \
        STRAND=NONE AS=false 
done

