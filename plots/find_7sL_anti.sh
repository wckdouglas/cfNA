#!/bin/bash

PROJECT_DIR=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map
DATA_DIR=/stor/work/Lambowitz/Data/NGS/JA18152

for SAMPLENAME in $(ls $PROJECT_DIR | egrep 'Q[cC][fF][0-9]+.*001')
do
    SAMPLE_FOLDER=$PROJECT_DIR/$SAMPLENAME
    BAM_FILE=$SAMPLE_FOLDER/Combined/primary.sorted.bam
    FQ=$DATA_DIR/${SAMPLENAME}.fastq.gz
    ID_FILE=$SAMPLE_FOLDER/7sl_anti.id
    echo zcat $REF/hg19/new_genes/all_annotation.bed.gz \
        \| awk \''$4~/7SL2$|^U4$|7SK|/'\' \
        \| bedtools intersect -b - -a $BAM_FILE -S -f 0.8 \
        \| samtools view -f64 \
        \| cut -f1 \
        \>  $ID_FILE \
        \; zgrep --color=no -A3 -f $ID_FILE $FQ \
        \| egrep  -v '^--$' \
        \> $SAMPLE_FOLDER/7sl_R1.fq
done
