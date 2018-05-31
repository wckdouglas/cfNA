#!/bin/bash

BED_PATH=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/merged_bed
STRANDED_BED_PATH=$BED_PATH/stranded
OUT_PATH=$BED_PATH/MACS2
mkdir -p $OUT_PATH

for SAMPLE in unfragmented all
do
    for STRAND in fwd rvs
    do
        echo macs2 callpeak \
            --treatment $STRANDED_BED_PATH/${SAMPLE}.${STRAND}.bed.gz  \
            --control $BED_PATH/alkaline.bed.gz \
            --outdir $OUT_PATH \
            --name ${SAMPLE}.${STRAND} \
            --nomodel \
            --format BEDPE \
            --keep-dup all \
            --gsize hs --broad
    done
done
