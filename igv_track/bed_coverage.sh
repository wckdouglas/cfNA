#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map
BED_PATH=$PROJECT_PATH/merged_bed
COVERAGE_PATH=$BED_PATH/coverage
GENOME=$REF/hg19/genome/hg19_genome.fa.fai
mkdir -p $COVERAGE_PATH

for BED in $BED_PATH/stranded/*bed.gz $BED_PATH/alkaline.bed.gz
do
    SAMPLENAME=$(basename ${BED%.bed.gz})
    TEMP=$COVERAGE_PATH/$SAMPLENAME
    OUT_BW=${TEMP}.bigWig
    echo bedtools genomecov \
            -i $BED \
            -g $GENOME \
            -bga \
        \| sort -k1,1 -k2,2n \
        \> $TEMP \
        \; bedGraphToBigWig $TEMP $GENOME $OUT_BW \
        \; rm $TEMP
done
