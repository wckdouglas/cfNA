#!/bin/bash

BED_PATH=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/merged_bed
OUT_PATH=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/classifier
TEST_BED=$OUT_PATH/test.bed
TRAIN_BED=$OUT_PATH/train.bed
TEMP=$OUT_PATH/temp

# sample DNA
SAMPLE=10000000
HALF_SAMPLE=$(( $SAMPLE/2 ))
TEST_SAMPLE=$(( $SAMPLE/10 ))

zcat $BED_PATH/alkaline.no_sncRNA.bed.gz \
    | awk '($3-$2) > 100 && ($3-$2) < 200 {print $0,"DNA"}' OFS='\t'  \
    | shuf -n $HALF_SAMPLE \
    > $TRAIN_BED

zcat $BED_PATH/alkaline.no_sncRNA.bed.gz \
    | awk '($3-$2) < 100 {print $0,"DNA"}' OFS='\t'  \
    | shuf -n $HALF_SAMPLE \
    >> $TRAIN_BED 


# sample RNA
zcat $BED_PATH/unfragmented.bed.gz \
    | bedtools intersect -a - -b $REF/hg19/genome/sncRNA_x_protein.sorted.bed.gz -v \
    | awk '($3 - $2) < 100  {print $0,"RNA"}' OFS='\t' \
    | shuf -n $SAMPLE \
    >> $TRAIN_BED


cat $TRAIN_BED | shuf > $TEMP
head -n $TEST_SAMPLE $TEMP > ${TEST_BED}
tail -n +$TEST_SAMPLE $TEMP > ${TRAIN_BED}
rm $TEMP
