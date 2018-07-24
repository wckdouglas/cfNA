#!/bin/bash

BED_PATH=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/merged_bed
OUT_PATH=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/classifier
TEST_BED=$OUT_PATH/test.bed
TRAIN_BED=$OUT_PATH/train.bed
TEMP=$OUT_PATH/temp

# sample DNA
SAMPLE_SIZE=10000000
##from shendure
#zcat /stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/bed_files/SRR2130051.bed.gz \
#    | awk '($3-$2) < 100 {print $0,"DNA"}' OFS='\t'  \
#    | shuf -n $((10*$SAMPLE_SIZE)) \
#    > $TRAIN_BED
#
zcat /stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/bed_files/SRR2130051.bed.gz \
    | awk '($3-$2) > 100 && ($3-$2) < 400 {print $0,"DNA"}' OFS='\t'  \
    | shuf -n $SAMPLE_SIZE \
    >> $TRAIN_BED
#
zcat $BED_PATH/alkaline.no_sncRNA.bed.gz \
    | awk '($3-$2) > 100 && ($3-$2) < 400 {print $0,"DNA"}' OFS='\t'  \
    | shuf -n $SAMPLE_SIZE \
    >> $TRAIN_BED

zcat $BED_PATH/alkaline.no_sncRNA.bed.gz \
    | awk '($3-$2) < 100 {print $0,"DNA"}' OFS='\t'  \
    | shuf -n $((10*$SAMPLE_SIZE)) \
    >> $TRAIN_BED


# sample RNA
#    | bedtools intersect -a - -b $REF/hg19/genome/sncRNA_x_protein.sorted.bed.gz -v \
zcat $BED_PATH/unfragmented.bed.gz \
    | awk '($3 - $2) < 100  {print $0,"RNA"}' OFS='\t' \
    >> $TRAIN_BED
#
#
TOTAL=$(cat $TEMP | wc -l)
TEST_SAMPLE=$(expr $TOTAL / 50)
cat $TRAIN_BED | shuf > $TEMP
python validation_bed.py $TEMP $OUT_PATH 100000
#rm $TEMP
