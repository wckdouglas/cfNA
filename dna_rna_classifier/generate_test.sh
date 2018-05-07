#!/bin/bash

BED_PATH=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/bed_files
TEST_BED=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/classifier/test.bed

# sample DNA
SAMPLE=500000
for i in Na1 NA2
do
    zcat $BED_PATH/Qcf_${i}_R1_001.bed.gz \
        | awk '{print $0,"DNA"}' OFS='\t'  
done | shuf -n $SAMPLE > $TEST_BED 


# sample RNA
for i in 7 8 9 10
do
    zcat $BED_PATH/Qcf${i}_R1_001.bed.gz \
        | awk '($3 - $2) < 80  {print $0,"RNA"}' OFS='\t' 
done | shuf -n $SAMPLE >> $TEST_BED
