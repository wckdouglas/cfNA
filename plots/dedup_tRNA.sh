#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map

for FOLDER in $PROJECT_PATH/*001
do
    echo cat $FOLDER/rRNA_tRNA_premap/tRNA_rRNA.bam \
        \| bam_to_bed.py -i - --add_cigar \
        \| sort -k1,1 -k2,2n -k3,3n -T $FOLDER \
        \| deduplicate_bed.py -i - -o - -d'_'  --ct 6 \
        \| poisson_umi_adjustment.py -i - -o - --umi 6 \
        \> $FOLDER/rRNA_tRNA_premap/tRNA_rRNA.dedup.bed
done
