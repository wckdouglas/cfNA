#!/bin/bash

BAM_PATH=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/merged_bam/filtered_bam

for BAM in $BAM_PATH/*bam
do
    for EXPR in High Low Medium
    do
        echo sambamba index ${BAM}\
            \;geneBody_coverage.py \
            -i $BAM \
            -r $REF/hg19/new_genes/bone_marrow.${EXPR}.bed \
            -o ${BAM%.bam}.$EXPR
    done
done
