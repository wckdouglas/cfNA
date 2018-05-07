#!/bin/bash

PROJECTPATH=$WORK/cdw2854/cell_Free_nucleotides
PROJECTPATH=$SCRATCH/cell_Free_nucleotides/tgirt_map
BED_PATH=$PROJECTPATH/merged_bed
RESULT_PATH=$BED_PATH/genome_WPS
REF_PATH=$REF/hg19/genome
GENOME=$REF_PATH/hg19_genome.fa
PROGRAM=strandedGenomeWPS.py
#PROGRAM=strandedGenomeWPS.py
PYTHON=$(which python)

mkdir -p $RESULT_PATH
for BED in  $BED_PATH/*.bed.gz
do
    SAMPLENAME=$(basename ${BED%.bed.gz})
    echo $PYTHON -u $PROGRAM --inFile=$BED \
            --outprefix=$RESULT_PATH/${SAMPLENAME} \
            --genome=$GENOME \
            --window=50000  \
            --threads=4
done
