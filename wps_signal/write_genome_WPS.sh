#!/bin/bash

PROJECTPATH=$WORK/cdw2854/cell_Free_nucleotides
PROJECTPATH=$SCRATCH/cell_Free_nucleotides/tgirt_map
BED_PATH=$PROJECTPATH/bed_files
RESULT_PATH=$BED_PATH/genome_WPS
REF_PATH=$REF/hg19/genome
GENOME=$REF_PATH/hg19_genome.fa.fai
PROGRAM=strandedGenomeWPS.py
#PROGRAM=strandedGenomeWPS.py
PYTHON=$(which python)

mkdir -p $RESULT_PATH
for BED in  $BED_PATH/*.bed.gz
do
	for CHROM in $(seq 1 22) X Y
	do
		SAMPLENAME=$(basename ${BED%.bed.gz})
		echo $PYTHON -u $PROGRAM --inFile=$BED \
				--outprefix=$RESULT_PATH/${SAMPLENAME} \
				--genome=$GENOME \
				--window=50000 \
				--chromosome=chr$CHROM
	done
done | egrep 'IGG|200_|OQ|NN|NQ|QCF|S96|ON'  
