#!/bin/bash

PEAK_PATH=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/genome_peaks/merged_peak/adjusted_files
PEAK_FILES=$PEAK_PATH/*.bed
REGION_PATH=$PEAK_PATH/gene_loc
GENOME=$REF/hg19/genome/hg19_genome.genome
BED12=$REF/hg19/genome/genes.bed12
mkdir -p $REGION_PATH

for PEAK_FILE in $PEAK_FILES
do
	SAMPLENAME=$(basename ${PEAK_FILE%.bed})
	echo bedtools bedtobam \
			-g $GENOME \
			-i $PEAK_FILE \
		\| read_distribution.py -i /dev/stdin -r $BED12 \
		\> $REGION_PATH/${SAMPLENAME}.rseqc
done

