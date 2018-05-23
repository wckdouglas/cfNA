#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map
BAM_PATH=$PROJECT_PATH/merged_bam
COVERAGE_PATH=$BAM_PATH/coverage
GENOME=$REF/hg19/genome/hg19_size.tsv
THREADS=12
mkdir -p $COVERAGE_PATH 

for BAM in $BAM_PATH/*.bam
do
	SAMPLENAME=$(basename ${BAM%.bam})
	for STRAND in forward reverse
	do
		if [[ $STRAND == "forward" ]]
		## see http://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html
		then
			STRAND_ANNOT=reverse 
		elif [[ $STRAND == "reverse" ]]
		then
			STRAND_ANNOT=forward
		fi
		OUT_PREFIX=${SAMPLENAME}.${STRAND_ANNOT}
		BW=$COVERAGE_PATH/${OUT_PREFIX}.bigWig
		echo sambamba index -t $THREADS $BAM \
			\;bamCoverage \
			--numberOfProcessors $THREADS \
			-b $BAM \
			-o $BW \
			--filterRNAstrand $STRAND \
			--outFileFormat bigwig \
			--samFlagExclude 1024 \
			--binSize 1
	done
done | grep -v name_sort
