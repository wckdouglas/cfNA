#!/bin/bash

PROJECT_PATH=/scratch/02727/cdw2854/cell_Free_nucleotides/tgirt_map
BED_PATH=$PROJECT_PATH/bed_files
WHITELIST=$REF/hg19/genome/wgEncodeDacMapabilityConsensusExcludable.bed.gz
mkdir -p $BED_PATH

for SAMPLE in `ls $PROJECT_PATH | egrep '001$' `
do
	SAMPLENAME=$(basename $SAMPLE)
	for RNA in sncRNA all repeats
	do
		if [[ $RNA == "sncRNA" ]]
		then
			IN_NAME=$SAMPLE/Combined/primary_no_sncRNA_tRNA_rRNA.bam
			OUT_NAME=$BED_PATH/${SAMPLENAME}_no_sncRNA.bed.gz
		elif [[ $RNA == "all" ]]
		then
			IN_NAME=$SAMPLE/Combined/primary.bam
			OUT_NAME=$BED_PATH/${SAMPLENAME}.bed.gz
		elif [[ $RNA == "repeats" ]]
		then
			IN_NAME=$SAMPLE/Combined/primary_no_sncRNA_tRNA_rRNA_repeats.bam
			OUT_NAME=$BED_PATH/${SAMPLENAME}_no_sncRNA_repeats.bed.gz
		fi

		TMP_FOLDER=${OUT_NAME/.bed.gz/_TMP}
		if  echo $SAMPLENAME |  egrep -q 'L[12]'
		then
			DEMUL=" "
		else
			DEMUL="| deduplicate_bed.py --infile - --outfile - --threshold 1 -d '_' "
            DEMUL=" $DEMUL | poisson_umi_adjustment.py -i - -o - --umi 6 "
		fi

		echo mkdir -p $TMP_FOLDER \
			\; cat $PROJECT_PATH/$IN_NAME \
			\| bam_to_bed.py -i - \
			\| sort -k1,1 -k2,2n -k3,3n -k6,6 --temporary-directory=$TMP_FOLDER \
			\| bedtools intersect -v -a - \
				-b $WHITELIST \
            $DEMUL \
			\| sort -k1,1 -k2,2n -k3,3n --temporary-directory=$TMP_FOLDER \
			\| bgzip \
			\> $OUT_NAME \
			\; tabix -f $OUT_NAME \
			\; rm -rf $TMP_FOLDER
	done
done
