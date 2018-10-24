#!/bin/bash
PROJECT_PATH=$SCRATCH/cfNA
PROJECT_PATH=$WORK/cdw2854/cfNA


### local cloud
###
#PROJECT_PATH=$WORK/cdw2854/cell_Free_nucleotides


DATA_PATH=$PROJECT_PATH/data
RESULT_PATH=$PROJECT_PATH/tgirt_map
#RESULT_PATH=$PROJECT_PATH/tgirt_map_new_penalty
REF_PATH=$REF/hg19/genome
NEW_GENE_PATH=$REF/hg19/new_genes
LOG_PATH=$RESULT_PATH/log
UNIVEC=$REF/UniVec/UniVec_Core
THREADS=24




mkdir -p $LOG_PATH

for FQ1 in $DATA_PATH/*_R1_001.fastq.gz
do
	FQ2=${FQ1/_R1_/_R2_}
	SAMPLENAME=$(basename ${FQ1%_R1_001.fastq.gz})
    polyA=''
	if echo $SAMPLENAME | grep -q '450\|200\|^[INOQS]'
	then
		TTN="--TTN"
		UMI="--umi 6 --count_all"

#		UMI="--umi 6 "
    elif echo $SAMPLENAME | grep -q 'TEV[12]'
    then
        TTN="--TTN"
        UMI=" "
        polyA='--polyA'
	else
		TTN=' '
		UMI=' '
	fi


    if echo $SAMPLENAME | egrep -q 'L[12E]|TEV3'
    then
        TTN=' '
        UMI=' '
        polyA='--polyA'
    fi

	echo tgirt_count.py map -1 $FQ1 -2 $FQ2 \
		--outdir $RESULT_PATH \
        --samplename ${SAMPLENAME}_R1_001 \
        --univec $UNIVEC \
		--hisat_index $REF_PATH/hg19_genome \
		--bowtie2_index $REF_PATH/hg19_genome \
		--bedpath $NEW_GENE_PATH \
		--splicesite $NEW_GENE_PATH/splicesites.tsv \
		--rRNA_mt_index $NEW_GENE_PATH/rRNA_mt \
        --smRNA_index $NEW_GENE_PATH/smallRNA \
		-p $THREADS $UMI $TTN \
        --trim_aggressive ${polyA} \
		--repeats $REF_PATH/rmsk.bed.gz \
		--repeats_index $REF_PATH/repeats/all_rmsk_From_bed \
        --snakemake \
		2\>\&1 \
		\| tee $RESULT_PATH/log/${SAMPLENAME}.log
done |  egrep -v  'TEV|TeI|GsI|SRR|[TG]0|200|450|[NO][QN]' #| egrep 'IGG|200|OQ|NN|NQ|QCF|S96|ON'
#		--skip_trim  --skip_hisat --skip_premap --skip_bowtie --skip_post_process_bam --skip_remap \
#		--repeats_index $REF_PATH/rmsk \
#		--repeats_index $REF_PATH/repeat_mask/all_rmsk_From_bed \
#		--repeats_index $REF_PATH/hs36.combined.fa \
