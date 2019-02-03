#!/bin/bash
PROJECT_PATH=$SCRATCH/cfNA
#PROJECT_PATH=$WORK/cdw2854/cfNA


DATA_PATH=$PROJECT_PATH/data
RESULT_PATH=$PROJECT_PATH/tgirt_map
#RESULT_PATH=$PROJECT_PATH/tgirt_map_new_penalty
REF_PATH=/scratch/02727/cdw2854/ref/hg19_ref
GENE_PATH=$REF_PATH/genes
GENOME_PATH=$REF_PATH/genome
LOG_PATH=$RESULT_PATH/log
THREADS=12


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
    elif echo $SAMPLENAME | grep -q '[PT]EV[0-9]+'
    then
        TTN="--TTN"
        UMI=" "
    fi

    if echo $SAMPLENAME | egrep -q 'L[0-9E]+|TEV3'
    then
        TTN=' '
        UMI=' '
        polyA='--polyA'
    fi

	echo source activate tgirt_map \
        \; tgirt_count.py map -1 $FQ1 -2 $FQ2 \
		--outdir $RESULT_PATH \
        --samplename ${SAMPLENAME}_R1_001 \
        --univec $GENE_PATH/UniVec_Core_Ecoli \
		--hisat_index $GENOME_PATH/hg19_genome \
		--bowtie2_index $GENOME_PATH/hg19_genome \
		--bedpath $GENE_PATH \
		--splicesite $GENE_PATH/splicesites.tsv \
		--rRNA_mt_index $GENE_PATH/rRNA_mt \
        --smRNA_index $GENE_PATH/smallRNA \
        --multi 20 -p $THREADS $UMI $TTN \
        --trim_aggressive ${polyA} \
		--repeats $GENE_PATH/rmsk.bed.gz \
		--repeats_index $GENE_PATH/rmsk \
        --rerun \
		2\>\&1 \
		\| tee $RESULT_PATH/log/${SAMPLENAME}.log
done |  egrep -v  '[PT]TeI|GsI|SRR|[TG]0|200|450|[NO][QN]' #| egrep 'IGG|200|OQ|NN|NQ|QCF|S96|ON'
#		--skip_trim  --skip_hisat --skip_premap --skip_bowtie --skip_post_process_bam --skip_remap \
#		--repeats_index $REF_PATH/rmsk \
#		--repeats_index $REF_PATH/repeat_mask/all_rmsk_From_bed \
#		--repeats_index $REF_PATH/hs36.combined.fa \
