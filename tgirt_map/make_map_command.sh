#!/bin/bash
#PROJECT_PATH=$SCRATCH/cfNA
PROJECT_PATH=$WORK/cdw2854/cfNA


DATA_PATH=$PROJECT_PATH/data
RESULT_PATH=$PROJECT_PATH/tgirt_map
#RESULT_PATH=$PROJECT_PATH/tgirt_map_new_penalty
REF_PATH=$REF/hg19_ref
GENE_PATH=$REF_PATH/genes
GENOME_PATH=$REF_PATH/genome
LOG_PATH=$RESULT_PATH/log
THREADS=6


mkdir -p $LOG_PATH

ANALYSIS=No
ANALYSIS=ECOLI
ANALYSIS=PNAS
if [[ $ANALYSIS == "PNAS" ]]; then
    DATA_PATH=/stor/work/Lambowitz/cdw2854/cfNA/data/exosome
elif [[ $ANALYSIS == "ECOLI" ]]; then
    DATA_PATH=/stor/work/Lambowitz/Data/NGS/JA19062
fi


for FQ1 in $DATA_PATH/*_R1_001.fastq.gz
do
	FQ2=${FQ1/_R1_/_R2_}
	SAMPLENAME=$(basename ${FQ1%_R1_001.fastq.gz})
    polyA=''
	if echo $SAMPLENAME | egrep -q '450|200|^[INOQS]'
	then
		TTN="--TTN"
		UMI="--umi 6 --count_all"

#		UMI="--umi 6 "
    elif echo $SAMPLENAME | egrep -q '[PT]EV[0-9]+|^DC|^PC|^[0-9]+hr|F3'
    then
        TTN="--TTN"
        UMI=" "
    elif echo $SAMPLENAME | egrep -q '^EVD|^GC'
    then
        TTN=" "
        UMI=" "
    fi

    if echo $SAMPLENAME | egrep -q 'L[0-9E]+|TEV3'
    then
        TTN=' '
        UMI=' '
        polyA='--polyA'
    fi

	echo source activate tgirt_env \
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
