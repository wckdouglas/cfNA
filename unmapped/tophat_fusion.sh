#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map
REF_INDEX=$REF/hg19/genome/hg19_genome
MICROBE=$REF/microbiome/all_seqs

for SAMPLE_FOLDER in $PROJECT_PATH/*001
do
    OUT_DIR=$SAMPLE_FOLDER/unmapped
    mkdir -p $OUT_DIR
    SAMPLENAME=$(basename $SAMPLE_FOLDER)
    BAM_FILE=$SAMPLE_FOLDER/Bowtie/bowtie2.bam

    UNMAP_FQ1=$OUT_DIR/unmapped.1.fq.gz
    UNMAP_FQ2=$OUT_DIR/unmapped.2.fq.gz

    echo samtools fastq -f4 $BAM_FILE \
        \| bowtie2 -x $REF/Ecoli/k12_mg1655 --no-discordant --no-mixed  \
            --very-sensitive-local --score-min G,1,10 --mm \
           --interleaved - \
        \| samtools fastq -f4 - \
        \| deinterleave_fastq.py -1 $UNMAP_FQ1  -2 $UNMAP_FQ2 -i - \
        \; tophat2 --fusion-search  \
            --b2-very-sensitive \
            --no-coverage-search \
            --fusion-anchor-length 8 \
            --output-dir $OUT_DIR \
            $REF_INDEX \
            $UNMAP_FQ1 $UNMAP_FQ2 
done

