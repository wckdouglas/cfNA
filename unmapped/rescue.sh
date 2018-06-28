#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map
REF_INDEX=$REF/hg19/genome/hg19_genome
MICROBE=$REF/microbiome/all_seqs

bwa shm $REF_INDEX
for SAMPLE_FOLDER in $PROJECT_PATH/*001
do
    OUT_DIR=$SAMPLE_FOLDER/unmapped
    mkdir -p $OUT_DIR
    SAMPLENAME=$(basename $SAMPLE_FOLDER)
    BAM_FILE=$SAMPLE_FOLDER/Bowtie/bowtie2.bam
    echo seqtk mergepe $OUT_DIR/unmapped.1.fq.gz $OUT_DIR/unmapped.2.fq.gz \
        \| bwa mem -k15 -p $REF_INDEX  - \
        \| samtools view -b - \
        \| tee $OUT_DIR/rescued.bam \
        \| samtools view \
        \| awk \'' $2~/^77$|^141$/ {printf "@%s\n%s\n+\n%s\n", $1,$10,$11}'\' \
        \| bowtie2 --local --score-min G,1,10 -D 20 -R 3 -N 0 -L 8 -i S,1,0.50 \
            -x $MICROBE --interleaved - --mm \
        \| samtools view -b - \
        \> $OUT_DIR/microbe.bam
done

