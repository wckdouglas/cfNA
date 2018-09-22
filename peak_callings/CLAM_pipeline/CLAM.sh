#!/usr/bin/env python

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map

for SAMPLE_FOLDER in $PROJECT_PATH/*001
do
    BOWTIE2_BAM=$SAMPLE_FOLDER/Bowtie/bowtie2.bam
    HISAT2_BAM=$SAMPLE_FOLDER/Hisat/hisat.bam
    COMBINED_PATH=$SAMPLE_FOLDER/Combined
    COMBINED_BAM=$COMBINED_PATH/merged.bam
    CLAM_PATH=$COMBINED_PATH/CLAM
    echo cat $BOWTIE2_BAM \
            \| python tag_bam.py -i - -o - -t NH  \
            \| samtools cat $HISAT2_BAM - \
            \> $COMBINED_BAM \
        \; CLAM preprocessor -i $COMBINED_BAM -o $CLAM_PATH --read-tagger-method median \
        \; CLAM realigner -i $COMBINED_BAM -o $CLAM_PATH  \
        \; cat $CLAM_PATH/unique.sorted.bam \
            \| python tag_bam.py -i - -o - -t AS --tag_value 1 \
            \| sambamba sort -n -o $CLAM_PATH/unique.name_sorted.bam /dev/stdin \
        \; sambamba sort -n $CLAM_PATH/realigned.sorted.bam -o $CLAM_PATH/realigned.name_sorted.bam \
        \; samtools cat $CLAM_PATH/realigned.name_sorted.bam $CLAM_PATH/unique.name_sorted.bam \
            \| bam_to_bed.py -i - -t AS \
            \| sort -k1,1 -k2,2n -k3,3n -k6,6 \
            \| python dedup.py -i - -o - \
            \| bgzip \
            \> $CLAM_PATH/pseudo_fragment.bed.gz \
            \; tabix -f $CLAM_PATH/pseudo_fragment.bed.gz 
done | egrep  '[Nn][Aa][0-9]+|[aA]ll[0-9]+|[Qq][cC][fF][0-9]+|[DE][DE]|Exo'
