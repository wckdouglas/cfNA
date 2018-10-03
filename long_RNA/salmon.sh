#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map
SALMON_PATH=$PROJECT_PATH/salmon_result
SALMON_REF_PROTEIN=$REF/hg19/new_genes/salmon_proteins_idx
#SALMON_REF=$REF/hg19/new_genes/salmon_genes_ercc_idx
BAM_PATH=$SALMON_PATH/bam_files
GENE_MAP_PROTEIN=$REF/hg19/new_genes/proteins.gtf
GENE_MAP=$REF/hg19/new_genes/genes_ercc.gtf
THREADS=12
mkdir -p $SALMON_PATH $BAM_PATH

for SAMPLENAME in `ls $PROJECT_PATH | egrep --color=no '001$' | egrep --color=no 'L[12E]|Frag'`
do
    SAMPLE_FOLDER=$PROJECT_PATH/$SAMPLENAME
    BAM=$SAMPLE_FOLDER/Combined/primary_no_sncRNA_tRNA_rRNA.bam
    R1=$SAMPLE_FOLDER/Combined/primary_no_sncRNA_tRNA_rRNA.1.fq
    R2=$SAMPLE_FOLDER/Combined/primary_no_sncRNA_tRNA_rRNA.2.fq

    if echo $SAMPLENAME | egrep -q 'L[12E]'
    then
        LIBTYPE=IU
    else
        LIBTYPE=ISF
    fi


    echo samtools fastq -@ $THREADS -1 $R1 -2 $R2 $BAM \
        \; salmon quant --mates1 $R1 --mates2 $R2 \
            --threads $THREADS --index $SALMON_REF_PROTEIN -o $SALMON_PATH/$SAMPLENAME \
            --gcBias --seqBias --geneMap $GENE_MAP_PROTEIN \
            --libType $LIBTYPE \
            --writeMappings \
        \> $BAM_PATH/${SAMPLENAME}.bam
done
