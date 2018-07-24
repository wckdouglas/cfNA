#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map
TRIM_PATH=$PROJECT_PATH/Trim
SALMON_PATH=$PROJECT_PATH/salmon_result
SALMON_REF_PROTEIN=$REF/hg19/new_genes/salmon_proteins_idx
SALMON_REF=$REF/hg19/new_genes/salmon_genes_ercc_idx
BAM_PATH=$SALMON_PATH/bam_files
GENE_MAP_PROTEIN=$REF/hg19/new_genes/proteins.gtf
GENE_MAP=$REF/hg19/new_genes/genes_ercc.gtf
mkdir -p $SALMON_PATH $BAM_PATH

for R1 in `ls $TRIM_PATH/*.1.fq.gz | egrep --color=no 'L[12E]|Frag'`
do
    R2=${R1/.1./.2.}
    SAMPLENAME=$(basename ${R1%_R1_001.1.fq.gz})

    if echo $SAMPLENAME | egrep -q 'L[12E]'
    then
        LIBTYPE=IU
    else
        LIBTYPE=ISF
    fi

    #echo salmon quant --mates1 $R1 --mates2 $R2 \
    echo salmon quant -r $R1 \
            --threads 12 --index $SALMON_REF -o $SALMON_PATH/$SAMPLENAME \
            --gcBias --seqBias --geneMap $GENE_MAP \
            --libType $LIBTYPE \
            --writeMappings \
        \> $BAM_PATH/${SAMPLENAME}.bam
done
