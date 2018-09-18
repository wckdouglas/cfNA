#!/bin/bash


PROJECT_PATH=/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map_old/CLAM
BED_PATH=$PROJECT_PATH/BED_files
BW_PATH=$PROJECT_PATH/bigWig_files
GENOME_REF=$REF/hg19/genome/hg19_genome.fa
mkdir -p $BW_PATH

for BED in $BED_PATH/*bed.gz
do
    SAMPLENAME=$(basename ${BED%.bed.gz})
    if [[ $SAMPLENAME == "alkaline" ]]
    then
        echo python bed_coverage.py \
            $BED  \
            $GENOME_REF \
            $BW_PATH/$SAMPLENAME.bigWig \
            '+-'
    else
        for STRAND in fwd rvs
        do
            if [[ $STRAND == "fwd" ]]
            then
                S="+"
            else
                S="-"
            fi
            echo python bed_coverage.py \
                $BED  \
                $GENOME_REF \
                $BW_PATH/$SAMPLENAME.$STRAND.bigWig \
                $S
        done
    fi
done
