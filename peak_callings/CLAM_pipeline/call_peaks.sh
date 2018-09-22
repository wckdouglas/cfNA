#!/bin/bash

CLAM_PATH=/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/CLAM
BIGWIG_PATH=$CLAM_PATH/bigWig_files
OUT_PATH=$CLAM_PATH/peaks
mkdir -p $OUT_PATH

for SAMPLE in unfragmented untreated
do
    for STRAND in fwd rvs
    do
        if [[ $STRAND == "fwd" ]]
        then
            S=forward
        else
            S=reverse
        fi
        echo python peak_caller.py \
            --in_bigwig $BIGWIG_PATH/${SAMPLE}.${STRAND}.bigWig  \
            --control $BIGWIG_PATH/alkaline.bigWig \
            -o $OUT_PATH/${SAMPLE}.${STRAND}.bed \
            --threads 6 \
            --two_pass \
            -s $S 

    done
done
