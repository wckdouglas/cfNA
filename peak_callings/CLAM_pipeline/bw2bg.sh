#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map_old/CLAM
BW_PATH=$PROJECT_PATH/bigWig_files
BG_PATH=$PROJECT_PATH/bedgraph_files
mkdir -p $BG_PATH

for BW in $BW_PATH/*bigWig
do
    SAMPLENAME=$(basename ${BW%.bigWig})
    echo bigWigToBedGraph $BW $BG_PATH/${SAMPLENAME}.bedGraph
done
