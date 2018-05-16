#!/bin/bash

PROJECT_PATH=$SCRATCH/cell_Free_nucleotides/tgirt_map
WPS_PATH=$PROJECT_PATH/merged_bed/genome_WPS
PEAK_PATH=$WPS_PATH/genome_peaks
mkdir -p $PEAK_PATH

for WPS in $WPS_PATH/*bigWig
do
	SAMPLENAME=$(basename $WPS)
	PEAK_FILE=$PEAK_PATH/${SAMPLENAME/.bigWig/.bed}
	STRAND=$(echo $SAMPLENAME | rev |cut -d'.' -f2 |rev)
    TREATMENT=$(echo ${SAMPLENAME} | cut -d'.' -f1 )

    STRANDENESS=''
    if [[ $STRAND == "fwd" ]]
    then
        STRANDENESS='forward'
    elif [[ $STRAND == "rvs" ]]
    then
        STRANDENESS='reverse'
    fi

	echo python peak_caller.py \
		--in_bigwig=$WPS  \
		--out_bed=$PEAK_FILE \
        --control_bigwig=${WPS/$TREATMENT/alkaline} \
        --strand=$STRANDENESS \
        --two_pass 
done
