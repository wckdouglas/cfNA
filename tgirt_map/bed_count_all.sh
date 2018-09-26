#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map
COUNT_PATH=$PROJECT_PATH/Counts/all_counts
REF_BED_PATH=$REF/hg19/new_genes
mkdir -p $COUNT_PATH

for SAMPLE_NAME in $(ls $PROJECT_PATH | egrep '001$')
do
	SAMPLE_FOLDER=$PROJECT_PATH/$SAMPLE_NAME
    ALL_ALN_PATH=$SAMPLE_FOLDER/Combined
    SMALL_RNA_PATH=$SAMPLE_FOLDER/smallRNA
    rRNA_PATH=$SAMPLE_FOLDER/rRNA_mt

    #Count all genes
    for BED in $ALL_ALN_PATH/primary_no_sncRNA_tRNA_rRNA_repeats.bam \
                $ALL_ALN_PATH/sncRNA.bam \
                $SMALL_RNA_PATH/aligned.bam \
                $rRNA_PATH/aligned.bam \
                $ALL_ALN_PATH/repeats.bam 
    do
        if [[ $BED == *bam ]]
        then
            BAM_TO_BED="| samtools view -bF 4 -F 256 -F 2048 "
            BAM_TO_BED="$BAM_TO_BED  | bam_to_bed.py -i - --add_cigar --primary "
        else
            BAM_TO_BED=" "
        fi
        
        TEMP_FOLDER=$PROJECT_PATH/${SAMPLE_NAME}.$(basename $(dirname $BED)).$(basename $BED)_TEMP
        TEMP_BED=${BED%.bam}.bed.gz
        BASE_COMMAND="mkdir -p $TEMP_FOLDER ; cat $BED $BAM_TO_BED\
            | sort -k1,1 -k2,2n -k3,3n -k6,6 --temporary-directory=$TEMP_FOLDER \
            | bgzip \
            > $TEMP_BED "
        COMMAND=""

        for COUNT_TYPE in all dedup
        do
            if [[ $COUNT_TYPE == "dedup" ]]
            then
                DEDUP_COMMAND=" | deduplicate_bed.py -i - -d '_' -f 0  "
                COUNT_COMMAND=" | cut -f7-14 "
                if echo $BED | egrep -q 'sncRNA.bam$|aligned.bam$'
                then
                    COUNT_COMMAND=" | cut -f7-14 "
                fi

                if echo $BED | egrep -q 'genome-sim'
                then
                    DEDUP_COMMAND="| sort -k1,1 -k2,2n -k3,3n -k6,6n -u "
                fi

            else
                DEDUP_COMMAND=" "
                COUNT_COMMAND=" | cut -f8-15 "
                if echo $BED | egrep -q 'smallRNA/aligned.bam$|rRNA_mt/aligned.bam$'
                then
                    COUNT_COMMAND="| cut-f8-13 "
                fi
            fi
        
            SUM_COMMAND="| awk  '{print \$2,\$3,\$4,\$5,\$6,\$7,\$8, \$9, \$1}'   OFS='\t'  "
            if  echo $BED | grep -q "_no_sncRNA_"
            then
                REF_BED=$REF_BED_PATH/genes.bed
                OUT_PATH=$COUNT_PATH/counts
            elif echo $BED | egrep -q 'sncRNA.bam$'
            then
                REF_BED=$REF_BED_PATH/sncRNA_x_protein.bed
                OUT_PATH=$COUNT_PATH/sncRNA
            elif echo $BED | egrep -q 'smallRNA/aligned.bam$'
            then
                REF_BED=$REF_BED_PATH/smallRNA.bed
                OUT_PATH=$COUNT_PATH/smallRNA
            elif echo $BED | egrep -q 'rRNA_mt/aligned.bam$'
            then
                REF_BED=$REF_BED_PATH/rRNA_mt.bed
                OUT_PATH=$COUNT_PATH/rRNA_mt
            elif echo $BED | egrep -q 'repeats.bam$'
            then
                REF_BED=$REF/hg19/genome/rmsk.bed.gz
                OUT_PATH=$COUNT_PATH/repeats
            fi
            mkdir -p $OUT_PATH

            # UMI adjust
            if [[ $COUNT_TYPE == "dedup"  ]]
            then
                if echo $BED | egrep -q 'sncRNA.bam$|aligned.bam$' 
                then
                    ADJUST_UMI=' -t 0 --ct 6 |  poisson_umi_adjustment.py -i - -o - --umi 6 '
                    if echo $BED | egrep -q 'genome-sim'
                    then
                        ADJUST_UMI=' '
                    fi
                else
                    ADJUST_UMI=' '
                fi
            else
                ADJUST_UMI=' '
            fi


            for STRAND in sense antisense
            do
                if [[ $STRAND == "sense" ]]
                then
                    STRANDENESS=" -s "
                else
                    STRANDENESS=" -S "
                fi
                COMMAND="$COMMAND; zcat $TEMP_BED \
                    $DEDUP_COMMAND $ADJUST_UMI \
                    | bedtools intersect -a - -b $REF_BED $STRANDENESS -wao \
                    | bgzip \
                    | tee ${BED%.bam}.${COUNT_TYPE}.${STRAND}.intersected.bed.gz  \
                    | zcat \
                    $COUNT_COMMAND \
                    | sort --temporary-directory=$TEMP_FOLDER \
                    | uniq -c \
                     $SUM_COMMAND \
                    > $OUT_PATH/${SAMPLE_NAME}.${COUNT_TYPE}.${STRAND}.counts"
            done
        done
        echo "$BASE_COMMAND $COMMAND; rm -rf $TEMP_FOLDER"
    done
done | egrep 'Q[Cc][Ff]|genome'
