#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map
COUNT_PATH=$PROJECT_PATH/Counts/all_counts
REF_BED_PATH=$REF/hg19/new_genes
mkdir -p $COUNT_PATH

for SAMPLE_NAME in $(ls $PROJECT_PATH | egrep '001$')
do
	SAMPLE_FOLDER=$PROJECT_PATH/$SAMPLE_NAME
    ALL_ALN_PATH=$SAMPLE_FOLDER/Combined
    tRNA_PATH=$SAMPLE_FOLDER/tRNA
    rRNA_PATH=$SAMPLE_FOLDER/rRNA

    #Count all genes
    for BED in $ALL_ALN_PATH/primary_no_sncRNA_tRNA_rRNA_repeats.bam \
                $ALL_ALN_PATH/sncRNA.bam \
                $tRNA_PATH/tRNA_remap.bam \
                $rRNA_PATH/rRNA_remap.bam \
                $ALL_ALN_PATH/repeats.bam 
    do
        if [[ $BED == *bam ]]
        then
            BAM_TO_BED=" | bam_to_bed.py -i - --add_cigar  "
        else
            BAM_TO_BED=" "
        fi
        
        TEMP_FOLDER=$PROJECT_PATH/${SAMPLE_NAME}.$(basename $BED)_TEMP
        TEMP_BED=${BED%.bam}.bed.gz
        BASE_COMMAND="mkdir -p $TEMP_FOLDER; cat $BED $BAM_TO_BED\
            | sort -k1,1 -k2,2n -k3,3n -k6,6 --temporary-directory=$TEMP_FOLDER \
            | bgzip \
            > $TEMP_BED "
        COMMAND=$BASE_COMMAND

        for COUNT_TYPE in all dedup
        do
            if [[ $COUNT_TYPE == "dedup" ]]
            then
                DEDUP_COMMAND=" | deduplicate_bed.py -i - -d '_' -f 0  "
                COUNT_COMMAND=" | awk '{print \$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14}' OFS='\t'"
                if echo $BED | egrep -q 'tRNA_remap.bam$'
                then
                    COUNT_COMMAND="| awk '{print \$7,\$8,\$9,\$10,\$11,\$12}' OFS='\t'"
                fi

            else
                DEDUP_COMMAND=" "
                COUNT_COMMAND=" | awk '{print \$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15}' OFS='\t'"
                if echo $BED | egrep -q 'tRNA_remap.bam$'
                then
                    COUNT_COMMAND="| awk '{print \$8,\$9,\$10,\$11,\$12, \$13}' OFS='\t'"
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
            elif echo $BED | egrep -q 'tRNA_remap.bam$'
            then
                REF_BED=$REF_BED_PATH/tRNA_yRNA.count.bed
                OUT_PATH=$COUNT_PATH/tRNA
                SUM_COMMAND=" |awk  '{print \$2,\$3,\$4,\$5,\$6,\$7, \$1}'   OFS='\t'  "
            elif echo $BED | egrep -q 'rRNA_remap.bam$'
            then
                REF_BED=$REF_BED_PATH/rRNA.bed
                OUT_PATH=$COUNT_PATH/rRNA
            elif echo $BED | egrep -q 'repeats.bam$'
            then
                REF_BED=$REF/hg19/genome/rmsk.bed.gz
                OUT_PATH=$COUNT_PATH/repeats
            fi
            mkdir -p $OUT_PATH

            if [[ $COUNT_TYPE == "dedup"  ]]
            then
                if echo $BED | egrep -q 'sncRNA.bam$|tRNA_remap.bam$'
                then
                    ADJUST_UMI=' -t 0 --ct 6|  poisson_umi_adjustment.py -i - -o - --umi 6 '
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
                    | bedtools intersect -a - -b $REF_BED -f 0.1 $STRANDENESS -wao \
                    | bgzip \
                    | tee ${BED%.bam}.${COUNT_TYPE}.${STRAND}.intersected.bed.gz  \
                    | zcat \
                    $COUNT_COMMAND \
                    | sort \
                    | uniq -c \
                     $SUM_COMMAND \
                    > $OUT_PATH/${SAMPLE_NAME}.${COUNT_TYPE}.${STRAND}.counts"
            done
        done
        echo "$COMMAND; rm -rf $TEMP_FOLDER"
    done
done | egrep 'Q[Cc][Ff]|genome'
