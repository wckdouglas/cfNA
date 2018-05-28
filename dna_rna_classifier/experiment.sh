PROJECT_PATH=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map
BED_PATH=$PROJECT_PATH/bed_files
OUT_PATH=$PROJECT_PATH/classifier

DEEP_BED=$OUT_PATH/untreated.classified.bed.gz
LENGTH_BED=$OUT_PATH/length.classified.bed.gz 

deep-cfNA predict \
    --inbed $BED_PATH/Qcf_All1_R1_001.bed.gz \
    --genome $REF/hg19/genome/hg19_genome.fa \
    --out $DEEP_BED \
    --model_prefix $OUT_PATH/deef_cfNA    
echo 'Predicted NA'


zcat $BED_PATH/Qcf_All1_R1_001.bed.gz \
    | awk '$3-$2 < 200 {if(($3-$2) < 80) print $0, "RNA"; else print $0, "DNA"}' OFS='\t'\
    | bgzip \
    > $LENGTH_BED


for BED in $DEEP_BED $LENGTH_BED
do
    OUT_BAM=${BED%.bed.gz}.bam
    COUNT_FILE=${OUT_BAM%.bam}.counts
    for NA in  RNA DNA
    do
        TEMP_BAM=$OUT_PATH/${NA}.temp.bam
        zcat $OUT_PATH/untreated.classified.bed.gz \
            | awk '$3-$2 < 200 {printf "%s\t%s\t%s\t%s:%s\t%s\t%s\n", $1,$2,$3,$7,$4,$5,$6 }' OFS='\t' \
            | grep --color=no $NA \
            | bedToBam -i - -g $REF/hg19/genome/hg19_genome.fa.fai \
            | samtools addreplacerg -r ID:$NA - \
            | samtools view -b \
            > $TEMP_BAM
        echo Made $TEMP_BAM
    done

    samtools merge -@ 24 - $OUT_PATH/DNA.temp.bam $OUT_PATH/RNA.temp.bam \
       | sambamba sort -t 24 -p -o $OUT_BAM /dev/stdin
    echo Written $OUT_BAM 

    for NA in RNA DNA
    do
        rm $OUT_PATH/${NA}.temp.bam
    done


    samtools view -h $OUT_BAM \
        | awk '$1~"^@" || $0~"RG:Z:RNA"' \
        | samtools view -b \
        | bedtools bamtobed -i - \
        | bedtools coverage -counts -b - \
                -a $REF/hg19/new_genes/genes.bed \
        | sort -k9,9nr \
        > $COUNT_FILE
done
