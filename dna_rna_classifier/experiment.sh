PROJECT_PATH=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map
BED_PATH=$PROJECT_PATH/bed_files
OUT_PATH=$PROJECT_PATH/classifier

DEEP_BED=$OUT_PATH/untreated.classified.bed.gz
OUT_BAM=${DEEP_BED%.bed.gz}.bam
COUNT_FILE=${OUT_BAM%.bam}.counts
LENGTH_COUNT=${OUT_PATH%.bam}.length.counts
LENGTH_BED=$OUT_PATH/untreated.classified.length_annotated.bed.gz



deep-cfNA predict \
    --inbed $BED_PATH/Qcf_All1_R1_001.bed.gz \
    --genome $REF/hg19/genome/hg19_genome.fa \
    --out $DEEP_BED \
    --model_prefix $OUT_PATH/deef_cfNA    
echo 'Predicted NA'


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

zcat $DEEP_BED \
    | awk  '{if(($3-$2) < 80) printf "%s\tRNA\n",$0 ; else printf "%s\tDNA\n",$0;}' OFS='\t' \
    | bgzip \
    > $LENGTH_BED

zcat $LENGTH_BED \
    | awk '$(NF-1)=="RNA"' \
    | bedtools coverage -b - -a $REF/hg19/genome/genes.bed -counts \
    | sort -k9,9nr \
    > $COUNT_FILE

zcat $LENGTH_BED \
    | awk '$(NF)=="RNA"' \
    | bedtools coverage -b - -a $REF/hg19/genome/genes.bed -counts \
    | sort -k9,9nr \
    > $LENGTH_COUNT





