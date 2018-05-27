PROJECT_PATH=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map
BED_PATH=$PROJECT_PATH/bed_files
OUT_PATH=$PROJECT_PATH/classifier

deep-cfNA predict \
    --inbed $BED_PATH/Qcf_All1_R1_001.bed.gz \
    --genome $REF/hg19/genome/hg19_genome.fa \
    --out $OUT_PATH/untreated.classified.bed.gz \
    --model_prefix $OUT_PATH/deef_cfNA    
echo 'Predicted NA'

for NA in  RNA DNA
do
    zcat $OUT_PATH/untreated.classified.bed.gz \
        | awk '$3-$2 < 200 {printf "%s\t%s\t%s\t%s:%s\t%s\t%s\n", $1,$2,$3,$7,$4,$5,$6 }' OFS='\t' \
        | grep --color=no $NA \
        | bedToBam -i - -g $REF/hg19/genome/hg19_genome.fa.fai \
        | bam_umi_tag.py -i - \
        | samtools addreplacerg -r ID:$NA - \
        | samtools view -b \
        > $OUT_PATH/${NA}.classified.bam
    echo Made $OUT_PATH/${NA}.classified.bam
done

samtools merge -@ 24 - $OUT_PATH/DNA.classified.bam $OUT_PATH/RNA.classified.bam \
   | sambamba sort -t 24 -p -o $OUT_PATH/classified.bam /dev/stdin
echo Written $OUT_PATH/classified.bam 

for NA in RNA DNA
do
    rm $OUT_PATH/${NA}.classified.bam
done


samtools view -h $OUT_PATH/classified.bam \
    |  awk '$1~"^@" || $0~"RG:Z:RNA"' \
    | samtools view -b \
    | bedtools bamtobed -i - \
    | bedtools coverage -counts -b - \
            -a $REF/hg19/new_genes/genes.bed \
    | sort -k9,9nr \
    > $OUT_PATH/classified.counts
