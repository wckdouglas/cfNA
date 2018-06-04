WORK_DIR=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/classifier/

echo deep-cfNA train \
    --train_bed_positive $WORK_DIR/train_RNA.bed \
    --train_bed_negative $WORK_DIR/train_DNA.bed \
    --validation_bed $WORK_DIR/test.bed \
    --genome $REF/hg19/genome/hg19_genome.fa \
    --model_prefix $WORK_DIR/deep_cfNA \
    --epochs 10 --steps 50000

