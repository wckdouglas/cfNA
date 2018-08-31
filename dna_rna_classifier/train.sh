WORK_DIR=/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/classifier/

for model in full_model reduced_model
do
    if [[ $model == "reduced_model" ]]
    then
        OPTION="--reduced_model" 
        MODEL_PREFIX=deep_cfNA_reduced
    else
        OPTION=" "
        MODEL_PREFIX=deep_cfNA
    fi
    echo deep-cfNA train \
        --train_bed_positive $WORK_DIR/train_RNA.bed \
        --train_bed_negative $WORK_DIR/train_DNA.bed \
        --validation_bed $WORK_DIR/test.bed \
        --genome $REF/hg19/genome/hg19_genome.fa \
        --model_prefix $WORK_DIR/$MODEL_PREFIX \
        --epochs 30 --steps 50000  $OPTION
done

