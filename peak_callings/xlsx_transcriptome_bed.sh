csvtk xlsx2csv transcriptome_peaks.xlsx \
      | awk -F',' '{print $12,$13,$14,$1,$8,$9,$10,$11}' OFS='\t' \
      | sed 1d \
      > transcriptome_peaks.bed


cat transcriptome_peaks.bed \
      | python project_split_transcriptome.py -i - -o - -b $REF/hg19_ref/genes/genes.bed12 \
      > $SCRATCH/transcriptome_peaks.genomics.bed

cat $SCRATCH/transcriptome_peaks.genomics.bed \
    | bedtools intersect -b $REF/hg19_ref/genes/RBP_IDR.reformated.bed.gz -a - -wa \
    > $SCRATCH/transcriptome_peaks.genomics.wao.bed
