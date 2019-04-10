#!/bin/bash

set -x

curl 'https://www.encodeproject.org/metadata/type=Experiment&assay_title=eCLIP&limit=all/metadata.tsv' \
	| grep narrow \
	| grep -v bigBed \
	| grep hg19  \
	| awk '{print $(NF-3)}' \
	> encode_files.txt

OUT_FILE=$REF/hg19_ref/genes/RBP.bed
rm $OUT_FILE
touch $OUT_FILE
for FILE_URL in $(cat encode_files.txt| grep -v meta )
do
	echo curl -L $FILE_URL \| zcat \| awk \''$2~/^[0-9]+$/'\' \>\> $OUT_FILE
done
echo bedtools sort -i $OUT_FILE \
	\| bedtools merge -i - -s -c 4 -o collapse  -delim \',\' \
	\| awk \''{print $1,$2,$3,"RBP",0,$4,$NF}'\' OFS=\''\t'\' \
	\> ${OUT_FILE%.bed}.sorted.bed
SORTED_OUT=${OUT_FILE%.bed}.sorted.bed
#rm $OUT_FILE
#touch $OUT_FILE
#for FILE_URL in $(cat encode_files.txt| grep -v meta )
#do
#	curl -L $FILE_URL | zcat | awk '$2~/^[0-9]+$/' >> $OUT_FILE
#done
#bedtools sort -i $OUT_FILE -faidx $REF/hg19/genome/hg19_genome.fa.fai \
#	> $SORTED_OUT

NO_SMRNA_RBP=${OUT_FILE%.bed}.no_smRNA.bed
#bedtools intersect \
#		-v -a $SORTED_OUT \
#		-b $REF/hg19/genome/sncRNA_x_protein.bed \
#	| bedtools intersect -a - \
#		-v -b $REF/hg19/genome/genome_rRNA.bed \
#	> $NO_SMRNA_RBP
python clean_rbp.py $NO_SMRNA_RBP 

NO_RRNA_RBP=${OUT_FILE%.bed}.no_smRNA.no_TF.bed
#bedtools intersect -v \
#		-a $NO_SMRNA_RBP \
#		-b  $REF/hg19/genome/hg19_Tfbs_clusteredV3.bed \
#		> $NO_RRNA_RBP
python clean_rbp.py $NO_RRNA_RBP
