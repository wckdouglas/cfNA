REF_PATH=$REF/hg19
ANNOTATION_PATH=$REF_PATH/new_genes
GENOME_PATH=$REF_PATH/genome
GTF_LINK=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh37_mapping/gencode.v28lift37.annotation.gtf.gz
tRNA_REF=http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi19/hg19-tRNAs.tar.gz

#annotationes
curl $GTF_LINK |zcat > $ANNOTATION_PATH/genes.gtf
python gtf_to_bed.py $ANNOTATION_PATH/genes.gtf > $ANNOTATION_PATH/genes.bed

#tRNA
curl $tRNA_REF > $ANNOTATION_PATH/tRNA.tar.gz
mkdir -p $ANNOTATION_PATH/tRNA
tar zxvf $ANNOTATION_PATH/tRNA.tar.gz --directory $ANNOTATION_PATH/tRNA
python make_tRNA.py \
    $ANNOTATION_PATH/tRNA/hg19-tRNAs-detailed.ss \
    $ANNOTATION_PATH/tRNA.bed \
    $ANNOTATION_PATH/tRNA/nucleo_tRNA.fa
cat $ANNOTATION_PATH/tRNA.bed |cut -f1-8 >> $ANNOTATION_PATH/genes.bed
cat $ANNOTATION_PATH/genes.bed \
    | grep 'Mt_tRNA' \
    | bedtools getfasta  -fi $GENOME_PATH/hg19_genome.fa -bed - -s -name -tab \
    | tr ':' '\t' \
    | awk '{printf ">%s\n%s\n",$1,$NF}' \
    | sed 's/(-)//g' | sed 's/(+)//g' \
    > $ANNOTATION_PATH/tRNA/mt_tRNA.fa
cat $ANNOTATION_PATH/tRNA/mt_tRNA.fa $ANNOTATION_PATH/tRNA/nucleo_tRNA.fa > $ANNOTATION_PATH/tRNA.fa





#make rRNA
echo 'gi|23898|emb|X12811.1|  274     394     5S_rRNA 0       +       5S_rRNA 5S_rRNA
gi|555853|gb|U13369.1|HSU13369  3657    5527    18S_rRNA        0       +       18S_rRNA        18S_rRNA
gi|555853|gb|U13369.1|HSU13369  6623    6779    5.8S_rRNA       0       +       5.8S_rRNA       5.8S_rRNA
gi|555853|gb|U13369.1|HSU13369  7935    12969   28S_rRNA        0       +       28S_rRNA        28S_rRNA' \
| awk '{print $1,$2,$3,$4,$5,$6,"rDNA",$8}' OFS='\t' \
> $ANNOTATION_PATH/rRNA.bed
cat $ANNOTATION_PATH/rRNA.bed >> $ANNOTATION_PATH/genes.bed

python split_bed_for_count.py $ANNOTATION_PATH

#make rRNA tRNA
python get_rRNA_fa.py > $ANNOTATION_PATH/rRNA.fa
cat $ANNOTATION_PATH/tRNA.fa $ANNOTATION_PATH/rRNA.fa \
    > $ANNOTATION_PATH/tRNA_rRNA.fa 
echo made tRNA_rRNA fasta
bowtie2-build $ANNOTATION_PATH/tRNA_rRNA.fa $ANNOTATION_PATH/tRNA_rRNA
bowtie2-build $ANNOTATION_PATH/rRNA.fa $ANNOTATION_PATH/rRNA
bowtie2-build $ANNOTATION_PATH/tRNA.fa $ANNOTATION_PATH/tRNA

