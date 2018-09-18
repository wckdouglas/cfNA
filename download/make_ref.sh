REF_PATH=$REF/hg19
ANNOTATION_PATH=$REF_PATH/new_genes
GENOME_PATH=$REF_PATH/genome
GTF_LINK=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh37_mapping/gencode.v28lift37.annotation.gtf.gz
tRNA_REF=http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi19/hg19-tRNAs.tar.gz
piRNA=http://www.regulatoryrna.org/database/piRNA/download/archive/v1.0/bed/piR_hg19_v1.0.bed.gz

#annotationes
curl $GTF_LINK |zcat > $ANNOTATION_PATH/genes.gtf
hisat2_extract_splice_sites.py $ANNOTATION_PATH/genes.gtf > $ANNOTATION_PATH/splicesites.tsv
python gtf_to_bed.py $ANNOTATION_PATH/genes.gtf > $ANNOTATION_PATH/genes.bed


#piRNA
curl $piRNA \
    | zcat \
    | sort -k1,1 -k2,2n -k3,3n \
    | awk '{print $0, "piRNA","piRNA"}' OFS='\t' \
    | bgzip \
    > $ANNOTATION_PATH/piRNA.bed.gz

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
python get_rRNA_fa.py $GENOME_PATH/hg19_genome.fa > $ANNOTATION_PATH/rRNA.fa
echo 'gi|23898|emb|X12811.1|  274     394     5S_rRNA 0       +       5S_rRNA 5S_rRNA
gi|555853|gb|U13369.1|HSU13369  3657    5527    18S_rRNA        0       +       18S_rRNA        18S_rRNA
gi|555853|gb|U13369.1|HSU13369  6623    6779    5.8S_rRNA       0       +       5.8S_rRNA       5.8S_rRNA
gi|555853|gb|U13369.1|HSU13369  7935    12969   28S_rRNA        0       +       28S_rRNA        28S_rRNA' \
| awk '{print $1,$2,$3,$4,$5,$6,"rDNA",$8}' OFS='\t' \
> $ANNOTATION_PATH/rRNA.bed

echo 'MT-RNR1 0 953 MT-RNR1 0 + mt_rRNA MT-RNR1
MT-RNR2 0 1558 MT-RNR2 0 + mt_rRNA MT-RNR2' \
| awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' OFS='\t' \
>> $ANNOTATION_PATH/rRNA.bed

cat $ANNOTATION_PATH/rRNA.bed >> $ANNOTATION_PATH/genes.bed
zcat $ANNOTATION_PATH/piRNA.bed.gz >> $ANNOTATION_PATH/genes.bed

python split_bed_for_count.py $ANNOTATION_PATH


#make tRNA filter
cat $ANNOTATION_PATH/tRNA.bed $ANNOTATION_PATH/rmsk_tRNA.bed $REF_PATH/genome/tRNA.bed \
    | bedtools sort \
    | bedtools merge -s -o first -c 4,5,6,7,8\
    > $ANNOTATION_PATH/tRNA_comprehensive.bed

#yRNA
cat $ANNOTATION_PATH/genes.bed \
| grep --color=no 'RNY' \
| bedtools getfasta -bed - -fi $GENOME_PATH/hg19_genome.fa -s -name \
| seqkit seq -u \
| seqkit rmdup -w 1000 -s  \
> $ANNOTATION_PATH/yRNA.fa 

#make rRNA tRNA
cat $ANNOTATION_PATH/tRNA.fa $ANNOTATION_PATH/rRNA.fa $ANNOTATION_PATH/yRNA.fa\
    > $ANNOTATION_PATH/tRNA_rRNA_yRNA.fa 

cat $ANNOTATION_PATH/tRNA.fa $ANNOTATION_PATH/yRNA.fa > $ANNOTATION_PATH/tRNA_yRNA.fa
cat $ANNOTATION_PATH/genes.bed | awk '$4~"RNY|Y_RNA"' > $ANNOTATION_PATH/yRNA.bed
cat $ANNOTATION_PATH/yRNA.bed $ANNOTATION_PATH/tRNA.bed > $ANNOTATION_PATH/tRNA_yRNA.bed

echo made tRNA_rRNA fasta
bowtie2-build $ANNOTATION_PATH/tRNA_rRNA_yRNA.fa $ANNOTATION_PATH/tRNA_rRNA_yRNA
bowtie2-build $ANNOTATION_PATH/rRNA.fa $ANNOTATION_PATH/rRNA
bowtie2-build $ANNOTATION_PATH/tRNA.fa $ANNOTATION_PATH/tRNA
bowtie2-build $ANNOTATION_PATH/yRNA.fa $ANNOTATION_PATH/yRNA
bowtie2-build $ANNOTATION_PATH/tRNA_yRNA.fa $ANNOTATION_PATH/tRNA_yRNA
