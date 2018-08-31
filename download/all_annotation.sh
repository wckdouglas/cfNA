#!/bin/bash

REF_PATH=$REF/hg19
GENOME_PATH=$REF_PATH/genome
GENE_PATH=$REF_PATH/new_genes

cat $GENOME_PATH/RBP.bed \
| grep -v HepG2 \
| awk '$8 > 2' \
| tr '_' '\t' \
| awk '{print $1,$2,$3,$4,$10,$8, "RBP",$4}' OFS='\t' \
| cat - $GENE_PATH/genes.bed \
            $GENOME_PATH/rmsk.bed \
            $GENOME_PATH/hg19_refseq.bed \
            $GENE_PATH/piRNA.bed \
| sort -k1,1 -k2,2n -k3,3n \
| bgzip \
> $GENE_PATH/all_annotation_K562.bed.gz

cat $GENOME_PATH/RBP.bed \
| awk '$8 > 2' \
| tr '_' '\t' \
| awk '{print $1,$2,$3,$4,$10,$8, "RBP",$4}' OFS='\t' \
| cat - $GENE_PATH/genes.bed \
            $GENOME_PATH/rmsk.bed \
            $GENOME_PATH/hg19_refseq.bed \
            $GENE_PATH/piRNA.bed \
| sort -k1,1 -k2,2n -k3,3n \
| bgzip \
> $GENE_PATH/all_annotation.bed.gz
