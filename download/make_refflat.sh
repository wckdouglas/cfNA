#!/bin/bash

GENE_PATH=/stor/work/Lambowitz/ref/hg19/new_genes

cat $GENE_PATH/genes.gtf | grep --color=no 'protein_coding' > $GENE_PATH/proteins.gtf
gtfToGenePred $GENE_PATH/proteins.gtf /dev/stdout | awk '{print $1,$0}' OFS='\t' > $GENE_PATH/proteins.refflat

