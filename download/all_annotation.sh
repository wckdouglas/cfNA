#!/bin/bash

REF_PATH=$REF/hg19/genome

cat $REF_PATH/RBP.bed \
| grep -v HepG2 \
| tr '_' '\t' \
| awk '{print $1,$2,$3,$4,$7,$8, "RBP",$4}' OFS='\t' \
| cat - $REF_PATH/genes.bed $REF_PATH/rmsk.bed $REF_PATH/hg19_refseq.bed \
| sort -k1,1 -k2,2n -k3,3n \
| bgzip \
> $REF_PATH/all_annotation.bed.gz
