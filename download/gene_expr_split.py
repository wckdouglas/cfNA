#!/usr/bin/env python

import pandas as pd
from collections import defaultdict

gene_path = '/stor/work/Lambowitz/ref/hg19/new_genes'

def split_bone_marrow():
    gene_expr = '/stor/work/Lambowitz/cdw2854/EV_polyA/published_expr/rna_tissue.tsv'
    expr_df = pd.read_table(gene_expr)  \
            .query('Sample == "bone marrow"') \
            .assign(expr_level = lambda d: pd.qcut(d['Value'], 3, labels=False))  \
            .assign(expr_level = lambda d: d.expr_level.map({0:'Low',1:'Medium',2:'High'}))
    return expr_df


def make_gene_dict(expr_df):
    gene_dict = {}
    for i, row in expr_df.iterrows():
        gene_dict[row['Gene']] = row['expr_level']
    return  gene_dict


def  open_files(levels):
    file_dict = {}
    for l in levels:
        filename = gene_path + '/bone_marrow.%s.bed' %l
        file_dict[l] = open(filename, 'w')
    return file_dict
        
expr_df = split_bone_marrow()
levels = expr_df.expr_level.unique()
file_dict = open_files(levels)
gene_dict = make_gene_dict(expr_df)

expr_count = defaultdict(int)
skipped = 0
with open(gene_path + '/genes.bed12.bed', 'r') as inbed:
    for line in inbed:
        line = line.strip()
        gene_id = line.split('\t')[3]
        if gene_id in gene_dict.keys():
            expr = gene_dict[gene_id]
            expr_count[expr] += 1
            print(line, file = file_dict[expr])
        else:
            skipped +=1
    
for key, value in expr_count.items():
    print('%s: %i genes' %(key, value))

print('skipped %i genes ' %skipped)