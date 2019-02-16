import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
import numpy as np
import glob
import pandas as pd
import mygene
import os
from functools import reduce
from sequencing_tools.viz_tools import simpsons_palette
from sklearn.metrics.pairwise import cosine_similarity
from scipy.stats import pearsonr, spearmanr
from plotting_utils import *
from functools import lru_cache

def genes_annot():
    return pd.read_table('/stor/work/Lambowitz/ref/hg19_ref/genes/genes.annot',
                 names = ['gname','gtype', 'Name'])\
        .assign(Name = lambda d: d.Name.str.split('.', expand=True).iloc[:,0])

def published():
    gene_expr = '/stor/work/Lambowitz/cdw2854/cfNA/platelets/tissues/rna_tissue.tsv'
    expr_df = pd.read_table(gene_expr) \
        .pipe(pd.pivot_table, columns = 'Sample', 
              index = ['Gene','Gene name'], values = 'Value')\
        .reset_index()
    return expr_df

def TOP_df():
    excel = 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2441802/bin/gkn248_nar-02753-r-2007-File009.xls'
    return pd.read_excel(excel)

def TOP_genes():
    tdf = TOP_df()
    tg = []
    for r in tdf['Refseq ID'].str.strip(',').values:
        tg.extend(r.split(','))
    return tg

def TOP_gene_df():
    mg = mygene.MyGeneInfo()
    top_genes = TOP_genes()
    return pd.DataFrame(mg.querymany(top_genes, scopes='refseq'))\
        .pipe(lambda d: d[~pd.isnull(d.symbol)])
    
    
@lru_cache()
def tid_to_gid():
    gtf = '/stor/work/Lambowitz/ref/hg19_ref/genes/genes.gtf'
    tdf = []
    with open(gtf) as gt:
        for r in gt:
            if not r.startswith('#') and r.split('\t')[2] == 'transcript':
                gid = r.split('gene_id')[1].split(';')[0].strip(' "')
                tid = r.split('transcript_id')[1].split(';')[0].strip(' "')
                gname = r.split('gene_name')[1].split(';')[0].strip(' "')
                tdf.append((tid, gid, gname))
    return pd.DataFrame(tdf, columns=['target_id', 'gid','gname'])
tdf = tid_to_gid()

    
def read_kallisto(tpm_table):
    return pd.read_table(tpm_table) \
        .filter(['tpm','target_id','eff_length'])\
        .merge(tdf, how = 'left', on ='target_id') \
        .groupby(['gid','gname'], as_index=False)\
        .agg({'tpm':'sum','eff_length':'mean'}) \
        .rename(columns = {'tpm':os.path.basename(os.path.dirname(tpm_table))})
        
        
gene_cats = ['Others',"5' TOP", 'Ribosomal proteins', 'Histone', 'Blood','aaRS']
def label_gene(x):
    label = gene_cats[0]
    if re.search('^HIST',x):
        label = gene_cats[3]
    elif re.search('^RPL[0-9]+|^RPS[0-9]+', x) and '-' not in x:
        label = gene_cats[2]
    elif re.search('^HB[ABDEZMGQ][0-9]$|^HB[ABDEZMGQ]$|^ALB$|^B2M$|^FG[A-Z]|^S100', x):
        label = gene_cats[4]
    elif re.search('^[A-Z]ARS$', x):
        label = gene_cats[5]
    return label

gene_encoder = color_encoder()
gene_encoder.encoder = {g:col for g, col in zip(gene_cats, 
                        ['grey', 'skyblue','darkblue','#0ba518','#d82915','#edbc61'])}

def get_top_rna():
    top_df = TOP_gene_df()\
        .pipe(lambda d: d[~pd.isnull(d.symbol)])
    top_df.head()
    rp = top_df.pipe(lambda d: d[d.symbol.str.contains('^RP[LS]')]).symbol.unique().tolist()
    TOP_RNA = 'EEF1A1, EEF1B2, EEF1D, EEF1G, EEF2, EIF3A, EIF3E, EIF3F, EIF3H, EIF4B, HNRNPA1, NAP1L1, NPM1, PABPC1, RACK1,TPT1,VIM'
    TOP_RNA = TOP_RNA.split(',')
    TOP_RNA = map(lambda x: x.strip(), TOP_RNA)
    TOP_RNA = list(TOP_RNA)
    TOP_RNA.extend(rp)
    return TOP_RNA
TOP_RNA = get_top_rna()

def plot_heatmap(tpm_df, ax, var = 'Poly(A)-selected', selected = 'L[12]|Poly\(A\)|DNase|[Ff]rag', colored=True):
    top_n = 50
    top_df = tpm_df\
        .assign(Name = lambda d: d.gid.str.split('.', expand=True).iloc[:,0])\
        .rename(columns = {'Name':'Gene'})\
        .merge(published().filter(regex='Gene'))\
        .drop('Gene', axis=1)\
        .set_index('Gene name')\
        .nlargest(top_n, var)\
        .fillna(0)\
        .filter(regex=selected)\
        .transform(lambda x: np.log2(x+1))\
        .transpose() \
        .rename(index={'DNase I':'TGIRT-seq','Poly(A)-selected':'SMART-Seq'})
    sns.heatmap(top_df, ax = ax, cmap = 'inferno')    
    ax.set_xticks(np.arange(top_n) + 0.5)
    xts = ax.set_xticklabels(top_df.columns.values, fontsize=10)
    xts = ax.set_yticklabels(ax.get_yticklabels(), rotation=0, rotation_mode='anchor', ha='right')
    ax.collections[0].colorbar.set_label('$log_2$ TPM', 
                                         rotation=270, 
                                         va = 'bottom',
                                        size = 15)
    if colored:
        for yt in ax.get_xticklabels():
            if yt.get_text() in TOP_RNA:
                gene_label = "5' TOP"
            else:
                gene_label = label_gene(yt.get_text())
            
            color = gene_encoder.encoder[gene_label]
            if gene_label != "Other":
                yt.set_color(color)
    
