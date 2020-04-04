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


def get_tpm_df(return_files = False):
    #kallisto_path = work_path + '/cfNA/tgirt_map/kallisto_protein_result'
    kallisto_path = '/stor/work/Lambowitz/yaojun/cfNA/tgirt_map/kallisto_protein_result'
    sample_folders = glob.glob(kallisto_path + '/*')
    sample_folders.sort()
    sample_folders = filter(lambda x: re.search('MP|PP|[Qq][cC][fF]', x), sample_folders)
    kallisto_tpm = list(map(lambda x: x + '/abundance.tsv', sample_folders))
    if return_files:
        return kallisto_tpm
    tpm_dfs = map(read_kallisto, kallisto_tpm)
    tpm_dfs = map(lambda d: d.drop(['eff_length'], axis=1), tpm_dfs)
    tpm_df = reduce(lambda x,y: x.merge(y, how = 'outer', on = ['gname','gid']), tpm_dfs) \
        .pipe(pd.melt, id_vars = ['gname','gid'], var_name ='samplename', value_name = 'TPM') \
        .assign(prep = lambda d: d.samplename.map(label_sample)) \
        .groupby(['gname','gid','prep'], as_index=False)\
        .agg({'TPM': 'mean'}) \
        .assign(TPM = lambda d: d.groupby('prep')['TPM'].transform(lambda x: x/x.sum() * 1e6))\
        .pipe(pd.pivot_table, columns = 'prep', values='TPM', index=['gid','gname'], fill_value = 0)\
        .reset_index()
    return tpm_df



def make_gene_df(tpm_df):
    gene_df = tpm_df \
        .assign(gene_label = lambda d: d.gname.map(label_gene))\
        .assign(top = lambda d: np.where(d.gname.isin(TOP_RNA), 'is_top','not_top'))\
        .assign(gene_label = lambda d: np.where(d.top=='is_top', "5' TOP", d.gene_label))\
        .assign(gene_label = lambda d: d.gene_label.astype(pd.api.types\
                                    .CategoricalDtype(gene_cats)))\
        .assign(color = lambda d: d.gene_label.map(gene_encoder.encoder))
    return gene_df

def published():
    #gene_expr = work_path + '/cfNA/platelets/tissues/rna_tissue.tsv'
    gene_expr = '/stor/work/Lambowitz/yaojun/cfNA/platelets/tissues/rna_tissue.tsv'
    expr_df = pd.read_table(gene_expr) \
        .pipe(pd.pivot_table, columns = 'Sample', 
              index = ['Gene','Gene name'], values = 'Value')\
        .reset_index()
    return expr_df

def TOP_df():
    excel = 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2441802/bin/gkn248_nar-02753-r-2007-File009.xls'
    return pd.read_excel(excel)

@lru_cache(maxsize=1000000)
def TOP_genes():
    tdf = TOP_df()
    tg = []
    for r in tdf['Refseq ID'].str.strip(',').values:
        tg.extend(r.split(','))
    return tg

@lru_cache()
def TOP_gene_df():
    mg = mygene.MyGeneInfo()
    top_genes = TOP_genes()
    return pd.DataFrame(mg.querymany(top_genes, scopes='refseq'))\
        .pipe(lambda d: d[~pd.isnull(d.symbol)])
    
    
@lru_cache(maxsize=1000000)
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

@lru_cache()
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
    
def coloring_gene_dot(gt_df, gt, xn, yn, ax):

    offsets = {'HBZ': (-0.2, -0.4),
                'HBG2': (-0.2, 0.1),
                'HBG1': (0.5, 0.2),
                'HBD': (0.5, 0.2),
                'HBQ1': (0.3, 0.3),
                'HBM': (0.2,-0.6),
                'HBE1': (0, -0.3),
                'HBA2': (0.5, 0.1),
                'HBB': (0.3, -0.3),
                'HBA1': (-0.5, 0.3),
                'S100A9': (0.3, 0),
                'S100A8': (0.3, -0.2),
                'HBD': (0.5,-0.4),
                'HBM': (0.4, - 0.2),
                'HBG2': (1.2, -1),
                'HBE1': (0,-0.1)}
        
    if gt == 'Blood':
        for i, row in gt_df.pipe(lambda d: d[d.gname.str.contains('^HB[A-Z]$|^HB[A-Z][0-9]+$|^S100A[89]')]).iterrows():
            if row['gname'] in offsets.keys():
                xoffset, yoffset = offsets[row['gname']]
            else:
                xoffset, yoffset = 0,0
            if row[xn] > 0:
                x = np.log10(row[xn]+1)
                y = np.log10(row[yn]+1)
                ax.annotate(s = row['gname'],
                            xy = (x,y), 
                            xytext = (x+xoffset, y+yoffset),
                            fontsize=13, 
                            color = gene_encoder.encoder['Blood'],
                            arrowprops = {'arrowstyle':'-', 
                                        'color':'red'},
                            weight = 'bold')
    
    lgd = gene_encoder.show_legend(ax, loc='upper left', 
                               frameon=False, fontsize=18)
    lgd.set_title(title = '', prop={'size':18})


def plot_scatter_kallisto(gene_df, xn, yn, ax, 
                        marginal_ax = (None, None),
                        gene_label=False,
                        cor_value = True):
    
    ax_xmarginal, ax_ymarginal = marginal_ax
    for (gt, col), gt_df in gene_df.groupby(['gene_label','color']):
        xv = np.log10(gt_df[xn].fillna(0)+1)
        yv = np.log10(gt_df[yn].fillna(0)+1)
        
        alpha = 0.5 if gt == 'Others' else 1
        size = 20 if gt == 'Others' else 20
        ax.scatter(xv, 
                yv, 
                s = size,
                color = col, 
                alpha = alpha)
        if ax_xmarginal and ax_ymarginal:
            sns.kdeplot(xv, ax = ax_xmarginal, color = col, cut = 0)
            sns.kdeplot(yv, ax = ax_ymarginal, color = col, vertical=True, cut=0)

        if gene_label:
            coloring_gene_dot(gt_df, gt, xn, yn, ax)



    r, _ = pearsonr(np.log10(gene_df[xn]+1), np.log10(gene_df[yn]+1))
    if cor_value:
        ax.text(5,1, "Spearman's\n" + r"$\rho$ = %.2f" %(r), fontsize=18)
    #p.ax_marg_y.set_visible(False)
    #p.ax_marg_x.set_visible(False)

    lmax = 8
    lmin = -0.5

    if ax_xmarginal and ax_ymarginal:
        ax_ymarginal.set_ylim(lmin,lmax)
        ax_xmarginal.set_xlim(lmin, lmax)
        ax_xmarginal.legend().set_visible(False)
        ax_xmarginal.legend().set_visible(False)
        ax_xmarginal.set_xlabel('')
        ax_ymarginal.set_ylabel('')

    ax.set_xlim(lmin, lmax)
    ax.set_ylim(lmin, lmax)
    ax.set_xticks(np.arange(0,lmax,1))
    ax.set_yticks(np.arange(0,lmax,1))
    ax.plot([0,7],[0,7], color='red')
