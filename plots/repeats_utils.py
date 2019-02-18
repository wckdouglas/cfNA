
import pandas as pd
import numpy as np
from numpy import log, exp
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os 
import random
from scipy.stats import binom_test
from scipy.stats import beta
from scipy.special import betaln
from functools import partial
from sequencing_tools.fastq_tools import reverse_complement
from collections import defaultdict



def fix_simple_repeats(p_df):
    simple_repeat_count = p_df\
        .pipe(lambda d: d[d.gene_name.str.contains('Simple_repeat')]) \
        .assign(oligo = lambda d: d.gene_name.str.extract('([ACTG]+)', expand=True).iloc[:,0]) \
        .assign(rev_oligo = lambda d: d.oligo.map(reverse_complement)) \
        .groupby(['oligo','rev_oligo'], as_index=False)\
        .agg({'sense':'sum',
             'antisense':'sum'})
    
    repeat_count_dict = defaultdict(lambda: defaultdict(int))
    
    for i, row in simple_repeat_count.iterrows():
        repeat_count_dict[row['oligo']]['sense'] = row['sense']
        repeat_count_dict[row['oligo']]['antisense'] = row['antisense']
    
    
    return simple_repeat_count\
        .assign(rev_sense = lambda d: d.rev_oligo.map(lambda x: repeat_count_dict[x]['sense'])) \
        .assign(rev_antisense = lambda d: d.rev_oligo.map(lambda x: repeat_count_dict[x]['antisense'])) \
        .assign(combine_sense = lambda d: d.sense + d.rev_antisense) \
        .assign(combine_antisense = lambda d: d.antisense + d.rev_sense) \
        .assign(gene_id = 'Simple_repeats')\
        .assign(gene_name = lambda d: 'Simple_repeats:(' + d.oligo + ')n' ) \
        .filter(regex = 'gene|combi')\
        .rename(columns = {'combine_sense':'sense',
                          'combine_antisense': 'antisense'})


def get_repeat_df(df, sample_regex=None):
    p_df = df\
        .pipe(lambda d: d[~d.gene_name.str.contains('No features')])\
        .query('grouped_type == "Repeats"')\
        .filter(regex = 'gene|grouped|%s' %(sample_regex)) \
        .pipe(pd.melt, 
              id_vars = ['gene_id','gene_name','gene_type','grouped_type'],
              var_name = 'samplename', 
              value_name = 'read_count') \
        .assign(strand = lambda d: d.samplename.str.split(':',expand=True).iloc[:, -1]) \
        .groupby(['gene_id','gene_name','strand'], as_index=False)\
        .agg({'read_count':'sum'}) \
        .pipe(pd.pivot_table, 
              index=['gene_id','gene_name'], 
              columns = 'strand', 
              values='read_count') \
        .reset_index() \
        .pipe(lambda d: d[~d.gene_id.str.contains('rRNA|tRNA|RNA')]) \
        .query('(antisense + sense) > 0') \
        .pipe(lambda d: pd.concat([d.pipe(lambda d: d[d.gene_name.str.contains('Satellite')]),
                              fix_simple_repeats(d)]))
    return p_df
    
def model_df(p_df, ax, title=''):
    # construct prior
    cov_cutoff = 100
    hi_df = p_df\
        .query('(sense + antisense) > %i' %cov_cutoff)\
        .assign(percentage_sense = lambda d: d.sense/(d.sense + d.antisense))
    fitted_params = beta.fit(data = hi_df.percentage_sense.values, floc=0, fscale=1 )
    print(fitted_params)
    alpha0, beta0, loc, scale  = fitted_params
    
    hist = True
    bins = 10
    sns.distplot(hi_df.percentage_sense, 
                 label = 'High count repeats (>%i fragments)' %cov_cutoff, 
                 ax =ax, hist_kws={'alpha':0.5}, 
                 bins=bins, 
                 hist=hist)
    ls = np.linspace(0,1, 1000)
    ax.plot(ls,beta.pdf(ls, alpha0, beta0, loc, scale), 
                label = 'Fitted beta-binomial')
    #sns.distplot(np.random.beta(alpha0, beta0, size=hi_df.shape[0]))
    ax.legend(frameon=False, bbox_to_anchor=(0.6,1.1))
    ax.set_xlabel('Propotion of sense strand fragments')
    ax.set_ylabel('Density')
    ax.set_xlim(0, 1)
    ax.set_title(title, fontsize=15)
    sns.despine()
    return alpha0, beta0


def update_beta(alpha0, beta0, sense, antisense, param = 'success_rate'):
    if param == "success_rate":
        return (alpha0 + sense) / (alpha0 + beta0 + sense + antisense)
    
    elif param == "alpha":
        return sense + alpha0
    
    elif param == "beta":
        return antisense + beta0
                          
def cal_beta_binom(alpha0, beta0, sense,antisense):
    '''
    https://www.johndcook.com/blog/2015/03/31/bayes-factors-vs-p-values/

    logbf = cal_beta_binom(alpha0, beta0, sense, antisense)
    '''
    N = sense + antisense
    logbf =  betaln(alpha0, beta0) - betaln(alpha0 + sense, beta0 + antisense)
    return logbf


def update_empirical_bayesr(p_df, alpha0, beta0):
    return p_df \
        .assign(average = lambda d: d.sense/(d.sense + d.antisense))\
        .assign(eb_estimate = lambda d: update_beta(alpha0, beta0, 
                                                d.sense, d.antisense, 
                                                param = 'success_rate')) \
        .assign(alpha1 = lambda d: update_beta(alpha0, beta0, 
                                                d.sense, d.antisense, 
                                                param = 'alpha'))\
        .assign(beta1 = lambda d: update_beta(alpha0, beta0, 
                                                d.sense, d.antisense, 
                                                param = 'beta')) 


vectorized_cdf = np.vectorize(beta.cdf)
def calulate_probability(posterior_df):
    ceil = 1e12
    ceil = 1e100
    return posterior_df \
        .assign(log_bf = lambda d: list(map(cal_beta_binom, d['NaOH: alpha1'], d['NaOH: beta1'], d['DNase: sense'], d['DNase: antisense']))) \
        .assign(pval = lambda d: vectorized_cdf(d['NaOH: eb_estimate'], d['DNase: alpha1'], d['DNase: beta1'] )) \
        .sort_values('pval')\
        .assign(qval = lambda d: d.pval.cumsum())


def fill_prior(d, priors):
    alpha0, beta0 = priors
    for col in d.columns[d.columns.str.contains('eb_estimate|alpha|beta')]:
        if 'eb' in col:
            d[col] = d[col].fillna(alpha0/(alpha0+beta0))
        elif 'alpha' in col:
            d[col] = d[col].fillna(alpha0)
        elif 'beta' in col:
            d[col] = d[col].fillna(beta0)
    return d


def plot_salmonTE(exp_df, treatment, ax):
    x = 'pvalue'
    plot_df = exp_df \
        .assign(pcolor = lambda d: np.where(d['padj']<0.05,'red', 'gray')) \
        .assign(alpha = lambda d: np.where(d[x]<0.05, 0.1,0.9))\
        .assign(log_padj = lambda d: -d[x].transform(lambda x: np.log2(x))) 
    for color, color_df in plot_df.groupby('pcolor'):
        alpha = 0.5 if color!="red" else 1
        color_df.plot.scatter('log2FoldChange','log_padj', 
                     ax = ax, alpha=alpha,
                     color = color)
    ax.set_title(treatment, size=15)
    
    xjust = 1
    ys = 0
    for i, row in plot_df\
                .query('pcolor=="red"') \
                .nlargest(5, 'log2FoldChange')\
                .iterrows():
        if xjust == 1:
            xs = 0.1
            ys = 0
            xjust = 2
        elif xjust == 2:
            xs = -0.1
            ys = 0
            xjust = 3
        elif xjust == 3:
            ys = 0.1
            xs = 0
            xjust = 1
        ax.text(row['log2FoldChange'] * (1+xs), row['log_padj'] * (1 + ys), 
                row['name'], color = 'red', fontsize=15)
    ax.set_xlabel(r'$log_2$ Fold change ($\frac{%s}{NaOH}$)' %treatment)
    ax.set_ylabel(r'$log_2$(p-value)')
    if treatment == "DNase":
        ax.set_xlim(plot_df.log2FoldChange.min() * 1.1, 
                   plot_df.log2FoldChange.max() * 1.4)
    ax.vlines(x = 0, ymin=0,ymax = 300, color = 'gray', alpha=0.5, linestyles=':')
    ax.set_ylim(0, plot_df.log_padj.max() * 1.1)
       
