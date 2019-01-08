
#!/usr/bin/env python

from matplotlib import use as mpl_use
mpl_use('agg')
import pandas as pd
import numpy as np
import re
import sys
from repeats_utils import get_repeat_df
import pymc3 as pm
from functools import partial
import logging
from multiprocessing import Pool
from scipy.stats import mode
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
logger = logging.getLogger(__name__)
EPSILON=1e12


def ab_test(obs, return_p, dnase_sense, dnase_antisense, naoh_sense, naoh_antisense):
    '''
    Modeling empirical beta distribytion and use that as prior
    add new evidence (NaOH or DNase)
    sample delta strand
    calulate bayes factor for DNase > NaOH by at least 10%
    '''
    with pm.Model() as model:
    
        # fit beta binom
        naoh_empirical_alpha = pm.Exponential('alpha', 1)
        naoh_empirical_beta = pm.Exponential('beta', 1)
        beta_binom_prior = pm.Beta('beta_prior',  naoh_empirical_alpha, naoh_empirical_beta, observed = obs)
        alpha = pm.Normal('alpha1', mu = naoh_empirical_alpha, sd = 1)
        beta = pm.Normal('beta1', mu = naoh_empirical_beta, sd = 1)
        
        #inference
        dnase_sense = pm.Beta('dnase_sense', 
                          alpha = alpha + dnase_sense,  
                          beta = beta + dnase_antisense)
        naoh_sense = pm.Beta('naoh_sense', 
                          alpha = alpha + naoh_sense,  
                          beta = beta + naoh_antisense)
        diff = pm.Deterministic('delta', dnase_sense - naoh_sense)
        step = pm.NUTS()
        
        progressbar = not return_p
        trace = pm.sample(1000, step, tune=1000,progressbar=True, cores=24)
        
    if return_p:
        delta = trace['delta']
        h1 = np.sum(delta >= 0.1)
        h0 = np.sum(delta <= 0)

        p_h0 = h0/len(delta)
        p_h1 = h1/len(delta)
        bf = p_h1/p_h0 if p_h0 > 0 else EPSILON
        return bf, delta.mean()
    else:
        return trace

    
def construct_empirical_bayes(naoh_df, coverage_threshold):
    # set up Empirical bayes
    high_counts = naoh_df\
        .query('sense + antisense > %i' %coverage_threshold)
    obs = high_counts.sense/(high_counts.sense+high_counts.antisense)
    logger.info('Using %i samples for constructing empirical bayes' %len(obs))
    return obs


def main():
    # read in count dataframe for sample
    count_path = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/Counts/all_counts'
    count_table = count_path + '/spreaded_all_counts.tsv'
    simple_repeat_tab = count_path + '/simple_repeats.feather'

    dnase_sample_regex = 'Q[cC][fF][0-9]+'
    naoh_sample_regex = 'Q[cC][fF]_[Nn][aA]'
    coverage_threshold = 100

    df = pd.read_table(count_table)\
        .query('grouped_type == "Repeats"')\
        .reset_index(drop=True)

    # merge reverse complement and split treatment
    dnase_df = get_repeat_df(df, sample_regex = dnase_sample_regex)
    naoh_df = get_repeat_df(df, sample_regex = naoh_sample_regex)
    obs = construct_empirical_bayes(naoh_df, coverage_threshold)


    # make data frame
    sample_df = pd.concat([dnase_df.assign(treatment = 'DNase'),
                         naoh_df.assign(treatment='NaOH')])\
        .pipe(pd.melt, id_vars = ['treatment','gene_name','gene_id']) \
        .assign(variable = lambda d: d.treatment + ': ' + d.variable) \
        .pipe(pd.pivot_table, index=['gene_name', 'gene_id'],
                columns = 'variable', values = 'value', fill_value = 0) \
        .reset_index()
    logger.info('DataFrame ready!')

    # do test
    pval_func = partial(ab_test, obs, True)
    res = map(pval_func,
            sample_df['DNase: sense'],
            sample_df['DNase: antisense'],
            sample_df['NaOH: sense'],
            sample_df['NaOH: antisense'])

    bf, delta_sense = zip(*res)
    sample_df['bayes_factor'] = bf
    sample_df['delta'] = delta_sense
    sample_df.to_feather(simple_repeat_tab)
    logger.info('Written %s' %simple_repeat_tab)


if __name__ == '__main__':
    main()