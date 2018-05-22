#!/usr/bin/env python

import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests
import os
import sys
import numpy as np

intable = sys.argv[1]

def fdr_wrapper(pval):
    cut, adj_p, a, b = multipletests(pval, alpha=0.05, method='fdr_bh')
    return cut

df = pd.read_table(intable, names=['chrom','start','end','name','score','strand',
                                     'fold','pval','padj','null2']) \
        .pipe(lambda d: d[fdr_wrapper(np.power(10,-d.pval))])

out_table = intable.replace('.bed','.filtered.bed')
df.to_csv(out_table, sep='\t', index=False, header=False)
os.system('bgzip -f %s' %out_table)
os.system('tabix -f %s.gz' %out_table)
