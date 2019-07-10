#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys
import re

seqs = []
for line in sys.stdin:
    if not line.startswith('>'):
        seq = line.strip()
        if len(seq.replace('-','')):
            seqs.append(list(seq))


seqs = np.array(seqs)
rows, cols = seqs.shape

concensus_seq = ''
for pos in range(cols):
    bases = seqs[:,pos]
    #bases = bases[bases!='-']
    b, bcount = np.unique(bases, return_counts=True)
    cb = b[bcount.argmax()]
    #if cb == '-':
    #    cb = str(b[bcount.argmax()][0])
    #elif len(b) > 1:
    #    cb = str(b[bcount == np.sort(bcount)[-2]][0])

    concensus_seq += cb
cs = concensus_seq.replace('-','')
cs = re.sub('ACA$','CCA', cs)
print(cs)





