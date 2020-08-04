#!/usr/bin/env python

import RNA
import random
from numpy import zeros
import dask.dataframe as dd
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('Folding test')

class FoldingTest():
    def __init__(self, n=1000):
        '''
        initiate test class,
        n = how many random sequence to generate
        '''
        self.n = n
        self.random_energies = zeros(self.n)
        self.energy = 0
        
    def testSeq(self, seq: str):
        '''
        testing sequence fold statstics significant by
        folding n random sequences and see how many times 
        we get a energy lower than the original folding energy
        
        input:
            seq: sequence to test
            
        output:
            p: pvalue of getting a fold energy lower than the original observed sequence
        '''
        seq = seq.upper()
        self.energy = self._GetFold(seq)
        self.random_energies = zeros(self.n)
        nucleotides = list(seq)
        for i in range(self.n):
            random.shuffle(nucleotides)
            _seq = ''.join(nucleotides)
            _energy = self._GetFold(_seq)
            self.random_energies[i] = _energy
        self.p_value = len(self.random_energies[self.random_energies <= self.energy])/self.n
        return self.p_value 
        
    def _GetFold(self, seq: str):
        fold, energy = RNA.fold(seq)
        return energy


if __name__ == '__main__':
    import pandas as pd
    import dask.dataframe as dd

    def test_fold(seq):
        logger.info('Testing %s' %seq)
        ft = FoldingTest()
        return ft.testSeq(seq)

    merged_peak = pd.read_excel('Sup_file_061620.xlsx', sheet_name = 'MACS2 peaks (Genomic)') \
        .pipe(dd.from_pandas, npartitions = 100) \
        .assign(p_val = lambda d: d['Peak sequence'].map(test_fold)) \
        .compute(num_workers=24, scheduler='processes')
    
    filename = 'folding_test.csv'
    merged_peak.filter(['ID','Peak sequence','p_val']).to_csv(filename, index=False)
    logger.info('Written %s' %filename)
