#!/usr/bin/env python

from pandas import DataFrame
from numpy import where
import os
from shutil import copyfile
import glob

def run_TE_test(TE_path, experiment = 'untreated'):

    assert(experiment in ['untreated','dnase'])
    sample_regex = '[aA]ll' if experiment == "untreated" else  'Q[cC][fF][0-9]+'

    table_name = TE_path + '/phenotype.csv'
    sampleID = glob.glob(TE_path + '/Q*001')
    SampleID = map(os.path.basename, sampleID)
    DataFrame({'SampleID': list(SampleID)}) \
        .pipe(lambda d: d[d.SampleID.str.contains('N[aA]|%s' %sample_regex )])\
        .assign(phenotype = lambda d: where(d.SampleID.str.contains('N[aA]'), 'control', 'RNA')) \
        .sort_values('phenotype')\
        .to_csv(table_name,index=False)
    copyfile(table_name, table_name.replace('.csv',experiment+'.csv'))

    salmonTE = os.environ['WORK'] + '/cdw2854/src/SalmonTE/SalmonTE.py'

    command = '{salmonTE} test --inpath={input} '\
            '--outpath={output} '\
            '--tabletype=csv --figtype=png'\
            .format(salmonTE = salmonTE,
                    input = TE_path,
                    output = TE_path + '/test_' + experiment)
    print(command)
    os.system(command)


run_TE_test('untreated')
run_TE_test('dnase')
