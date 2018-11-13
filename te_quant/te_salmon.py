import glob
import re
import os
from shutil import copyfile
from pandas import DataFrame
from numpy import where


PROJECT_PATH='/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
SAMPLE_FOLDERS = glob.glob(PROJECT_PATH + '/*001')
SAMPLE_FOLDERS = filter(lambda x: re.search('Q[cC][fF]', x), SAMPLE_FOLDERS)
SAMPLE_FOLDERS = filter(lambda x: re.search('[aA]ll|N[aA]|L[0-9E]+|[fF]rag|Phos|Qcf11',x), SAMPLE_FOLDERS)
SAMPLENAMES = map(os.path.basename, SAMPLE_FOLDERS)

SALMON_TE = '/stor/work/Lambowitz/cdw2854/src/SalmonTE/SalmonTE.py'
SALMON_TE_RESULT_PATH = PROJECT_PATH + '/salmonTE'
FQ_PATH = SALMON_TE_RESULT_PATH + '/fastq'
SAMPLE_FQ1 = FQ_PATH + '/{SAMPLENAME}_R1.fq.gz'
SAMPLE_FQ2 = FQ_PATH + '/{SAMPLENAME}_R2.fq.gz'
REPEAT_FQ1 = PROJECT_PATH + '/{SAMPLENAME}/repeats/repeats.1.fq.gz'
REPEAT_FQ2 = PROJECT_PATH + '/{SAMPLENAME}/repeats/repeats.2.fq.gz'
SALMON_TE_SAMPLE_RESULT_PATH = PROJECT_PATH + '/salmonTE/{SAMPLENAME}'
SALMON_DE_OUT = SALMON_TE_RESULT_PATH + '/test_{TREATMENT}/results.csv'
SALMON_PHENOTYPE_TABLE = SALMON_TE_RESULT_PATH + '/phenotype_{TREATMENT}.csv'
TREATMENTS = ['dnase', 'untreated']



def do_test():
    for experiment in TREATMENTS:
        sample_regex = '[aA]ll' if experiment == "untreated" else  'Q[cC][fF][0-9]+'

        table_name = SALMON_TE_RESULT_PATH + '/phenotype.csv'
        sampleID = glob.glob(SALMON_TE_RESULT_PATH + '/Q*' )
        sampleID = filter(lambda x: 'bam' not in x, sampleID)
        SampleID = list(map(os.path.basename, sampleID ))
        print(SampleID)
        DataFrame({'SampleID': SampleID} ) \
            .pipe(lambda d: d[d.SampleID.str.contains('N[aA]|%s' %sample_regex )])\
            .assign(phenotype = lambda d: where(d.SampleID.str.contains('N[aA]'), 'control', 'RNA' ) ) \
            .sort_values('phenotype')\
            .to_csv(table_name,index=False )
        copyfile(table_name, table_name.replace('.csv', '_'+experiment+'.csv'))


        command = '{salmonTE} test --inpath={input} '\
                '--outpath={output} '\
                '--tabletype=csv --figtype=png'\
                .format(salmonTE = SALMON_TE,
                        input = SALMON_TE_RESULT_PATH,
                        output = SALMON_TE_RESULT_PATH + '/test_' + experiment)
        print(command)
        os.system(command)


def TE_quant(FQs):

    command = SALMON_TE + \
        ' quant --reference=hs '\
        ' --outpath={SALMON_TE_OUT} '\
        ' --num_threads=12  {FQs}' \
        .format(SALMON_TE_OUT = SALMON_TE_RESULT_PATH,
                FQs = ' '.join(FQs))

    print(command)
    os.system(command)
        

def make_fq(SAMPLENAME):
    in_FQ1 = REPEAT_FQ1.format(SAMPLENAME = SAMPLENAME)
    in_FQ2 = REPEAT_FQ2.format(SAMPLENAME = SAMPLENAME)

    out_FQ1 = SAMPLE_FQ1.format(SAMPLENAME = SAMPLENAME.replace('_R1_001',''))
    out_FQ2 = SAMPLE_FQ2.format(SAMPLENAME = SAMPLENAME.replace('_R1_001',''))

    for _in, _out in zip([in_FQ1, in_FQ2], 
                         [out_FQ1, out_FQ2]):
        command = 'cp {in_f} {out_f}'.format(in_f = _in, out_f = _out)
        print(command)
        os.system(command)

    return '{} {}'.format(out_FQ1, out_FQ2)


#FQs = [make_fq(SAMPLENAME) for SAMPLENAME in SAMPLENAMES]
#TE_quant(FQs)
do_test()
