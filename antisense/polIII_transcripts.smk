import os
import glob
import re
import pandas as pd
import numpy as np
from multiprocessing import Pool

wildcard_constraints:
    SAMPLENAME = '[A-Za-z0-9_]+',
    RNA_TYPE = '[A-Za-z_]+',

PROJECT_PATH = '/stor/work/Lambowitz/cdw2854/cfNA'
FOLDER_PATH = PROJECT_PATH + '/tgirt_map'
SAMPLE_FOLDERS = glob.glob(FOLDER_PATH + '/Q*001')
SAMPLENAME_REGEX = '[Q][cC][fF].*001$'
SAMPLENAMES = map(os.path.basename, SAMPLE_FOLDERS)
SAMPLENAMES = filter(lambda x: re.search(SAMPLENAME_REGEX,x ), SAMPLENAMES)
SAMPLENAMES = filter(lambda x: re.search('Q[cC][fF][0-9]+|[nN][aA][0-9]+|[aA]ll[0-9]+',x), SAMPLENAMES)
SAMPLENAMES = list(SAMPLENAMES)

TEMPLATE_FOLDER = FOLDER_PATH + '/{SAMPLENAME}'
INTERSECTED_FOLDER = TEMPLATE_FOLDER + '/count_temp/intersected'
PolIII_FOLDER = INTERSECTED_FOLDER + '/PolIII'
INTERSECTED_FILE = INTERSECTED_FOLDER + '/{RNA_TYPE}.dedup.bed.gz'
PolIII_FILE = PolIII_FOLDER + '/{RNA_TYPE}.counts'
PolIII_COUNT = PolIII_FOLDER + '/{RNA_TYPE}.tsv'
STRAND_COUNT_FILE = TEMPLATE_FOLDER + '/'
RNA_TYPES = ['counts','rRNA_mt','small_RNA','sncRNA'] 
COUNT_TABLE = FOLDER_PATH + '/Counts/polIII_table.feather'
THREADS = 24

''' 
Pol III transcripts:
- Identification of RNA polymerase III-transcribed genes in eukaryotic genomes
tRNAs, 5S rRNA, U6 snRNA, SRP (7SL) RNA, RNase P and RNase MRP RNAs, vault RNAs, Y RNAs, 7SK RNA
'''
Pol_III_TRANSCRIPTS = ['TR[A-Z]-','5S_rRNA','RNU[1-6]','7SL','RMRP','RPPH1','Y-RNA','vaultRNA','7SK']
Pol_III_TRANSCRIPTS_REGEX = '|'.join(Pol_III_TRANSCRIPTS)


def read_polIII(tablename):
    samplename = tablename.split('/')[-5]
    df = pd.read_table(tablename,
                names = ['tid','read_strand',
                        'ref_strand','ttype','read_count'])\
        .assign(strandeness = lambda d: np.where(d.read_strand == d.ref_strand, 'Sense','Antisense'))\
        .pipe(pd.pivot_table,
                columns = 'strandeness',
                index = ['tid','ttype'],
                values = 'read_count',
                fill_value = 0)\
        .assign(samplename = samplename)\
        .reset_index()
    return df


rule all:
    input:
        COUNT_TABLE


rule count_polIII:
    input:
        polIII_FILE = expand(PolIII_FILE, 
                SAMPLENAME = SAMPLENAMES,
                RNA_TYPE = RNA_TYPES)

    output:
        TAB = COUNT_TABLE

    run:
        p = Pool(THREADS)
        dfs = p.map(read_polIII, input.polIII_FILE)
        p.close()
        p.join()
        df = pd.concat(dfs) \
            .reset_index(drop=True)
        df.to_feather(output.TAB)
        print('Written %s' %output.TAB)
        

rule extract_polIII:
    input:
        INTERSECTED_FILE

    params:
        TRANSCRIPTS = Pol_III_TRANSCRIPTS_REGEX,
        TEMP = PolIII_FOLDER 

    output:
        PolIII_FILE

    shell:
        'zcat {input} '\
        '| egrep --color=no "{params.TRANSCRIPTS}" '\
        "| awk ' ($3-$2) < (($9-$8)*1.1) {{print $10, $6, $12, $13 }}' OFS='\\t' "\
        '| sort  -T {params.TEMP} '\
        "| uniq -c "\
        "| awk '{{print $2,$3,$4,$5,$1}}' OFS='\\t' "\
        "> {output}"
        
