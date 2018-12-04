
import os
import pandas as pd

PLATELETS_PATH = '/stor/work/Lambowitz/cdw2854/cfNA/platelets'
TISSUE_PATH = PLATELETS_PATH + '/tissues'
GENE_PATH = os.environ['REF'] + '/hg19_ref/genes'
KALLISTO_INDEX = GENE_PATH + '/genes_protein.kallisto_idx'

# VARIABLES
PLATELETS_DATA_SRR =  ['SRR590742' + str(i) for i in range(3,9)]
PLATELET_FQ = PLATELETS_PATH + '/data/{SRR}.fastq.gz'
UNZIPPED_PLATELET_FQ = PLATELET_FQ.replace('.gz','')
KALLISTO_PATH = PLATELETS_PATH + '/kallisto'
KALLISTO_RESULT_DIR = KALLISTO_PATH + '/{SRR}' 
TPM_TABLE = KALLISTO_RESULT_DIR + '/abundance.tsv'
COMBINED_TPM = KALLISTO_PATH + '/platelets_tpm.feather'
RNA_TABLE = TISSUE_PATH + '/rna_{RNA_SOURCE}.tsv'
RNA_SOURCES = ['celline','tissue']
COMBINED_RNA = TISSUE_PATH + '/published.feather'


def gene_annot():
    return pd.read_table(GENE_PATH + '/genes.annot',
                names = ['target_id','gid','gtype'])\
        .assign(gid = lambda d: d.gid.str.split('.', expand=True).iloc[:,0])


rule all:
    input:
        COMBINED_TPM,
        COMBINED_RNA


rule combine_tissue:
    input:
        TABS = expand(RNA_TABLE, RNA_SOURCE = RNA_SOURCES)
    
    output:
        TAB = COMBINED_RNA

    run:
        dfs = []
        for TAB in TABS:
            df = pd.read_table(TAB)
            dfs.append(df)

        pd.concat(dfs,sort=True)\
            .reset_index(drop=True)\
            .to_feather(output.TAB)


rule make_tissue_table:
    input:
    
    params:
        RNA_SOURCE = '{RNA_SOURCE}'

    output:
        TAB = RNA_TABLE,

    shell:
        'curl https://www.proteinatlas.org/download/rna_{params.RNA_SOURCE}.tsv.zip '\
        '> {output.TAB}.zip'\
        '; unzip -p {output.TAB}.zip > {output.TAB}'


rule make_tpm_table:
    input:
        TPMS = expand(TPM_TABLE, SRR = PLATELETS_DATA_SRR)
    
    output:
        TAB = COMBINED_TPM
    
    run:
        df = []
        for TPM in input.TPMS:
            df.append(pd.read_table(TPM, usecols = [0,4])\
                .merge(gene_annot())\
                .groupby('gid', as_index=False)\
                .agg({'tpm':'sum'}))

        pd.concat(df)\
            .groupby(['gid'], as_index=False)\
            .mean() \
            .rename(columns = {'gid':'Gene'}) \
            .to_feather(output.TAB)


rule kallisto:
    input:
        FQ = PLATELET_FQ

    params:
        INDEX = KALLISTO_INDEX,
        OUT = KALLISTO_RESULT_DIR,
        THREADS = 12
    
    output:
        FILE = TPM_TABLE

    shell:
        'kallisto quant --bias '\
        '--output-dir {params.OUT} '\
        '--index {params.INDEX} '\
        '--threads {params.THREADS} ' \
        '-l 100 -s 20 --single '\
        '{input.FQ}'


rule zip_fastq:
    input:
        FQ = UNZIPPED_PLATELET_FQ

    output:
        ZIPPED_FQ = PLATELET_FQ

    shell:
        'gzip {input.FQ}'


rule download_fastq:
    input:
    
    params:
        SRR = lambda w: w.SRR,
        TMP = PLATELETS_PATH + '/data'

    output:
        FQ = UNZIPPED_PLATELET_FQ
    
    shell:
        'fasterq-dump {params.SRR} --temp {params.TMP} -O {params.TMP}'

    

    
