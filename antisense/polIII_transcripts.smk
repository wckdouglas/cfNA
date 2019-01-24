import os
import glob
import re

wildcard_constraints:
    SAMPLENAME = '[A-Za-z0-9_]+',
    RNA_TYPE = '[A-Za-z]+'

PROJECT_PATH = '/stor/work/Lambowitz/cdw2854/cfNA'
FOLDER_PATH = PROJECT_PATH + '/tgirt_map'
SAMPLE_FOLDERS = glob.glob(FOLDER_PATH + '/Q*001')
SAMPLENAME_REGEX = '[Q][cC][fF].*001$'
SAMPLENAMES = map(os.path.basename, SAMPLE_FOLDERS)
SAMPLENAMES = filter(lambda x: re.search(SAMPLENAME_REGEX,x ), SAMPLENAMES)

TEMPLATE_FOLDER = FOLDER_PATH + '/{SAMPLENAME}'
INTERSECTED_FILE = TEMPLATE_FOLDER + '/count_temp/intersected/{RNA_TYPE}.dedup.bed.gz'
STRAND_COUNT_FILE = TEMPLATE_FOLDER + '/'

''' 
Pol III transcripts:
- Identification of RNA polymerase III-transcribed genes in eukaryotic genomes
tRNAs, 5S rRNA, U6 snRNA, SRP (7SL) RNA, RNase P and RNase MRP RNAs, vault RNAs, Y RNAs, 7SK RNA
'''
Pol_III_TRANSCRIPTS = ['TR[A-Z]-','5S_rRNA','RNU6','7SL','RMRP','RPPH1','Y-RNA','vaultRNA','7SK']
Pol_III_TRANSCRIPTS_REGEX = '|'.join(Pol_III_TRANSCRIPTS)

rule :