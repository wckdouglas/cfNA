import glob
import re
import os


PROJECT_PATH = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
SAMPLE_FOLDERS = glob.glob(PROJECT_PATH + '/*001')
SAMPLE_NAMES = map(os.path.basename, SAMPLE_FOLDERS)
SAMPLE_NAMES = filter(lambda x: re.search('Q[cC][fF]',x), SAMPLE_NAMES)
SAMPLE_NAMES = list(SAMPLE_NAMES)
SAMPLE_FOLDER = PROJECT_PATH + '/{SAMPLE_NAME}'
REPEAT_FOLDER = SAMPLE_FOLDER + '/repeats'
FQ1 = REPEAT_FOLDER + '/repeats.1.fq.gz'
FQ2 = REPEAT_FOLDER + '/repeats.2.fq.gz'
STR_BAM = REPEAT_FOLDER + '/STR.bam'
STR_BED = REPEAT_FOLDER + '/STR.bed'
REF_PATH = '/stor/work/Lambowitz/ref/hg19_ref/STR'
REF_INDEX = REF_PATH + '/STRdecoys.fasta'
THREADS = 24


rule all:
    input:
        expand(STR_BED, SAMPLE_NAME = SAMPLE_NAMES)
    

rule make_bed:
    input:
        BAM = STR_BAM
    
    output:
        BED = STR_BED
    
    shell:
        'samtools view -F4 -F2048 -bF256 {input.BAM} '\
        '| bam_to_bed.py -i -  -o {output.BED}'


rule mapping:
    input:
        FQ1 = FQ1, 
        FQ2 = FQ2,

    params:
        THREADS = THREADS,
        INDEX = REF_INDEX

    output:
        BAM = STR_BAM
    
    shell:
        'bowtie2  -p {params.THREADS} -x {params.INDEX} '\
        '--no-discordant --no-mixed --dovetail --mm '\
        ' -1 {input.FQ1} -2 {input.FQ2} '\
        '| samtools view -b@ {params.THREADS} - '\
        '> {output.BAM}'
   