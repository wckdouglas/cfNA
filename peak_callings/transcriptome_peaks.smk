import glob
import os
import re

wildcard_constraints:
    STRAND = 'rvs|fwd',
    SAMPLE_NAME = '[A-Za-z0-9_\-]+',
    TREATMENT = '[A-Za-z0-9_\-]+'

REF_PATH = '/stor/work/Lambowitz/ref/hg19_ref'
BED12 = REF_PATH + '/genes/genes.bed12'
TRANSCRIPTOME_FA = REF_PATH + '/genes/transcriptome.fa'
TRANSCRIPTOME_FA_IDX = TRANSCRIPTOME_FA + '.1.bt2'
JUNCTION_SITE = TRANSCRIPTOME_FA.replace('.fa','.junction.bed')
PROJECT_PATH= '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
TRANSCRIPTOME_PROJECT_PATH = PROJECT_PATH + '/transcriptome'
COMBINED_BED = TRANSCRIPTOME_PROJECT_PATH + '/{TREATMENT}.bed.gz'
COMBINED_BED_IDX = COMBINED_BED + '.tbi'
STRANDED_BED = TRANSCRIPTOME_PROJECT_PATH + '/{TREATMENT}.{STRAND}.bed.gz'
STRANDED_COV = TRANSCRIPTOME_PROJECT_PATH + '/coverage/{TREATMENT}.{STRAND}.bigWig'
COV = TRANSCRIPTOME_PROJECT_PATH + '/coverage/{TREATMENT}.bigWig'
MACS2_PATH = TRANSCRIPTOME_PROJECT_PATH + '/macs2'
MACS_PEAK = MACS2_PATH + '/{TREATMENT}.{STRAND}_peaks.narrowPeak'
GENOMICS_MACS_PEAK = MACS2_PATH + '/{TREATMENT}.{STRAND}_peaks_genomics.narrowPeak.gz'
GENOMICS_MACS_PEAK_TBI = MACS2_PATH + '/{TREATMENT}.{STRAND}_peaks_genomics.narrowPeak.gz.tbi'
INTERSECTED = MACS2_PATH + '/{TREATMENT}.{STRAND}.bed'
BAM = PROJECT_PATH + '/{SAMPLE_NAME}/Combined/primary_no_sncRNA_tRNA_rRNA_repeats.bam'
FQ1_TEMPLATE = BAM.replace('.bam','.1.fq.gz')
FQ2_TEMPLATE = FQ1_TEMPLATE.replace('.1.fq.gz','.2.fq.gz')
SAMPLE_NAMES, = glob_wildcards(BAM)
SAMPLE_NAMES = filter(lambda x: not re.search('L[0-9E]+',x), SAMPLE_NAMES)
SAMPLE_NAMES = list(SAMPLE_NAMES)
OUT_FOLDER = PROJECT_PATH + '/{SAMPLE_NAME}/transcriptome'
BAM_FILE = OUT_FOLDER + '/aligned.bam'
BED_FILE = OUT_FOLDER + '/aligned.bed.gz'
BED_IDX = OUT_FOLDER + '/aligned.bed.gz.tbi'
TREATMENTS = ['unfragmented','fragmented','polyA',
            'alkaline', 'all','exonuclease',
            'EV','RNP','EV-RNP',
            'MNase_EV','MNase_RNP','MNase_EV-RNP','high_salt']
TREATMENT_REGEXES = ['Q[Cc][Ff][0-9]+|Exo|[DE][DE]', 'Frag', 'L[12]', 
                'N[aA][0-9]+', 'All','Exo|[DE][DE]',
                '^MPF4','^MPF10','^MPCEV',
                '^PPF4','^PPF10','^PPCEV','[Hh][sS][0-9]+']
STRANDS = ['fwd', 'rvs']
regex_dict = {t:tr for t, tr in zip(TREATMENTS, TREATMENT_REGEXES)}
def get_bed(wildcards):
    regex = regex_dict[wildcards.TREATMENT]
    samplenames  = filter(lambda x: re.search(regex, x), SAMPLE_NAMES)
    beds = [BED_FILE.format(SAMPLE_NAME = x) for x in samplenames]
    return beds



COVERAGE_COMMAND = 'mkdir -p {params.TEMP_DIR}'\
        ';bedtools genomecov -bga -i {input.BED} -g {params.GENOME}.fai '\
        '| sort -k1,1 -k2,2n -T {params.TEMP_DIR} '\
        '> {params.TEMP} '\
        '; bedGraphToBigWig {params.TEMP} {params.GENOME}.fai {output.COV_FILE} '\
        '; rm -rf {params.TEMP} {params.TEMP_DIR}'


rule all:
    input:
        expand(INTERSECTED, 
            TREATMENT = 'unfragmented', 
            STRAND = STRANDS),
        expand(GENOMICS_MACS_PEAK_TBI, 
            TREATMENT = 'unfragmented', 
            STRAND = 'fwd'),
        JUNCTION_SITE,
        expand(COV, TREATMENT = 'alkaline'),
        expand(STRANDED_COV, TREATMENT = 'unfragmented', STRAND = STRANDS)
    
rule index_peaks:
    input:
        GENOMICS_MACS_PEAK

    output:
        GENOMICS_MACS_PEAK_TBI

    shell:
        'tabix -p bed -f {input}'


rule genomic_peaks:
    input:
        MACS_PEAK

    params:
        BED12 = BED12

    output:
        GENOMICS_MACS_PEAK

    shell:
        ' cat {input} '\
        '| python project_transcriptome.py -i - -o - -b {params.BED12} '\
        '| sort -k1,1 -k2,2n -k3,3n '\
        '| bgzip '\
        '> {output}'


rule annotate_peak:
    input:
        PEAK = MACS_PEAK,
        BED = STRANDED_BED

    output:
        INTERSECTED

    shell:
        'bedtools coverage -a {input.PEAK} -b {input.BED} -counts '\
        '> {output}'
        

rule coverage:
    input:
        BED = COMBINED_BED

    params:
        TEMP_DIR = COV + '_temp',
        TEMP = COV + '.bedGraph',
        GENOME = TRANSCRIPTOME_FA

    output:
        COV_FILE = COV

    shell:
        COVERAGE_COMMAND
        
rule coverage_stranded:
    input:
        BED = STRANDED_BED

    params:
        TEMP_DIR = STRANDED_COV + '_temp',
        TEMP = STRANDED_COV + '.bedGraph',
        GENOME = TRANSCRIPTOME_FA

    output:
        COV_FILE = STRANDED_COV

    shell:
        COVERAGE_COMMAND


rule macs2:
    #call strand specific peaks
    input:
        STRANDED_BED = STRANDED_BED,
        CONTROL_BED = COMBINED_BED.format(TREATMENT = 'alkaline')

    params:
        OUT_PATH = MACS2_PATH,
        STRAND = '{STRAND}',
        TREATMENT = '{TREATMENT}'
        
    output:
        MACS_PEAK

    shell:
        'macs2 callpeak '\
        '--treatment {input.STRANDED_BED} ' \
        '--outdir {params.OUT_PATH} ' \
        '--name {params.TREATMENT}.{params.STRAND} '\
        '--nomodel  --format BEDPE --keep-dup all '\
        ' --qvalue 0.05 '
  #      '--control {input.CONTROL_BED} '\




rule stranded_bed:
    input:
        BED = COMBINED_BED,
        IDX = COMBINED_BED_IDX
    
    params:
        STRAND = lambda w: '+' if w.STRAND=="fwd" else '-',
        TMP = STRANDED_BED + '_tmp'

    output:
        STRANDED_BED
    
    shell:
        'mkdir -p {params.TMP} '\
        '; zcat {input.BED} '\
        "| awk '$6 == \"{params.STRAND}\" && $3-$2 < 1000' "\
        '| sort -k1,1 -k2,2n -k3,3n -T {params.TMP} '\
        '| bgzip '\
        '> {output} '  \
        '; rm -rf {params.TMP}'


rule index_combined_bed:
    input:
        COMBINED_BED
    
    output:
        COMBINED_BED_IDX
    
    shell:
        'tabix -p bed -f {input}'
        

rule combine_bed:
    input:
        BEDS = lambda w: get_bed(w)

    params:
        TMP = COMBINED_BED + '_tmp'

    output:
        COMBINED_BED
    
    shell:
        'mkdir -p {params.TMP} '\
        '; zcat {input.BEDS} '\
        '| sort -k1,1 -k2,2n -k3,3n -T {params.TMP} '\
        '| bgzip > {output} '\
        '; rm -rf {params.TMP}'


rule index_sample_bed:
    input:
        BED_FILE

    output:
        BED_IDX

    shell:
        'tabix -p bed -f {input}'


rule bam_to_bed:
    input:
        BAM_FILE

    params:
        ID = '{SAMPLE_NAME}',
        TMP_FOLDER = BED_FILE + '_tmp'

    output:
        BED_FILE

    shell:
        'mkdir -p {params.TMP_FOLDER}'\
        '; samtools view -b -F 256 -F 2048 -F4 {input} '\
        '| bam_to_bed.py -i - -o - --add_cigar '\
        '| sort -k1,1 -k2,2n -k3,3n -k6,6 --temporary-directory={params.TMP_FOLDER}'\
        "| deduplicate_bed.py --infile - --outfile - --threshold 1 -d '_' --ct 6 "\
        '| poisson_umi_adjustment.py -i - -o - --umi 6 --prefix {params.ID} '\
        '| sort -k1,1 -k2,2n -k3,3n --temporary-directory={params.TMP_FOLDER} '\
        '| bgzip '\
        '> {output}'\
        '; rm -rf {params.TMP_FOLDER}'


rule bowtie_mapping:
    input:
        IDX = TRANSCRIPTOME_FA_IDX,
        FQ1 = FQ1_TEMPLATE,
        FQ2 = FQ2_TEMPLATE,

    params:
        IDX = TRANSCRIPTOME_FA
    output:
        BAM_FILE

    shell:
        'bowtie2 '\
        '--mm  --local -N 1 -D 20 -L 20 -X 1000 --no-mixed --no-discordant '\
        ' -x {params.IDX} -1 {input.FQ1} -2 {input.FQ2} '\
        '| samtools view -b '\
        '> {output}'
        
rule make_fq:
    input:
        BAM 
    
    output:
        FQ1 = FQ1_TEMPLATE, 
        FQ2 = FQ2_TEMPLATE
    
    shell:
        'samtools fastq -1 {output.FQ1} -2 {output.FQ2} {input}'

rule junction_ref:
    input:

    params:
        GTF = REF_PATH + '/genes/genes.gtf'

    output:
        JUNCTION_SITE
    
    shell:
        'gtfToGenePred {params.GTF} /dev/stdout '\
        '| python genePred_to_transcriptome_bed.py '\
        '> {output}'
    


rule index_transcript:
    input:
        TRANSCRIPTOME_FA

    output:
        TRANSCRIPTOME_FA_IDX

    shell:
        'bowtie2-build {input} {input}'

rule make_transcript:
    input:

    params:
        GENOME = REF_PATH + '/genome/hg19_genome.fa',
        GTF = REF_PATH + '/genes/genes.gtf'

    output:
        TRANSCRIPTOME_FA

    shell:
        'gffread -w {output} -g {params.GENOME} {params.GTF}'
