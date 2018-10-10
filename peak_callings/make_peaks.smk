import glob
import os


PROJECT_PATH= os.environ['WORK'] + '/cdw2854/cfNA/tgirt_map'
SAMPLE_NAMES = glob.glob(PROJECT_PATH + '/Q*001')
SAMPLE_NAMES = list(map(os.path.basename, SAMPLE_NAMES))
BED_PATH = PROJECT_PATH + '/bed_files'
MERGED_BED_PATH = BED_PATH + '/merged_bed'
COV_PATH = MERGED_BED_PATH + '/coverage'
STRANDED_BED_PATH = MERGED_BED_PATH + '/stranded'
MACS2_PATH = MERGED_BED_PATH + '/MACS2'
ANNOTATED_PEAK_PATH = MACS2_PATH + '/annotated' 
BED_TEMPLATE = BED_PATH + '/{SAMPLENAME}.bed.gz'
BAM_TEMPLATE = PROJECT_PATH + '/{SAMPLENAME}/Combined/primary.bam'
MERGED_BED_TEMPLATE = MERGED_BED_PATH + '/{TREATMENT}.bed.gz'
STRANDED_BED_TEMPLATE = STRANDED_BED_PATH + '/{TREATMENT}.{STRAND}.bed.gz'
MACS2_PEAK_TEMPLATE = MACS2_PATH + '/{TREATMENT}.{STRAND}_peaks.narrowPeak' 
STRANDED_COV_FILE_TEMPLATE = COV_PATH + '/{TREATMENT}.{STRAND}.bigWig'
UNSTRANDED_COV_FILE_TEMPLATE = COV_PATH + '/{TREATMENT}.bigWig'

# set up treatments
TREATMENT = ['unfragmented','fragmented','polyA',
            'alkaline', 'all','exonuclease']
TREATMENT_REGEX = ['Q[Cc][Ff][0-9]+|Exo|[DE][DE]', 'Frag', 'L[12]', 
                'N[aA][0-9]+', 'All','Exo|[DE][DE]']
STRANDS = ['fwd', 'rvs']
regex_dict = {t:tr for t, tr in zip(TREATMENT, TREATMENT_REGEX)}
def get_bed(wildcards):
    regex = regex_dict[wildcards.TREATMENT]
    samplenames  = filter(lambda x: re.search(regex, x), SAMPLE_NAMES)
    beds = [BED_TEMPLATE.format(SAMPLENAME = x) for x in samplenames]
    return beds

SAMPLES = []
for T in TREATMENT:
    SAMPLES.extend(list(filter(lambda x: re.search(regex_dict[T], x), SAMPLE_NAMES)))

# set up contstraints
wildcard_constraints:
    SAMPLENAME = 'Q.*',
    TREATMENT = '|'.join(TREATMENT),
#    STRAND = '^rvs$|^fwd$',

# Run commands
rule all:
    input:
#        expand(BED_TEMPLATE, SAMPLENAME =  SAMPLE_NAMES),
#        expand(MERGED_BED_TEMPLATE, TREATMENT = ['alkaline', 'unfragmented', 'all']),
#        expand(STRANDED_BED_TEMPLATE, 
#                TREATMENT = ['unfragmented','all'], 
#                STRAND = STRANDS),
        expand(MACS2_PEAK_TEMPLATE, 
                TREATMENT = ['unfragmented','all'], 
                STRAND = STRANDS),
        
rule macs2:
    #call strand specific peaks
    input:
        STRANDED_BED = STRANDED_BED_TEMPLATE,
        CONTROL_BED = MERGED_BED_TEMPLATE.format(TREATMENT = 'alkaline')

    params:
        OUT_PATH = MACS2_PATH,
        STRAND = '{STRAND}',
        TREATMENT = '{TREATMENT}'
        
    output:
        MACS2_PEAK_TEMPLATE

    shell:
        'macs2 callpeak '\
        '--treatment {input.STRANDED_BED} ' \
        '--control {input.CONTROL_BED} '\
        '--outdir {params.OUT_PATH} ' \
        '--name {params.TREATMENT}.{params.STRAND} '\
        '--nomodel  --format BEDPE --keep-dup all '\
        '--gsize hs --qvalue 0.05 '
    

rule split_strand:
    #splitting bed file into postive/negative strand bed file
    #also filter out full length exons
    input:
        BED = MERGED_BED_TEMPLATE
    
    params:
        OUT_PREFIX = STRANDED_BED_PATH + '/{TREATMENT}'

    output:
        expand(STRANDED_BED_TEMPLATE.replace('{TREATMENT}','{{TREATMENT}}'), STRAND = STRANDS)
    
    shell:
        'python process_bed.py {input.BED} {params.OUT_PREFIX}'
    

rule merged_bed:
    #merging fragment files
    input:
        BEDS = lambda wildcards: get_bed(wildcards)

    output:
        MERGED_BED = MERGED_BED_TEMPLATE

    shell:
        'zcat {input.BEDS} '\
        '| sort -k1,1 -k2,2n -k3,3n '\
        '| bgzip '\
        '> {output.MERGED_BED}'\
        '; tabix -p bed -f {output.MERGED_BED}'



rule make_bed:
    input:
        BAM = BAM_TEMPLATE

    params:
        TMP_FOLDER = BED_PATH + '/{SAMPLENAME}_TMP',
        WHITELIST = os.environ['REF'] + '/hg19/genome/wgEncodeDacMapabilityConsensusExcludable.bed.gz',
        SAMPLENAME = '{SAMPLENAME}'

    output:
        BED = BED_TEMPLATE
    
    shell:
        'mkdir -p {params.TMP_FOLDER}; ' \
        'cat {input.BAM} ' \
	    '| bam_to_bed.py --in_bam - --add_cigar '\
        '| awk \'$7!~"N"\' ' \
	    '| sort -k1,1 -k2,2n -k3,3n -k6,6 --temporary-directory={params.TMP_FOLDER} '\
	    '| bedtools intersect -v -a - '\
        ' -b {params.WHITELIST} '\
        "| deduplicate_bed.py --infile - --outfile - --threshold 1 -d '_' --ct 6 " \
        "| poisson_umi_adjustment.py -i - -o - --umi 6 --prefix {params.SAMPLENAME} " \
        "| sort -k1,1 -k2,2n -k3,3n --temporary-directory={params.TMP_FOLDER} "\
        '| bgzip > {output.BED} '\
        '; tabix -f {output.BED} '\
	    '; rm -rf {params.TMP_FOLDER}'


COVERAGE_COMMAND = 'bedtools genomecov -bga -i {input.MERGED_BED} -g {params.GENOME}'\
        '| sort -k1,1 -k2,2n '\
        '> {params.TEMP} '\
        '; bedGraphToBigWig {params.TEMP} {params.GENOME} {output.COV_FILE} '\
        '; rm {params.TEMP} '
rule bed_coverage_strand:
    input:
        BED = STRANDED_BED_TEMPLATE

    params:
        GENOME = os.environ['REF'] + '/hg19/genome/hg19_genome.fa.fai'
        TEMP = STRANDED_COV_FILE_TEMPLATE.replace('bigWig','.bedGraph')

    output:
        COV_FILE = STRANDED_COV_FILE_TEMPLATE
    
    shell:
        COVERAGE_COMMAND


rule bed_coverage_unstranded:
    input:
        BED = MERGED_BED_TEMPLATE
    
    params:
        GENOME = os.environ['REF'] + '/hg19/genome/hg19_genome.fa.fai'
        TEMP = UNSTRANDED_COV_FILE_TEMPLATE.replace('bigWig','.bedGraph')

    output:
        COV_FILE = UNSTRANDED_COV_FILE_TEMPLATE
    
    shell:
        COVERAGE_COMMAND
