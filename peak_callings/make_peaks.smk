import glob
import os


MITO_INDEX= os.environ['REF'] + '/hg19/genome/chrM.fa'
ECOLI_INDEX= os.environ['REF'] + '/Ecoli/BL21_DE3.fa'
PROJECT_PATH= os.environ['WORK'] + '/cdw2854/cfNA/tgirt_map'
SAMPLE_NAMES = glob.glob(PROJECT_PATH + '/*001')
SAMPLE_NAMES = list(map(os.path.basename, SAMPLE_NAMES))
BED_PATH = PROJECT_PATH + '/bed_files'
MERGED_BED_PATH = BED_PATH + '/merged_bed'
COV_PATH = MERGED_BED_PATH + '/coverage'
STRANDED_BED_PATH = MERGED_BED_PATH + '/stranded'
MACS2_PATH = MERGED_BED_PATH + '/MACS2'
PEAK_READ_COUNT_PATH =MACS2_PATH + '/EV_count'
EV_COUNT_FILE = PEAK_READ_COUNT_PATH + '/{TREATMENT}.{FILTER}.ev_count.tsv'
GENOME_BAM = PROJECT_PATH + '/merged_bam/{TREATMENT}.bam'
SPLICED_TABLE = PROJECT_PATH + '/merged_bam/{TREATMENT}.spliced.tsv.gz'
SPLICED_EXON_TABLE = PROJECT_PATH + '/merged_bam/{TREATMENT}.spliced_exon.bed.gz'
ANNOTATED_PEAK_PATH = MACS2_PATH + '/annotated' 
ANNOTATED_PEAK = ANNOTATED_PEAK_PATH + '/{TREATMENT}.{FILTER}.tsv'
CONFIDENT_PEAK = MACS2_PATH + '/{TREATMENT}.{FILTER}.confident.tsv'
ANNOTATION_TABLE = os.environ['REF'] + '/hg19/new_genes/all_annotation.bed.gz' 
BED_TEMPLATE = BED_PATH + '/{SAMPLENAME}.bed.gz'
BAM_TEMPLATE = PROJECT_PATH + '/{SAMPLENAME}/Combined/primary.bam'
FILTER_BAM_TEMPLATE = PROJECT_PATH + '/{SAMPLENAME}/Combined/primary.chrM_filtered.bam'
MERGED_BED_TEMPLATE = MERGED_BED_PATH + '/{TREATMENT}.bed.gz'
STRANDED_BED_TEMPLATE = STRANDED_BED_PATH + '/{TREATMENT}.{FILTER}.{STRAND}.bed.gz'
MACS2_PEAK_TEMPLATE = MACS2_PATH + '/{TREATMENT}.{FILTER}.{STRAND}_peaks.narrowPeak' 
STRANDED_COV_FILE_TEMPLATE = COV_PATH + '/{TREATMENT}.{STRAND}.bigWig'
UNSTRANDED_COV_FILE_TEMPLATE = COV_PATH + '/{TREATMENT}.bigWig'
PEAK_FA = ANNOTATED_PEAK_PATH + '/{TREATMENT}.{RNA_TYPE}.fa'
CMSCAN_PEAK = ANNOTATED_PEAK_PATH + '/{TREATMENT}.{RNA_TYPE}.cmscan'
CMTBLOUT_PEAK = ANNOTATED_PEAK_PATH + '/{TREATMENT}.{RNA_TYPE}.tblout'
INTRON_TAB = ANNOTATED_PEAK_PATH + '/{TREATMENT}.intron.bed'
FOLD_FILE = ANNOTATED_PEAK_PATH + '/{TREATMENT}.{FILTER}.fold.fa'
GENOME = os.environ['REF'] + '/hg19_ref/genome/hg19_genome.fa'
EXON_ANNOTATION = os.environ['REF'] + '/hg19_ref/genes/exons.bed'
RNA_TYPES = ['Long_RNA','others']
FILTERS = ['filtered','unfiltered']
THREADS = 24

# set up treatments
TREATMENT = ['unfragmented','fragmented','polyA',
            'alkaline', 'all','exonuclease',
            'EV','RNP','EV-RNP',
            'MNase_EV','MNase_RNP','MNase_EV-RNP']
TREATMENT_REGEX = ['Q[Cc][Ff][0-9]+|Exo|[DE][DE]', 'Frag', 'L[12]', 
                'N[aA][0-9]+', 'All','Exo|[DE][DE]',
                '^MPF4','^MPF10','^MPCEV',
                '^PPF4','^PPF10','^PPCEV']

STRANDS = ['fwd', 'rvs']
TESTED_TREATMENT = ['unfragmented','all','MNase_EV','MNase_RNP','MNase_EV-RNP','EV','RNP','EV-RNP']
EV_LIBS = list(filter(lambda x: re.search('EV|RNP', x), TESTED_TREATMENT ))
regex_dict = {t:tr for t, tr in zip(TREATMENT, TREATMENT_REGEX)}
def get_bed(wildcards):
    regex = regex_dict[wildcards.TREATMENT]
    samplenames  = filter(lambda x: re.search(regex, x), SAMPLE_NAMES)
    beds = [BED_TEMPLATE.format(SAMPLENAME = x) for x in samplenames]
    return beds


def peak_filter(wildcards):
    TERM = wildcards.RNA_TYPE.replace('_',' ')
    if TERM == 'others':
        TERM = ''
    return '$pileup >= 4 && $sample_count >= 5 && $sense_gtype=="{TERM}"'.format(TERM = TERM)

SAMPLES = []
for T in TREATMENT:
    SAMPLES.extend(list(filter(lambda x: re.search(regex_dict[T], x), SAMPLE_NAMES)))
   

wildcard_constraints:
    RNA_TYPE = '[a-zA-Z_]+',
    FILTER = 'filtered|unfiltered',
    SAMPLENAME = 'Q.*|[MP]P.*',
    TREATMENT = '[-_A-Za-z]+',
    STRAND = 'rvs|fwd',


# Run commands
rule all:
    input:
        expand(CMTBLOUT_PEAK, TREATMENT = ['unfragmented'], RNA_TYPE = RNA_TYPES), 
        expand(CMSCAN_PEAK, TREATMENT = ['unfragmented'], RNA_TYPE = RNA_TYPES), 
        expand(STRANDED_COV_FILE_TEMPLATE, STRAND = STRANDS, TREATMENT = TESTED_TREATMENT),
        UNSTRANDED_COV_FILE_TEMPLATE.format(TREATMENT = 'alkaline'),
        expand(ANNOTATED_PEAK, 
                TREATMENT = ['unfragmented','EV','RNP','EV-RNP','MNase_EV',
                                            'MNase_RNP','MNase_EV-RNP'], 
                FILTER = FILTERS),
        expand(FOLD_FILE, TREATMENT = ['unfragmented'], FILTER = ['filtered']),
        expand(PEAK_FA, TREATMENT = ['unfragmented'], RNA_TYPE = ['Long_RNA','RBP']),
        expand(INTRON_TAB, TREATMENT = TESTED_TREATMENT),
        expand(EV_COUNT_FILE, TREATMENT = EV_LIBS, FILTER=FILTERS),

rule EV_counting:
    input:
        PEAK_FILE = CONFIDENT_PEAK.format(TREATMENT = 'unfragmented',
                                        FILTER = '{FILTER}'),
        BED = MERGED_BED_TEMPLATE
    
    output:
        EV_COUNT_FILE

    shell:
        'python EV_peaks.py {input.PEAK_FILE} {input.BED} '\
        '> {output}'

rule filter_peak:
    input:
        ANNOTATED_PEAK
    
    output:
        CONFIDENT_PEAK
    
    shell:
        'cat {input} '\
        "| csvtk filter2 -t -f '$sample_count >= 5 '"\
        '> {output} '


rule make_IGV:
    input:

    params:
        SESS = '/Users/wckdouglas/Desktop/plasma_bam/coverages/igv_session.xml'

    output:

    shell:
        'cat {input}'\
        '| sed 1d'\
        '| bedtools igv -i - -slop 100 '\
        '| -sess {params.SESS} -name -img png -path {output} '
        
rule intron:
    input:
        TAB = ANNOTATED_PEAK.replace('{FILTER}','unfiltered')

    params:
        INDEPENDENT_INTRON = '/stor/work/Lambowitz/ref/hg19/genome/independent_intron.bed'

    output:
        INTRON_TAB
        
    shell:
        'cat {input.TAB} '\
        '| csvtk cut -t -f chrom,start,end,peakname,pileup,strand ' \
        '| bedtools intersect -f 0.8 -F 0.8 -wb -a - -b {params.INDEPENDENT_INTRON} '\
        '> {output}'

rule fold:
    input:
        ANNOTATED_PEAK
    
    threads: THREADS
    params:
        GENOME = GENOME,
        OUTPREFIX = ANNOTATED_PEAK_PATH + '/unannotated_fold'

    output:
        FOLD_FILE

    shell:
        'cat {input} '\
        '| csvtk filter2 -t -f \'$sense_gname=="" && $pileup>3 && $sample_count > 5\'  '\
        '| csvtk cut -t -f1,2,3,4,5,13 '\
        '| sed 1d '\
        '| bedtools getfasta -fi {params.GENOME} -bed - -s '\
        '| RNAfold --jobs={threads} --noPS ' \
        "| cut -d' ' -f1 "\
        '> {output}'



rule macs2:
    #call strand specific peaks
    input:
        STRANDED_BED = STRANDED_BED_TEMPLATE,
        CONTROL_BED = MERGED_BED_TEMPLATE.format(TREATMENT = 'alkaline')

    params:
        OUT_PATH = MACS2_PATH,
        STRAND = '{STRAND}',
        TREATMENT = '{TREATMENT}',
        FILTER = '{FILTER}'
        
    output:
        MACS2_PEAK_TEMPLATE

    shell:
        'macs2 callpeak '\
        '--treatment {input.STRANDED_BED} ' \
        '--control {input.CONTROL_BED} '\
        '--outdir {params.OUT_PATH} ' \
        '--name {params.TREATMENT}.{params.FILTER}.{params.STRAND} '\
        '--nomodel  --format BEDPE --keep-dup all '\
        '--gsize hs --qvalue 0.05 '

       
rule split_strand:
    #splitting bed file into postive/negative strand bed file
    #also filter out full length exons
    input:
        BED = MERGED_BED_TEMPLATE,
        SPLICED_TABLE = SPLICED_EXON_TABLE.format(TREATMENT = 'unfragmented')
    
    params:
        OUT_PREFIX = STRANDED_BED_PATH + '/{TREATMENT}'

    output:
        expand(STRANDED_BED_TEMPLATE.replace('{TREATMENT}','{{TREATMENT}}'), STRAND = STRANDS, FILTER=FILTERS),
    
    shell:
        'python process_bed.py {input.BED} {params.OUT_PREFIX} {input.SPLICED_TABLE}'
    


rule merged_bed:
    #merging fragment files
    input:
        BEDS = lambda wildcards: get_bed(wildcards)

    params:
        TMP_DIR = MERGED_BED_TEMPLATE + '_TMP'

    output:
        MERGED_BED = MERGED_BED_TEMPLATE

    shell:
        'mkdir {params.TMP_DIR}'\
        ';zcat {input.BEDS} '\
        '| sort -k1,1 -k2,2n -k3,3n -T {params.TMP_DIR} '\
        '| bgzip > {output.MERGED_BED}'\
        '; tabix -p bed -f {output.MERGED_BED}'\
        '; rm -rf {params.TMP_DIR}'\


rule filter_bam:
    input:
        BAM = BAM_TEMPLATE

    params:
        MITO_INDEX = MITO_INDEX,
        ECOLI_INDEX = ECOLI_INDEX,

    output:
        BAM = FILTER_BAM_TEMPLATE

    shell:
        'cat {input.BAM} '\
        '| python ~/ngs_qc_plot/exogenous_filter.py '\
        ' -i - -o {output.BAM}  -x {params.MITO_INDEX} --nm 0.2 '


rule make_bed:
    input:
        BAM = lambda w: BAM_TEMPLATE.replace('{SAMPLENAME}',w.SAMPLENAME) \
                        if re.search('[aA]ll|[nN][aA]', w.SAMPLENAME) \
                        else FILTER_BAM_TEMPLATE.replace('{SAMPLENAME}',w.SAMPLENAME)

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
        '; tabix -p bed -f {output.BED} '\
	    '; rm -rf {params.TMP_FOLDER}'


COVERAGE_COMMAND = 'mkdir -p {params.TEMP_DIR}'\
        ';bedtools genomecov -bga -i {input.BED} -g {params.GENOME}'\
        '| sort -k1,1 -k2,2n -T {params.TEMP_DIR} '\
        '> {params.TEMP} '\
        '; bedGraphToBigWig {params.TEMP} {params.GENOME} {output.COV_FILE} '\
        '; rm -rf {params.TEMP} {params.TEMP_DIR}'

rule bed_coverage_strand:
    input:
        BED = STRANDED_BED_TEMPLATE.replace('{FILTER}','unfiltered')

    params:
        GENOME = GENOME + '.fai',
        TEMP = STRANDED_COV_FILE_TEMPLATE.replace('.bigWig','.bedGraph'),
        TEMP_DIR = STRANDED_COV_FILE_TEMPLATE + 'TMP'

    output:
        COV_FILE = STRANDED_COV_FILE_TEMPLATE
    
    shell:
        COVERAGE_COMMAND


rule bed_coverage_unstranded:
    input:
        BED = MERGED_BED_TEMPLATE
    
    params:
        GENOME = GENOME + '.fai',
        TEMP = UNSTRANDED_COV_FILE_TEMPLATE.replace('.bigWig','.bedGraph'),
        TEMP_DIR = UNSTRANDED_COV_FILE_TEMPLATE + 'TMP'

    output:
        COV_FILE = UNSTRANDED_COV_FILE_TEMPLATE
    
    shell:
        COVERAGE_COMMAND


rule PEAK_TO_FA:
    input:
        ANNOTATED_PEAK = ANNOTATED_PEAK.replace('{FILTER}','filtered')

    params:
        FILTER_TERM = lambda w: peak_filter(w),
        GENOME = GENOME

    output:
        FA = PEAK_FA

    shell:
        'cat {input.ANNOTATED_PEAK} '\
        "| csvtk filter2 -t -f '{params.FILTER_TERM}' "\
        '| csvtk cut -t -f chrom,start,end,peakname,score,strand '\
        '| sed 1d '\
        "| awk '{{printf \"%s\\t%s\\t%s\\t%s_%s:%s-%s\\t%s\\t%s\\n\", $1,$2-20,$3+20,$4,$1,$2,$3,$5,$6}}'"\
        '| bedtools getfasta -fi {params.GENOME} -bed - -s -name '\
        '| seqtk seq -U '\
        '> {output.FA}'

rule SCAN_FA:
    input:
        FA = PEAK_FA.replace('{RNA_TYPE}','{RNA_TYPE, [a-zA-Z_]+}')

    params:
        SCAN_REF = os.environ['REF'] + '/Rfam/Rfam.cm',
        THREADS = THREADS

    output:
        TBL_OUT = CMTBLOUT_PEAK,
        CMSCAN = CMSCAN_PEAK

    shell:
        'cmscan -o {output.CMSCAN} --tblout {output.TBL_OUT} '\
        ' --cpu {params.THREADS} '\
        '{params.SCAN_REF} {input.FA}'
        

rule peak_anntation:
    input:
        PEAK_FILES = expand(MACS2_PEAK_TEMPLATE\
                .replace('{TREATMENT}','{{TREATMENT}}')\
                .replace('{FILTER}','{{FILTER}}'), 
            STRAND = STRANDS),
        SPLICED_TABLE = SPLICED_TABLE.format(TREATMENT = 'unfragmented')
    
    params:
        ANNOTATION_TABLE = ANNOTATION_TABLE,
        BED_PATH = STRANDED_BED_PATH

    output:
        ANNOTATED_PEAK = ANNOTATED_PEAK

    shell:
        'python macs_peaks.py '\
        '{output.ANNOTATED_PEAK} {params.ANNOTATION_TABLE} '\
        '{params.BED_PATH} {input.SPLICED_TABLE} {input.PEAK_FILES} '


rule spliced_exon:
    input:
        SPLICED_TABLE

    params:
        EXON_ANNOTATION = EXON_ANNOTATION,
        TMP_DIR = os.path.dirname(SPLICED_EXON_TABLE)
 
    output:
        SPLICED_EXON_TABLE

    shell:
        'bedtools intersect -a {params.EXON_ANNOTATION} -b {input} -s '\
        '| sort -k1,1 -k2,2n -k3,3n -T {params.TMP_DIR} '\
        '| bgzip '\
        '> {output}'\
        '; tabix -f -p bed {output}'


rule find_splice:
    input:
        GENOME_BAM

    output:
        SPLICED_TABLE
    
    shell:
        'python fetch_junctions.py {input} {output}'
