import glob

wildcard_constraints:
    TREATMENT = '[A-Za-z]+'

BAM_PATH = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/merged_bam/small_rna'
SUFFIX ='.smallRNA.total.nameSorted.bam'
BAMS = glob.glob(BAM_PATH + '/*' + SUFFIX)
TREATMENTS = list(map(lambda x: os.path.basename(x).split('.')[0], BAMS))
BAM_TEMPLATE = BAM_PATH + '/{TREATMENT}' + SUFFIX
FRAMGNET_TEMPLATE = BAM_TEMPLATE.replace('.bam','.tsv.gz')
INDEX_TEMPLATE = FRAMGNET_TEMPLATE.replace('.gz','.gz.tbi')

rule all:
    input:
        expand(INDEX_TEMPLATE, TREATMENT = TREATMENTS)


rule index:
    input:
        FRAMGNET_TEMPLATE

    output:
        INDEX_TEMPLATE
    
    shell:
        'tabix -fp bed {input}'


rule extract:
    input:
        BAM_TEMPLATE

    params:
        TEMP = BAM_PATH

    output:
        FRAMGNET_TEMPLATE

    shell:
        'cat {input} '\
        '| filter_soft_clip.py -i - -o - --pe  '\
        '| bam_to_bed.py -i - '\
        "| grep '^TR' --color=no "\
        "| awk '$6==\"+\"' "\
        '| cut -f1,2,3  '
        '| sort -k1,1 -k2,2n -k3,3n -T {params.TEMP}'\
        '| uniq -c '\
        "| awk '{{print $2,$3,$4,$1}}' OFS='\\t'" \
        '| bgzip ' \
        '> {output}'
