merged_bam = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/merged_bam/{TREATMENT}.nameSorted.bam'
kallisto_bam = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/kallisto_all_result/bam_files/{TREATMENT}_kallisto.bam'
sncRNA = '/stor/work/Lambowitz/ref/hg19_ref/genes/sncRNA_x_protein.bed'
out_path = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/merged_bam/viz_bam'
sncRNA_bam = out_path + '/{TREATMENT}_sncRNA.bam'
out_bam = out_path + '/{TREATMENT}_viz.bam'
TREATMENTS = ['unfragmented','polyA',
            'EV','RNP','RNP-EV',
            'MNase_EV','MNase_RNP','MNase_EV-RNP','high_salt']


rule all:
    input:
        expand(out_bam, TREATMENT = TREATMENTS)

rule merge:
    input:
        sncBAM = sncRNA_bam,
        kallisto = kallisto_bam

    output:
        BAM = out_bam.replace('{TREATMENT}','{TREATMENT,[a-zA-Z-_]+}')


    shell:
        'samtools merge - {input.sncBAM} {input.kallisto} '\
        '| sambamba sort -t 24 -o {output.BAM} /dev/stdin '

rule make_smRNA:
    input:
        BAM = merged_bam

    params:
        sncRNA = sncRNA

    output:
        BAM = sncRNA_bam

    shell:
        'bedtools pairtobed -abam {input.BAM} '\
        ' -b {params.sncRNA} -type both '\
        '| samtools sort -@ 24 -o - - '
        '> {output.BAM} '


        
