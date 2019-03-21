#!/usr/bin/env python

import os
import pysam
import subprocess
import pandas as pd
import numpy as np
from functools import lru_cache, partial
from multiprocessing import Pool, Manager
from pybedtools import BedTool, set_tempdir
import random
import RNA
from sequencing_tools.stats_tools import p_adjust
from sequencing_tools.fastq_tools import reverse_complement
from concensus_seq import concensus
set_tempdir(os.environ['SCRATCH'])

class peak_analyzer:
    def __init__(self):
        self.fa = pysam.Fastafile('/stor/work/Lambowitz/ref/hg19_ref/genome/hg19_genome.fa')
        self.gene_bed = '/stor/work/Lambowitz/ref/hg19_ref/genes/genes.bed'

    def fetch_seq(self, chrom, start, end, strand):
        seq = self.fa.fetch(chrom, int(start), int(end))
        return seq if strand == "+" else reverse_complement(seq)

    
    def abstract_structure(self, structure):
        s = subprocess.check_output(['RNAshapes','-D',structure,'-t','5'])
        return s.decode().strip()


    def filter_gene(self, peak_tab, gene_regex='^HB[ABGEMQP]$|^HB[ABGEMQP][0-9]+$|TMSB4'):
        final_columns = peak_tab.columns 
        needed_columns = ['chrom','start', 'end','peakname','pileup','strand']
        extra_columns = set(final_columns) - set(needed_columns)
        needed_columns.extend(list(extra_columns))
        peak_bed = peak_tab.filter(needed_columns)
        peaks = BedTool().from_dataframe(peak_bed)

        bed = pd.read_table(self.gene_bed, 
                            names = ['chrom','start','end', 'gene_name',
                                'gene_score','strand','gene_type','gene_id'])\
            .pipe(lambda d: d[d.gene_name.str.contains(gene_regex)])
        gene = BedTool().from_dataframe(bed)
        return peaks.intersect(b = gene, v=True) \
            .to_dataframe(names = needed_columns)

def label_sense(picked_type_sense, picked_type_anti):
    label = '.'
    if picked_type_sense not in ['.','Unannotated']:
        label = 'Sense'
    elif picked_type_anti not in ['.','Unannotated'] :
        label = "Antisense"
    else:
        label = 'Unannotated'
    return label


def __fold_random__(sequence_bases):
    random.shuffle(sequence_bases)
    f, e = RNA.fold(''.join(sequence_bases))
    return e


@lru_cache(maxsize=128)
def get_folding_pvalue(simulation, seq):
    folded, energy = RNA.fold(seq)
    sequence_bases = list(seq)
            
    random_energies = [__fold_random__(sequence_bases) for sim in range(simulation)]
    random_energies = np.array(random_energies)
    return len(random_energies[random_energies < energy])/simulation


def concensus_module():
    workdir = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
    fwd = workdir + '/bed_files/merged_bed/coverage/unfragmented.fwd.bigWig'
    rvs = workdir + '/bed_files/merged_bed/coverage/unfragmented.rvs.bigWig'
    bam = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/merged_bam/dedup/unfragmented.chrM_filter.dedup.bam'
    return concensus(bam, coverage_files = [fwd, rvs])

def main():
    WORK_DIR = os.environ['WORK'] + '/cdw2854/cfNA/tgirt_map/bed_files/merged_bed/MACS2/annotated'
    ANNOT_PEAK_FILE = WORK_DIR + '/unfragmented.filtered.tsv'
    peak_analyze = peak_analyzer()
    concensus_analyzer = concensus_module()


    p = Pool(24)
    fold_p = partial(get_folding_pvalue, 2000)
    peak_df = pd.read_table(ANNOT_PEAK_FILE) \
        .query('sample_count >= 5 & pileup >= 4') \
        .fillna('.')\
        .query('sense_gtype == "Long RNA" | sense_gtype =="."')\
        .query('antisense_gtype != "Repeats" | antisense_gtype != "RBP"' )\
        .pipe(peak_analyze.filter_gene) \
        .assign(concensus_seq = lambda d: list(map(concensus_analyzer.find_concensus, d.chrom, d.start, d.end, d.strand)))\
        .assign(seq = lambda d: list(map(peak_analyze.fetch_seq, d.chrom, d.start, d.end, d.strand)))\
        .assign(fold = lambda d: list(map(lambda x: RNA.fold(x)[0], d.seq)))\
        .assign(concensus_fold = lambda d: list(map(lambda x: RNA.fold(x)[0], d.concensus_seq)))\
        .assign(fold = lambda d: d.fold.map(peak_analyze.abstract_structure)) \
        .assign(concensus_fold = lambda d: d.concensus_fold.map(peak_analyze.abstract_structure)) \
        .assign(folding_pval = lambda d: p.map(fold_p, d.concensus_seq)) \
        .assign(peak_name = lambda d: d.chrom +':'+ d.start.astype(str) +'-'+ d.end.astype(str)) \
        .assign(gname = lambda d: np.where(d.sense_gname != ".", 
                                        d.sense_gname,
                                        d.antisense_gname)) \
        .assign(is_sense = lambda d: list(map(label_sense, d.sense_gname, d.antisense_gname))) \
        .filter(['peak_name','is_sense','gname', 'strand','pileup','sample_count', 'seq', 'fold', 'folding_pval'])
    p.close()
    p.join()
    peak_df.to_csv(WORK_DIR + '/supp_tab.tsv', index=False, sep='\t')
        


if __name__ == '__main__':
    main()
