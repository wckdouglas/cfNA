#!/usr/bin/env python

import sys
import os
import re
import pysam
import subprocess
import pandas as pd
import numpy as np
from functools import lru_cache, partial
from multiprocessing import Pool, Manager
from pybedtools import BedTool, set_tempdir, set_bedtools_path
import random
import RNA
from operator import itemgetter
from sequencing_tools.stats_tools import p_adjust
from sequencing_tools.fastq_tools import reverse_complement
from concensus_seq import concensus
from exon_coverage import ExonFilter
import pyBigWig as pbw
from tblout_parser import read_tbl
from sequencing_tools.bam_tools import get_strand
import pyximport
pyximport.install()
from junction_function import get_junction
import dask.dataframe as dd
import pyranges as pr
import io
from collections import Counter
set_bedtools_path('/stor/work/Lambowitz/cdw2854/src/miniconda3/bin')
set_tempdir('/stor/scratch/Lambowitz/cdw2854')

class GeneMapper():
    def __init__(self):

        self.HB_genes = '''chr11,5246694,5250625,HBB,-
chr11,5253908,5264887,HBD,-
chr11,5269309,5271089,HBG1,-
chr11,5274420,5526835,HBG2,-
chr11,5289575,5526882,HBE1,-
chr16,202686,204502,HBZ,+
chr16,203891,216767,HBM,+
chr14,50037702,50065882,RPS29,-
chr16,222875,223709,HBA2,+
chr16,226679,227521,HBA1,+
chr16,230452,231180,HBQ1,+
chr17,16294748,16295196,FTLP12,+
chr17,21747827,21748211,FTLP13,+
chr19,3870989,3872105,FTLP5,+
chr19,49468558,49470135,FTL,+
chr2,85132749,85133799,TMSB10,+
chr20,44603888,44604714,FTLP,-
chr20,49457152,49457286,TMSL6,-
chr4,69048010,69078188,FTLP10,+
chr7,106809406,106842974,HBP1,+
chr9,131104432,131105049,TMSB4XP4,+
chrX,12993227,12995346,TMSB4X,+
chrX,30648415,30648942,FTLL2,+
chrX,101768604,101771712,TMSB15A,-
chrX,103173479,103221285,TMSB15B,+
chrY,15815447,15817904,TMSB4Y,+'''
        self.hb_gene_df = pd.read_csv(io.StringIO(self.HB_genes),
                            names = ['Chromosome','Start', 'End', 'genes','Strands'])
        self.hb_gene_ranges = pr.PyRanges(self.hb_gene_df)
        self.gene_df = pd.read_csv('/stor/work/Lambowitz/ref/hg19_ref/genes/protein.bed',
                                   usecols = [0,1,2,3,5],
                                   sep='\t',
                                   names = ['Chromosome','Start','End', 'genes','Strands'])
        self.gene_ranges = pr.PyRanges(self.gene_df)

    def test_gene(self, chrom, start, end, return_name=False):
        genes =  self.gene_ranges[chrom, int(start):int(end)]
        answer = 'Not HB'
        HB_name = ''
        if genes:
            HB_name = genes.genes.values[0]
            answer = 'HB'
        
        return answer if not return_name else HB_name

class RBP():
    '''
    Argonaute-associated short introns are a novel
    class of gene regulators
    Hansen, et al. 2016
    '''
    def __init__(self):
        self.RBP = '/stor/work/Lambowitz/ref/hg19_ref/genes/RBP.collapsed.bed.gz'
        self.RBP = pr.PyRanges(pd.read_csv(self.RBP, sep='\t',
                                            names = ['Chromosome',
                                                    'Start',
                                                    'End','RBP', 
                                                    'Strand']))
        

    def search(self, chrom, start, end, strand):
        hits = self.RBP[chrom,strand, int(start):int(end)]
        return 'RBP' if hits else '.'


class Mirtron():
    '''
    http://mirtrondb.cp.utfpr.edu.br/download.php
    '''
    def __init__(self):
        tsv = '/stor/work/Lambowitz/cdw2854/novel_RNA/ref/mirtronDB.tsv'
        self.mirtron = pd.read_csv(tsv, sep='\t', skiprows=2) \
            .rename(columns={'chromosome':'Chromosome',
                            'start':'Start',
                            'end':'End',
                            'strand':'Strand'})\
            .assign(Chromosome = lambda d: 'chr'+d.Chromosome )\
            .filter(['Chromosome','Start','End','Strand','name'])
        self.mirtron = pr.PyRanges(self.mirtron)
    
    def search(self, chrom, start, end, strand):
        hits = self.mirtron[chrom, strand, int(start):int(end)]
        return 'Mirtron' if hits else '.'

class ArgoTron():
    '''
    Argonaute-associated short introns are a novel
    class of gene regulators
    Hansen, et al. 2016
    '''
    def __init__(self):
        link = 'https://media.nature.com/original/nature-assets/ncomms/2016/160513'\
                '/ncomms11538/extref/ncomms11538-s3.xlsx'
        self.argotron = pr.PyRanges(pd.read_excel(link))

    def search(self, chrom, start, end, strand):
        hits = self.argotron[chrom,strand, int(start):int(end)]
        return 'Agotron' if hits else '.'


class TrnaLookAlike():
    '''
    https://cm.jefferson.edu/data-tools-downloads/trna-lookalikes/
    '''

    def __init__(self):
        self.HTML = 'https://cm.jefferson.edu/wordpress/wp-content/uploads/2014/09/tRNA-Lookalikes.txt'
        self.tRNA = pd.read_csv(self.HTML,sep='\t', skiprows=1)\
            .rename(columns={'chromosome':'Chromosome',
                        'FROM_position_inclusive':'Start', 
                        'TO_position_inclusive':'End',
                        'strand':'Strand'})\
            .assign(Chromosome = lambda d:'chr'+ d.Chromosome)
        self.tRNA_lookalike = pr.PyRanges(self.tRNA) 

    def search(self, chrom, start, end, strand):
        hits = self.tRNA_lookalike[chrom,strand, int(start):int(end)]
        return 'tRNA-lookalike' if hits else '.'


class GnomadVar():
    '''
    Genomad vcf snp
    '''
    def __init__(self):
        gnomad = '/stor/work/Lambowitz/ref/hg19_ref/genomad/genomad_grch37.vcf.bgz'
        self.gnomad = pysam.Tabixfile(gnomad)

    def search(self, chrom, start, end):
        try:
            variants = set()
            for variant in self.gnomad.fetch(chrom.strip('chr'), start, end):
                fields = variant.split('\t')
                chrom,pos = fields[1], fields[2]
                variants.add(chrom+':'+pos)
            return len(variants)
        except ValueError:
            return 0


class Mirtrons():
    def __init__(self):
        mirtron_table = '/stor/work/Lambowitz/cdw2854/novel_RNA/ref/mirtron.tsv'
        with open(mirtron_table, 'r') as tab:
            self.mirtrons = set(mir.replace('hsa-mir-','MIR').strip() for mir in tab)


class NameConversion():
    def __init__(self):
        self.encoder = {'RP11-958N24.2': 'PKD1P4-NPIPA8',
                'RP11-1212A22.4':'PKD1P5',
                'RP11-1186N24.5':'PKD1P6',
                'RP11-185O17.3':'RPS2P55',
                'RP11-365K22.1':'RPL13AP25',
                'RP11-111F5.5':'FGF7P6',
                'RP11-1212A22.1':'PKD1P4',
                'RP11-958N24.1':'PKD1P3',
                'RP11-163G10.3':'TMSB4X',
                'AC019188.1':'TMSB4XP8',
                'AP001340.2':'RPL13AP7',
                'RP11-1212A22.4':'PKD1P5'}

    def convert(self, x):
        if x in self.encoder.keys():
            return self.encoder[x]
        else:
            return x

class GenicIntersect():
    def __init__(self):
        self.exons = '/stor/work/Lambowitz/ref/hg19_ref/genes/exons.gencode.bed.gz'
        self.introns = '/stor/work/Lambowitz/ref/hg19_ref/genes/introns.gencode.bed.gz'
        self.pseudogene = '/stor/work/Lambowitz/ref/hg19_ref/genes/small_pseudogenes.bed.gz'
        self.tabix = pysam.Tabixfile(self.introns)
        bam = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/merged_bam/dedup/unfragmented.chrM_filter.dedup.bam'
        self.bam = pysam.Samfile(bam)
       
    def compute_psi(self, chrom, peak_start, peak_end, strand):
        aln_count = 0
        spliced_count = 0
        peak_start, peak_end = int(peak_start), int(peak_end)
        peak_size = peak_end - peak_start
        overlapped_threshold = peak_size * 0.1
        for aln in self.bam.fetch(chrom, peak_start, peak_end):
            if (aln.is_read1 and get_strand(aln) == strand)  or \
                    (aln.is_read2 and get_strand(aln) != strand):
                aln_count += 1
                if 'S' in aln.cigarstring and \
                        self.__test_splice__(aln, peak_start, peak_end, overlapped_threshold):
                    spliced_count+=1
                    
        psi = spliced_count/aln_count if aln_count!=0 else 0
        return psi
    
    def __test_splice__(self, spliced_aln, peak_start, peak_end, th):
        introns = get_junction(spliced_aln)
        coordinates = (intron.split(':')[1].split('_')[0].split('-') for intron in introns)
        overlapped = ( min(peak_end, int(coor[1])) - max(peak_start, int(coor[0])) for coor in coordinates)
        return any(overlap > th  for overlap in overlapped)
        

    
    def intersect(self, bed):
        needed_columns = ['chrom','start','end','peakname', 'pileup','strand']
        self.columns = needed_columns + list(set(bed.columns.tolist()) - set(needed_columns))
        out_columns = self.columns + ['filenum','gchrom','gstart','gend',
                                'gname','gscore','gstrand','gtype','overlap']
        self.inbed = bed\
            .filter(self.columns)
        b_files = [self.exons, self.introns, self.pseudogene] 
        return self.inbed \
            .pipe(BedTool().from_dataframe)\
            .intersect(b = b_files,
                      wao = True)\
            .to_dataframe(names = out_columns)


    def fulllength_intron(self, bed):
        needed_columns = ['chrom','start','end','peakname','score', 'strand']
        cols = needed_columns + list(set(bed.columns.tolist()) - set(needed_columns))
        bed = bed.filter(cols) 
        new_cols = np.append(cols, ['intron_chrom',
                                'intron_start',
                                'intron_end',
                                'intron_name',
                                'intron_score', 
                                'intron_strand','intron_type','overlapped'])
        return BedTool.from_dataframe(bed)\
            .intersect(b = self.introns,
                    f = 0.9, F=0.9, s = True, wao = True)\
            .to_dataframe(names = new_cols) \
            .assign(fulllength_intron = lambda d: np.where(d.intron_name != '.',
                                                'Excised full-length intron',
                                                '.')) \
            .filter(cols + ['fulllength_intron'])\
            .drop_duplicates()\
            .reset_index(drop=True)
            
    def check_intron(self, pchrom, ps, pe, gs, ge):
        rt = 'Exon'
        for intron in self.tabix.fetch(pchrom, ps, pe):
            fields = intron.strip().split('\t')
            istart, iend = int(fields[1]), int(fields[2]) 
            right_continous_exon_intron = (istart - 3) < ge < (istart +3) 
            left_continous_exon_intron = (iend - 3) < gs < (iend + 3)
            right_not_peak_included = iend > pe
            left_not_peak_included = istart < ps
            overlapped = min(pe, iend) - max(ps, istart)
            if (right_continous_exon_intron and right_not_peak_included) \
                    or (left_continous_exon_intron and  left_not_peak_included):
                if overlapped > 10:
                    rt = 'Exon-intron'
        return rt               


    def categorize(self, row):
        gt = row['gtype']
        if row['gtype'] == 'intron':
            if row['overlap_score'] > 0.8:
                gt = 'Full-length intron'
            elif row['psi'] < 0.1:
                gt = 'Within intron'
            elif row['psi'] > 0.2:
                gt =self.check_intron(row['chrom'],
                        row['start'], row['end'],
                        row['gstart'], row['gend'])    
        elif 'pseudo' in row['gtype']:
            gt = 'Pseudogene'
        elif row['gtype'] == 'exon':
            gt =self.check_intron(row['chrom'],
                        row['start'], row['end'],
                        row['gstart'], row['gend'])

        return gt


    def resolve_genic(self, df):
        gbed = df \
            .assign(peak_size = lambda d: d.end - d.start)\
            .assign(gsize = lambda d: d.gend - d.gstart)\
            .set_index('peakname')\
            .pipe(dd.from_pandas, npartitions=24)\
            .groupby(self.columns)\
            .apply(select_annotation, meta = {'gstart':'f8',
                                            'gend':'f8',
                                            'overlap_score':'f8',
                                            'gtype':'f8',
                                            'goverlap':'f8',
                                            'peak_overlap':'f8'})\
            .compute(scheduler='processes')
        gbed.index = gbed.index.droplevel(-1)
        gbed = gbed.reset_index()\
            .assign(gtype = lambda d: [self.categorize(row) for i, row in d.iterrows()])   \
            .assign(gtype = lambda d: np.where((d.gtype=='Exon-intron') & (d.peak_overlap.astype(float) > 0.6),
                                            'Exon',
                                            d.gtype))
        return gbed


    def Annotate(self, df):
        trna=TrnaLookAlike()
        argotron = ArgoTron()
        mirtron = Mirtron()
        rbp = RBP()

        return df\
            .assign(psi = lambda d: list(map(self.compute_psi, d.chrom ,d.start, d.end, d.strand)) )\
            .pipe(self.intersect)\
            .pipe(self.resolve_genic) \
            .assign(is_tRNA_lookalike = lambda d: list(map(trna.search, d.chrom, d.start, d.end, d.strand))) \
            .assign(is_agotron = lambda d: list(map(argotron.search, d.chrom, d.start, d.end, d.strand))) \
            .assign(is_mirtron = lambda d: list(map(mirtron.search, d.chrom, d.start, d.end, d.strand))) \
            .assign(is_rbp = lambda d: list(map(rbp.search, d.chrom, d.start, d.end, d.strand))) \
            .pipe(self.fulllength_intron) 


 
def select_annotation(df):
    return df \
        .assign(peak_overlap = lambda d: d.overlap / (d.peak_size))\
        .assign(goverlap = lambda d: d.overlap / (d.gsize))\
        .assign(overlap_score = lambda d:  d.peak_overlap * d.goverlap )\
        .fillna(0.000001) \
        .nlargest(1,'overlap_score')  \
        .filter(['gstart','gend','overlap_score','gtype','goverlap','peak_overlap']) 



class mRNAFilter():
    '''
    if the peak is on exon?
    is the peak also called in transcriptome?
    '''
    def __init__(self):
        ref_path = '/stor/work/Lambowitz/ref/hg19_ref/genes'
        exons = ref_path + '/gencode.exon.bed.gz'
        self.exons = pysam.Tabixfile(exons)
        transcriptom_peaks = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/transcriptome/macs2/unfragmented.fwd_peaks_genomics.narrowPeak.gz'
        self.transcriptome_peaks = pysam.Tabixfile(transcriptom_peaks)

    def search(self, chrom, start, end, attribute = 'exon'):
        if attribute == 'exon':
            it = self.exons
        elif attribute == 'transcriptome':
            it = self.transcriptome_peaks
        return 'yes' if any(True for frag in it.fetch(chrom, start, end) if  self.__overlap_test__(start, end, frag)) else 'no'

    def __overlap_test__(self, start, end, fragment):
        fields = fragment.split('\t')
        return int(fields[1]) < start and int(fields[2]) > end


class WPS:
    '''
    Window protection scores
    '''
    def __init__(self, start, end, window = 20):
        self.start = int(start)
        self.end = int(end)
        self.window = window
        self.half_window = int(self.window/2)
        self.range = self.end - self.start + self.window
        self.WPS_array = np.zeros(self.range)

    def add_fragment(self, frag_start, frag_end):
        frag_start, frag_end = int(frag_start), int(frag_end)
        frag_size = frag_end - frag_start
        self.fragment_holder = np.zeros(frag_size + self.window)
        self.fragment_holder[self.half_window:-self.half_window] = 1
        self.fragment_holder[self.fragment_holder != 1] = -1
        offset = self.start - (frag_start - self.half_window)  # how far is the aln start relative to peak start
        if offset > 0:
            wps = self.fragment_holder[offset:]
            end = len(wps) if len(wps) < self.range else self.range
            self.WPS_array[:end] += wps[:end]
        else:
            baseShifted = abs(offset)
            end = self.range if baseShifted + len(self.fragment_holder) > self.range else baseShifted + len(self.fragment_holder)
            alignedBases = self.range + offset
            wps = self.fragment_holder[:alignedBases]
            self.WPS_array[baseShifted:end] += wps 
        



class PeakAnalyzer:
    def __init__(self,
            genome = '/stor/work/Lambowitz/ref/hg19_ref/genome/hg19_genome.fa',
            RNAshapes = '/stor/work/Lambowitz/cdw2854/src/miniconda3/bin/RNAshapes',
            gene_bed = '/stor/work/Lambowitz/ref/hg19_ref/genes/genes.bed.gz',
            sample_bed = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/bed_files/merged_bed/unfragmented.bed.gz',
            tblout = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/bed_files/merged_bed/MACS2/annotated/unfragmented.tblout',
            phyloP = '/stor/work/Lambowitz/ref/hg19_ref/phyloP/hg19.100way.phastCons.bw'):
        self.fa = pysam.Fastafile(genome)
        self.gene_bed = gene_bed
        self.RNAshapes = RNAshapes
        self.sample_bed = pysam.Tabixfile(sample_bed)
        self.phyloP = pbw.open(phyloP)
        self.tblout = tblout
        self.mRNA_filter = mRNAFilter()
        self.gene_bed_columns = ['chrom','start','end', 'gene_name',
                                'gene_score','strand','gene_type','gene_id']
        self.gene_bed_tab = pysam.Tabixfile(gene_bed)

    def fetch_seq(self, chrom, start, end, strand):
        seq = self.fa.fetch(chrom, int(start), int(end))
        return seq if strand == "+" else reverse_complement(seq)

    
    def add_rfam(self, peak_df):
        tbl = read_tbl(self.tblout)\
            .pipe(lambda d: d[d['E-value'] < 0.01])\
            .groupby('query name', as_index=False)\
            .apply(lambda d: d.nsmallest(1,'E-value'))\
            .filter(regex='name')\
            .assign(peak_name = lambda d: d['query name'].str.extract('_(chr[0-9XY]+:[0-9]+\-[0-9]+)\(')) \
            .drop(['query name'],axis=1)\
            .rename(columns = {'target name':'rfam'})
        return peak_df\
                .merge(tbl, on = 'peak_name', how = 'left') \
                .fillna('.')

    
    def abstract_structure(self, structure):
        if set(structure) == {'(', ')', '.'}:
            s = subprocess.check_output([self.RNAshapes,'-D',structure,'-t','5']).decode()
        else:
            s = ''
        return s.strip()

    def RNAfold(self, seq, return_energy=False):
        fold = ''
        energy = 0
        if seq != '':
            fold, energy = RNA.fold(seq)
        return energy if return_energy else fold

    def full_length(self, chrom, start, end, strand):
        frag_count = 0
        fulllength = 0
        start, end = int(start), int(end)
        for frag in self.sample_bed.fetch(chrom, start, end):
            fields = frag.split('\t')
            frag_strand = fields[5]
            if frag_strand == strand:
                frag_count += 1
                if start -5 < int(fields[1]) < start + 5 and end -5 < int(fields[2]) < end + 5:
                    fulllength += 1

        if frag_count == 0:
            return 0
        return fulllength


    def mirtron_filter(self, peak_tab, remove_intron=True):
        final_columns = peak_tab.columns 
        needed_columns = ['chrom','start', 'end','peakname','pileup','strand']
        extra_columns = set(final_columns) - set(needed_columns)
        needed_columns.extend(list(extra_columns))
        needed_columns.extend(['intron_chrom',
                                'intron_start',
                                'intron_end',
                                'intron_name',
                                'intron_score',
                                'intron_strand',
                                'intron',
                                'overlapped'])
        
        intersected = BedTool()\
            .from_dataframe(peak_tab.filter(needed_columns))\
            .intersect('/stor/work/Lambowitz/ref/hg19_ref/genes/introns.gencode.bed.gz', 
                    f= 0.5,F=0.9, wao = True, s=True)\
            .to_dataframe(names = needed_columns)  \
            .pipe(lambda d: pd.concat([d.query('intron_chrom == "."'),
                                        d\
                                            .query('intron_chrom != "."')\
                                            .assign(fulllength = lambda d: list(map(self.full_length, d.intron_chrom, d.intron_start, d.intron_end, d.strand)))\
                                            .assign(intron_chrom = lambda d: np.where(d.fulllength < 1, '.',d.intron_chrom))\
                                            .drop('fulllength', axis=1)]))  

        if remove_intron:
            return intersected \
                .query('intron_chrom=="."')\
                .pipe(lambda d: d.drop(list(d.filter(regex='^intron|overlap')), axis=1))\
                .drop_duplicates()
        else:
            return intersected\
                .assign(is_intron = lambda d: np.where(d.intron_chrom!=".", 'full-length intron','.'))\
                .pipe(lambda d: d.drop(list(d.filter(regex='^intron|overlap')), axis=1))\
                .drop_duplicates()

    def filter_mRNA(self, peak_tab):
        return peak_tab \
            .assign(is_exon = lambda d: [self.mRNA_filter.search(chrom, start, end, attribute = 'exon') for chrom, start, end in zip(d.chrom, d.start, d.end)]) \
            .assign(is_transcriptome_peak = lambda d: [self.mRNA_filter.search(chrom, start, end, attribute = 'transcriptome') for chrom, start, end in zip(d.chrom, d.start, d.end)])\
            .pipe(lambda d: d[(d.is_exon=="no") | ((d.is_exon=="yes") & (d.is_transcriptome_peak=="yes"))])

    def filter_gene(self, peak_tab, gene_regex='^HB[ABDGEMQP]$|^HB[ABDGEMQP][0-9]+$|TMSB4X|FTLP3$|RPS29$|RPL6P27'):
        final_columns = peak_tab.columns 
        needed_columns = ['chrom','start', 'end','peakname','pileup','strand']
        extra_columns = set(final_columns) - set(needed_columns)
        needed_columns.extend(list(extra_columns))
        peak_bed = peak_tab.filter(needed_columns)
        peaks = BedTool().from_dataframe(peak_bed)

        bed = pd.read_csv(self.gene_bed, sep='\t',
                            names = self.gene_bed_columns)\
            .pipe(lambda d: d[d.gene_name.str.contains(gene_regex)])
        gene = BedTool().from_dataframe(bed)
        return peaks.intersect(b = gene, v=True) \
            .to_dataframe(names = needed_columns)
    

    def add_genename(self, chrom, start, end, strand):
        gene_name = ''
        max_overlap = 0
        start, end = int(start), int(end)
        for gene in self.gene_bed_tab.fetch(chrom, start, end):
            fields = gene.strip().split('\t')
            name = fields[3]
            gstart, gend = int(fields[1]), int(fields[2])
            overlap = min(gend,end) - max(start, gstart) 
            if overlap > max_overlap and strand == fields[5]:
                max_overlap = overlap
                gene_name = name
        return gene_name
        
    def cloverleaf(self, seq, return_df=False):
        if not seq:
            return 'N/A'

        cmd = 'echo '+seq+ \
            '| '+ self.RNAshapes + ' -e 20  -m "[[][][]]"' 
        ps = subprocess.Popen(cmd, shell=True, 
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT)
        res = ps.communicate()[0].decode() 
        if 'not found' in res:
            return 'N/A'
        else:
            res_list = re.sub('\s+','\t',res.strip().split('\n')[1]).split('\t')
            return res_list[1]
        
    def uncharacterized_peak(self, peak_df, remove_intron=True):
        rows = []
        for i, row in peak_df.iterrows():
            is_intron = re.search('[iI]ntron', row['gtype'])
            is_sense = row['sense_gtype'] != '.'
            not_rbp = row['sense_gtype'] != 'RBP'
            is_antisense = row['sense_gtype'] == 'Antisense'
            is_unannotated = row['sense_gtype'] == '.' and row['antisense_gtype'] == "."
            if (is_sense and is_intron and not_rbp) or is_antisense or is_unannotated :
                rows.append(row)
        return pd.DataFrame(rows)

    
    def __add_wps__(self, chrom, start, end, strand):
        window = 10
        wps = WPS(start, end, window = window)
        for fragment in self.sample_bed.fetch(chrom, int(start), int(end)):
            fields = fragment.strip().split('\t')
            frag_chrom, frag_start, frag_end, frag_strand = itemgetter(0,1,2,5)(fields)
            if strand == frag_strand:
                wps.add_fragment(frag_start, frag_end)
        max_WPS = wps.WPS_array[wps.WPS_array > 0.8* wps.WPS_array.max()]
        return len(max_WPS)/(len(wps.WPS_array) - window)


    def wps_function(self, peak_df):
        return peak_df \
            .assign(max_wps = lambda d: list(map(self.__add_wps__, d.chrom, d.start, d.end, d.strand)))


    def add_phyloP(self, chrom, start, end):
        return self.phyloP.values(chrom, int(start), int(end), numpy=True).mean()



    def __max_percentage__(self, vector):
        vector = np.array(list(vector))
        return vector.max()/vector.sum()


    def __flush_ends__(self, chrom, start, end, strand):
        start_positions = Counter()
        end_positions = Counter()
        for frag in self.sample_bed.fetch(chrom, start, end):
            frag_start, frag_end, frag_strand = itemgetter(1,2,5)(frag.rstrip().split('\t'))
            if frag_strand == strand:
                start_positions[frag_start] += 1
                end_positions[frag_end] += 1
        start_percentage = self.__max_percentage__(start_positions.values())
        end_percentage = self.__max_percentage__(end_positions.values())
        return start_percentage, end_percentage

    def add_flush_end(self, df):
        res = map(self.__flush_ends__, df.chrom, df.start, df.end, df.strand)
        max_start, max_end = zip(*res)
        df.loc[:, 'max_start_fraction'] = max_start
        df.loc[:, 'max_end_fraction'] = max_end
        return df
        


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
    if seq != '':
        folded, energy = RNA.fold(seq)
        sequence_bases = list(seq)
                
        random_energies = [__fold_random__(sequence_bases) for sim in range(simulation)]
        random_energies = np.array(random_energies)
        return len(random_energies[random_energies <= energy])/simulation
    else:
        return 1

def mirtron(peak_path):
    peak_file = peak_path + '/mirtron_peak.tsv'
    print('Read %s' %peak_file)
    return pd.read_csv(peak_file) \
        .assign(peakname = lambda d: d.peakname + '_mirtron')
    


def concensus_module():
    workdir = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
    fwd = workdir + '/bed_files/merged_bed/coverage/unfragmented.fwd.bigWig'
    rvs = workdir + '/bed_files/merged_bed/coverage/unfragmented.rvs.bigWig'
    bam = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/merged_bam/dedup/unfragmented.chrM_filter.dedup.bam'
    return concensus(bam, coverage_files = [fwd, rvs])

def main(remove_intron=False, exon=False):
    #remove_intron=False
    #exon=True
    WORK_DIR = os.environ['WORK'] + '/cdw2854/cfNA/tgirt_map/bed_files/merged_bed/MACS2/annotated'
    ANNOT_PEAK_FILE = WORK_DIR + '/unfragmented.filtered.tsv'
    peak_analyze = PeakAnalyzer()
    concensus_analyzer = concensus_module()
    exon_filter = ExonFilter()
    gi = GenicIntersect()
    name_conversion = NameConversion()
    gene_mapper = GeneMapper()
    gnomad = GnomadVar()
    trna=TrnaLookAlike()
    argotron = ArgoTron()
    rbp = RBP()

    out_columns = ['peak_name','is_sense','gname','sense_gtype', 'strand','pileup','sample_count', 'seq', 'energy',
                 'concensus_seq', 'concensus_fold', 'concensus_energy', #'concensus_folding_pval',
                 'force_cloverleaf','is_exon','is_transcriptome_peak','rfam',  'is_mirtron','protein_gene',
                 'log10p', 'max_start_fraction', 'max_end_fraction','num_var', 'gtype', 'is_tRNAlookalike', 'is_argotron','is_rbp'] 
    if not remove_intron:
        out_columns.append('is_intron')

    fold_p = partial(get_folding_pvalue, 2001)
    peak_df = pd.read_csv(ANNOT_PEAK_FILE, sep='\t') \
        .pipe(lambda d: pd.concat([d, mirtron(WORK_DIR)])) \
        .query('sample_count >= 5 & pileup >= 5') \
        .fillna('.')\
        .assign(peak_name = lambda d: d.chrom +':'+ d.start.astype(str) +'-'+ d.end.astype(str)) \
        .pipe(peak_analyze.add_rfam)\
        .pipe(peak_analyze.mirtron_filter, remove_intron=remove_intron) \
        .assign(sense_gtype = lambda d: np.where((d.sense_gtype=="RBP") & (d.is_intron =="full-length intron"),
                                                    'full-length intron (RBP)',
                                                    d.sense_gtype))\
        .query('sense_gtype != "RBP" & sense_gtype != "Repeats"')\
        .assign(psi = lambda d: list(map(gi.compute_psi, d.chrom ,d.start, d.end, d.strand)) )\
        .pipe(gi.intersect) \
        .pipe(gi.resolve_genic) \
        .assign(is_sense = lambda d: list(map(label_sense, d.sense_gname, d.antisense_gname))) 

    if exon:
        peak_df = peak_df\
            .query('is_sense=="Sense"') \
            .pipe(lambda d: d[d.gtype.isin(['Exon',"Pseudogene"])]) 

    else:
        peak_df = peak_df\
            .pipe(lambda d: d[~((d.gtype == "Pseudogene") & (d.is_sense=="sense"))]) \
            .pipe(lambda d: d[~((d.gtype == "Exon") & (d.is_sense=="Sense"))])

    peak_df = peak_df \
        .assign(concensus_seq = lambda d: list(map(concensus_analyzer.find_concensus, d.chrom, d.start, d.end, d.strand)))\
        .assign(seq = lambda d: list(map(peak_analyze.fetch_seq, d.chrom, d.start, d.end, d.strand)))\
        .assign(concensus_seq = lambda d: np.where(d.is_intron==".",d.concensus_seq,d.seq)) \
        .assign(fold = lambda d: d.seq.map(peak_analyze.RNAfold))\
        .assign(energy = lambda d: d.seq.map(lambda x: peak_analyze.RNAfold(x, return_energy=True)))\
        .assign(concensus_fold = lambda d: d.concensus_seq.map(peak_analyze.RNAfold))\
        .assign(concensus_energy = lambda d: d.concensus_seq.map(lambda x: peak_analyze.RNAfold(x, return_energy=True)))\
        .assign(fold = lambda d: d.fold.map(peak_analyze.abstract_structure)) \
        .assign(concensus_fold = lambda d: d.concensus_fold.map(peak_analyze.abstract_structure)) \
        .assign(force_cloverleaf = lambda d: d.concensus_seq.map(peak_analyze.cloverleaf)) \
        .assign(sense_gname = lambda d: d.sense_gname.map(name_conversion.convert)) \
        .assign(gname = lambda d: np.where(d.sense_gname != ".", 
                                        d.sense_gname,
                                        d.antisense_gname)) \
        .assign(is_mirtron = lambda d: np.where(d.peakname.str.contains('mirtron'), 'mirtron','.'))\
        .assign(protein_gene = lambda d: [gene_mapper.test_gene(chrom, start, end, return_name=True) \
                                    for chrom, start, end in zip(d.chrom, d.start, d.end)]) \
        .pipe(peak_analyze.add_flush_end)\
        .assign(is_argotron = lambda d: list(map(argotron.search, d.chrom, d.start, d.end, d.strand)))\
        .assign(is_tRNAlookalike = lambda d: list(map(trna.search, d.chrom, d.start, d.end, d.strand))) \
        .assign(is_rbp = lambda d: list(map(rbp.search, d.chrom, d.start, d.end, d.strand))) \
        .assign(num_var = lambda d: list(map(gnomad.search, d.chrom, d.start, d.end)))\
        .filter(out_columns)\
        .rename(columns = {'pileup':'fragment_count'}) \
        .pipe(lambda d: d[~d.peak_name.isin(['chrX:44654061-44654153', 
                                            'chr11:74457163-74457239',
                                            'chr16:31173304-31173386'])])  
        #.assign(concensus_folding_pval = lambda d: p.map(fold_p, d.concensus_seq)) \
        #.pipe(peak_analyze.filter_mRNA)\
        #.pipe(peak_analyze.filter_gene)\
    out_table = WORK_DIR + '/supp_tab.tsv'
    if not remove_intron:
        out_table = WORK_DIR + '/supp_tab_intron.tsv'
        if exon:
            out_table = out_table.replace('_intron.tsv','_exon.tsv')
    else:
        if exon:
            out_table = out_table.replace('.tsv','_exon.tsv')
    peak_df.to_csv(out_table, index=False, sep='\t')
    print('Written: ', out_table)
        


if __name__ == '__main__':
    main(remove_intron=False, exon=True)
    main(remove_intron=False, exon=False)
