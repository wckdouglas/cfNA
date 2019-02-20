import pandas as pd
import numpy as np
from sequencing_tools.viz_tools import okabeito_palette, color_encoder, simpsons_palette
from sequencing_tools.stats_tools import p_adjust
from scipy.special import ndtr
from collections import defaultdict
from sequencing_tools.fastq_tools import reverse_complement
from sequencing_tools.bam_tools import get_strand
import RNA
from multiprocessing import Pool
import random
import pysam
import glob
import re
from pybedtools import BedTool
from plotting_utils import figure_path
import seaborn as sns
import mappy
from tblout_parser import read_tbl
from bwapy import BwaAligner
import io


pileup_cutoff = 4
sample_cutoff = 5
project_path = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
peak_path = project_path + '/bed_files/merged_bed/MACS2/annotated'

def only_choice(row):
    # cases where only 1 potential RNA is found
    if row['gstrand'] == row['strand']:
        row['picked_RNA_sense'] = row['gname']
        row['picked_type_sense'] = row['gtype']
        row['picked_RNA_anti'] = '.'
        row['picked_type_anti'] = '.'
    else:
        row['picked_RNA_sense'] = '.'
        row['picked_type_sense'] = '.'
        row['picked_RNA_anti'] = row['gname']
        row['picked_type_anti'] = row['gtype']
    return row
    
def peak_info_table(row):
    peak_info = pd.DataFrame({'RNA': row['gname'].split(','),
                           'strand': row['gstrand'].split(','),
                           'peak_type': row['gtype'].split(',')}) \
        .assign(peak_rank = lambda d: d.peak_type.map(rank_type))\
        .assign(strand = lambda d: d.strand.astype(str))\
        .sort_values('peak_rank') \
        .reset_index()
    return peak_info
    

full_length_mRNA = '^HIST|^FT[LH]|^RP[LS]'
full_length_regex = re.compile(full_length_mRNA)
def rank_peaks(row):
    # multiple peak annotations
    peak_info = peak_info_table(row)
    
    for strand in ['-','+']:
        strand_peak_info = peak_info[(peak_info['strand'] == strand) | (peak_info['peak_type'].str.contains('Low|Simple|DNA'))]
        #print(strand, peak_info)
        stranding = '_sense' if row['strand'] == strand else '_anti'
        
        picked_RNA = 'picked_RNA' + stranding
        picked_type = 'picked_type' + stranding
        
        if strand_peak_info.shape[0] > 0:    
            row[picked_RNA] = strand_peak_info['RNA'].values[0]
            row[picked_type] = strand_peak_info['peak_type'].values[0]
        else:
            row[picked_RNA] = '.'
            row[picked_type] = '.'
            
            
        # correct for full length RNA
        full_mRNA_df = strand_peak_info[strand_peak_info.RNA.str.contains(full_length_mRNA)]
        if full_mRNA_df.shape[0] > 0 and row[picked_type] == "piRNA":
            row[picked_RNA] = full_mRNA_df.RNA.values[0]
            row[picked_type] = 'protein_coding'
    return row
            
        
def peak_assignment(args):
    # picked best represented peak type
    i, row = args
    
    RNA = ''
    peak_type = ''
    
    if ',' not in row['gstrand']:
        row = only_choice(row)
    else:
        # cases where several potential RNA
        row = rank_peaks(row)
    return row
    

lrna_regex = 'lincR|protein|pseudo|TR|proces|sense_intr'\
            'prime|IG|antisen|lncRNA|sense_ov|TEC'   
def merge_type(x):
    if re.search('LINE|Satellite|Simple_repeat|SINE|Unknown'
                 '|Low_complexity|LTR|^DNA$|^DNA\?$|RC|Other', x):
        return 'Repeats'

    elif x == ".":
        return 'Unannotated'
    
    elif re.search(lrna_regex, x):
        return 'Long RNA'
        
    elif re.search('rRNA|rDNA', x):
        return 'rRNA'
    elif re.search('misc|guid|scRN|srpRNA', x):
        return 'misc RNA'
    else:
        return x


def label_sense(picked_type_sense, picked_type_anti):
    label = '.'
    if picked_type_sense not in ['.','Unannotated']:
        label = 'Sense'
    elif picked_type_anti not in ['.','Unannotated'] :
        label = "Antisense"
    else:
        label = 'Unannotated'
    return label

def load_peaks_old(filename):
    peak_df = pd.read_table(filename) 
    p = Pool(24)
    peaks = pd.DataFrame(p.map(peak_assignment, peak_df.iterrows()))
    p.close()
    p.join()
    peaks = peaks \
        .assign(pvalue = lambda d: np.power(10, -d.log10p))\
        .assign(picked_type_sense = lambda d: np.where((d.picked_RNA_sense.str.contains('HY')) & \
                                                        (d.picked_type_sense == 'scRNA'), 
                                             'misc RNA',
                                             d.picked_type_sense))\
        .assign(picked_type_sense = lambda d: np.where(d.picked_RNA_sense.str.contains('srpRNA'), 
                                             'misc RNA',
                                             d.picked_type_sense))\
        .assign(picked_type_sense = lambda d: np.where(d.picked_RNA_sense == "BC200", 
                                                       'misc RNA', 
                                                       d.picked_type_sense))\
        .assign(merged_type = lambda d: d.picked_type_sense.map(merge_type)) \
        .assign(FDR = lambda d: p_adjust(d.pvalue) ) \
        .assign(is_sense = lambda d: list(map(label_sense, d.picked_RNA_sense, d.picked_RNA_anti))) \
        .fillna('.')
    return peaks

def load_peaks(filename):
    EPSILON = 1e-100
    peak_df = pd.read_table(filename) 
    if 'log10p' in peak_df.columns:
        peak_df['pvalue'] = np.power(10, -peak_df.log10p)

    elif 'zscore' in peak_df.columns:
        pvalue = ndtr(-peak_df.zscore) 
        pvalue[pvalue==0] = EPSILON
        peak_df['pvalue'] = pvalue
        peak_df['log10p'] = peak_df.pvalue.transform(lambda x: -np.log10(x))

    peak_df = peak_df \
        .assign(FDR = lambda d: p_adjust(d.pvalue) ) \
        .fillna('.')\
        .assign(is_sense = lambda d: list(map(label_sense, d.sense_gtype, d.antisense_gtype))) 
    return peak_df


def plot_peak_strand(peaks, ax):
    peaks\
        .query('pileup>=%i & sample_count >= %i' %(pileup_cutoff, sample_cutoff))\
        .groupby('is_sense')\
        .agg({'sense_gtype':'count'})\
        .reset_index() \
        .assign(index = lambda d: d['is_sense'] + '\n(' + d.sense_gtype.astype(str)+')')\
        .set_index('index')\
        .sort_values('sense_gtype', ascending=True)\
        .plot.pie(y='sense_gtype', 
                explode = [0,0.05, 0.2], 
                startangle = 180, ax = ax,
                colors = ['#2a7a1c','#40649e','#f7b707'])
    ax.set_ylabel('')
    ax.legend().set_visible(False)


def change_annotation(lab):
    if 'RBP' in lab:
        return lab.replace('RBP','Long RNA (RBP)')
    elif 'Long RNA' in lab:
        return lab.replace('Long RNA', 'Long RNA (narrow peak)')
    else:
        return lab



def plot_peak_pie(peaks, ax, ce, gtype='sense_gtype'):
    peak_pie= peaks\
        .query('pileup>=%i & sample_count >= %i' %(pileup_cutoff, sample_cutoff))\
        .groupby(gtype, as_index=False)\
        .agg({'pvalue':'count'}) \
        .assign(fraction = lambda d: d.pvalue.transform(lambda x: 100*x/x.sum())) \
        .assign(rtype = lambda d: d[gtype]+ ' (' + d.pvalue.astype(str) + ')')\
        .assign(merged_type = lambda d: d.rtype.map(change_annotation))\
        .set_index('merged_type')\
        .assign(explode = lambda d: (100-d.fraction)/100) \
        .assign(explode = lambda d: np.where(d.explode < 0.95,0, 
                                         np.where(d.explode < 0.99, 0.2,
                                                  np.where(d.explode < 0.994, 0.4, 0.8))))\
        .sort_values('pvalue', ascending=False)
    
    rna_types = list(map(lambda x: x.split('(')[0].strip(), peak_pie.rtype))
    colors = pd.Series(rna_types).map(ce.encoder)
    peak_pie.plot(kind = 'pie',
              y = 'fraction', 
              ax = ax,
              explode = peak_pie.explode,
              colors = colors)
    
    index = peak_pie.index
    ax.legend(bbox_to_anchor = (1.5,0.7), 
              labels = list(map(lambda x: x.split(' ')[0], index)), 
              ncol=2, 
              fontsize=12)
    
    ax.set_ylabel('')
    ax.legend().set_visible(False)

def plot_peak_bar(ax,peaks):
    for i, row in combined_peaks \
            .query('is_sense == "Sense" & pileup >= %i & sample_count >= %i' %(pileup_cutoff, sample_cutoff)) \
            .groupby(['annotation','merged_type'],as_index=False) \
            .agg({'pvalue':'count'}) \
            .pipe(pd.pivot_table, index = 'merged_type', 
                    columns ='annotation', values = 'pvalue')\
            .iterrows():
                
        alpha = 1 if i in ['Long RNA','RBP'] else 0.5
        ax.plot([1,2], [row['K562'], row['K562 + HepG2']],
               color = ce.encoder[i], alpha=alpha)
        ax.scatter([1,2], [row['K562'], row['K562 + HepG2']],
               color = ce.encoder[i], alpha=alpha)
        x = 3 if i in ['misc RNA','miRNA'] else 2.1 
        ax.text(x, row['K562 + HepG2'], s = i, color = ce.encoder[i], fontsize=15)
    ax.set_xticks([1,2])
    ax.set_xticklabels(['k562 only', 'K562 + HepG2'], 
                   rotation = 70, rotation_mode = 'anchor', ha = 'right')
    ax.set_xlim(0,3)
    ax.set_xlabel('ENCODE RNA binding-protein\nbinding site annotation')
    ax.set_ylabel('Number of peaks')
    sns.despine() 


def assign_rna_name(x):
    if re.search('RNY|Y_RNA|^HY', x):
        return 'Y-RNA'
    elif re.search('7SL|SRP', x):
        return '7SL RNA'
    elif re.search('7SK', x):
        return '7SK RNA'
    elif re.search('VTR', x):
        return 'Valut RNA'
    else:
        return x

def plot_repeats_RNA(peaks, ax, ce, rnatype="Repeats", top_n = 10):
    peaks\
        .query('sense_gtype == "%s"' %rnatype)\
        .query('pileup >= %i & sample_count >= %i' %(pileup_cutoff, sample_cutoff)) \
        .assign(RNA_name = lambda d: d.sense_gname.map(assign_rna_name)) \
        .groupby('RNA_name', as_index=False)\
        .agg({'chrom':'count'}) \
        .nlargest(top_n, 'chrom')\
        .plot.bar('RNA_name','chrom', color = ce.encoder[rnatype], ax = ax)
    ax.set_xlabel('')
    ax.set_ylabel('Peak count')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=70, rotation_mode='anchor', ha = 'right')
    ax.legend().set_visible(False)


def color_rbp(ax):
    '''
    example: ','.join(rbp_df.head(15).index)
    http://plasmaproteomedatabase.org/
    '''

    for xt in ax.get_xticklabels():
        if xt.get_text() not in ['IGF2BP1','LARP4','LIN28B']:
            xt.set_color('red')


def plot_rbp(peaks, ax, ce, top_n = 10):
    rbp_count = defaultdict(int)
    for i, row in peaks.query('pileup >= %i & sample_count >= %i' %(pileup_cutoff, sample_cutoff)).iterrows():
        added = set()
        for btype, bp in zip(row['sense_gtype'].split(','),
                        row['sense_gname'].split(',')):
            if btype == "RBP":
                if bp not in added:
                    rbp_count[bp] += 1
                    added.add(bp)
    
    rbp_df = pd.DataFrame\
        .from_dict(rbp_count, orient='index')\
        .sort_values(0, ascending=False) 
    rbp_df.to_csv(figure_path + '/rbp_table.tsv', sep='\t')
    
    rbp_df.head(top_n).plot.bar(ax = ax, color = ce.encoder['RBP'])
    ax.legend().set_visible(False)
    ax.set_xlabel('')#RNA-binding protein')
    ax.set_ylabel('Number of protected\nRNA binding site')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=70, rotation_mode='anchor', ha = 'right')
    sns.despine()
    color_rbp(ax)
    return rbp_df


Rfam_labs = {'RnaseP':'black',
            'tRNA': "#999999",
            'snoRNA':"#CC79A7", 
            'IsrR': "#D55E00",
            'miRNA': "#0072B2",
            'vRFE': "#F0E442",
            'Others':"#009E73", 
            'Other unannotated sncRNA':"#009E73", 
            'ToxI':"#56B4E9",
            'KRAS_3UTR':"#E69F00",
            'Hemaglobin':'red',
            'tRNA-like': '#245db7',
            'Excised structural intron RNA':'#f9bb00'}
rfam_ce = color_encoder()
rfam_ce.encoder = Rfam_labs
def group_annotation(x):
    lab = 'Other unannotated sncRNA'
    if re.search('tRNA', x):
        lab = 'tRNA'
        lab = 'tRNA-like'
#    elif re.search('RNaseP',x):
#        lab = Rfam_labs[0]
#    elif re.search('[sS][nN][oO]|[sS][nN][rR]|HACA', x):
#        lab = 'snoRNA'
#    elif x == 'IsrR':
#        lab = 'IsrR'
    elif re.search('mir|MIR', x):
        lab = 'miRNA'
        lab = 'Excised structural intron RNA'
    elif x == 'veev_FSE':
        lab = 'vRFE'
    return lab

def get_peak_rfam_annotation(peaks):
    cmscan_df = read_tbl(peak_path + '/unfragmented.Long_RNA.tblout') \
        .assign(peakname = lambda d: d['query name'].str.split('_chr', expand=True).iloc[:,0])\
        .merge(peaks.filter(['sense_gname','peakname']), on = 'peakname', how = 'right')\
        .assign(score = lambda d: d.score.fillna(0))\
        .fillna('NA')\
        .assign(strand = lambda d: np.where(d.strand=="+", 0, 1) )\
        .assign(score = lambda d: d.score.astype(float))\
        .groupby('peakname', as_index=False)\
        .apply(lambda d: d.pipe(lambda d1: d1[d1.strand==d1.strand.min()]).nlargest(1,'score'))\
        .assign(rfam_lab = lambda d: d['target name'].map(group_annotation))

    return {row['sense_gname']:row['rfam_lab'] for i, row in cmscan_df.iterrows()}
    

def pick_lp(d):
    return d \
        .pipe(lambda d: d[(d.log10p==d.log10p.max())])


    
def plot_long_RNA_peak(peaks, ax, ce, top_n = 10, y_val = 'log10p'):
    lp = peaks[peaks.sense_gtype.str.contains('Long RNA')] \
        .query('sample_count >= %i' %sample_cutoff)\
        .groupby('sense_gname', as_index=False)\
        .apply(pick_lp) \
        .nlargest(top_n, y_val)
    rfam_labs = get_peak_rfam_annotation(lp)
    rfam_labs['CASKIN2'] = 'Excised structural intron RNA'
    assert(y_val in ['log10p','pileup'])
    name_conversion = {'RP11-958N24.2': 'PKD1P4-NPIPA8',
                'RP11-1212A22.1':'NPIPA8',
                'RP11-958N24.2':'PKD1P4-NPIPA8',
                'RP11-958N24.1':'PKD1P3-NPIPA1',
                'RP11-958N24.1':'NPIPA1',
                'RP11-163G10.3':'TMSB4X',
                'AC019188.1':'TMSB4XP8',
                'AP001340.2':'RPL13AP7',
                'RP11-1212A22.4':'PKD1P5'}
    rev_name_conversion = {v:k for k,v in name_conversion.items()}

    lp = lp.assign(picked_RNA_sense = lambda d: d.sense_gname\
                .map(lambda x: name_conversion[x] if x in name_conversion.keys() else x))
    lp.plot.bar('picked_RNA_sense',
              y_val,
              color = ce.encoder['Long RNA'],
             ax = ax)
    ax.legend().set_visible(False)
    ax.set_xlabel('')
    if y_val == 'log10p':
        ax.set_ylabel('-$log_{10}$ p-value')
    else:
        ax.set_ylabel('Coverage')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=70, rotation_mode='anchor', ha = 'right')
    
    used_rfam = []
    for i, xt in enumerate(ax.get_xticklabels()):
        #if re.search('PKD1|ARHGAP30|NPIPA1|CASK', xt.get_text()): 
        #    # intron
        #    color = 'salmon'
        #elif re.search('B2M|FTH|RPL|RMRP|RPPH1|TMSB4|HBA|HBB|NRGN|PPBP|HIST|ALB|RPS' ,xt.get_text()):
        #    # full lengtj
        #    color = 'skyblue'
        #elif re.search('AC092156.2|CENPP|RP11-3N2.6|AC068137.4|RP11-193H5.8|RP11-1217F2.24|AC073869.8', xt.get_text()):
        #    color = 'purple'
        #color = 'black'
        gn = xt.get_text()
        if gn in rev_name_conversion.keys():
            gn = rev_name_conversion[gn]
        rfam = rfam_labs[gn]
        used_rfam.append(rfam)
        color = rfam_ce.encoder[rfam]
        xt.set_color(color)

    plot_ce = color_encoder()
    plot_ce.encoder = Rfam_labs.copy()
    plot_ce.encoder = {k:v for k,v in plot_ce.encoder.items() if k in used_rfam}
    plot_ce.show_legend(ax,bbox_to_anchor = (0.1,0.7), 
                        fontsize=15, frameon=False)


def plot_peak_number(peaks,ax, ce):
    for pt, pdf in peaks \
            .query('pileup >= %i' %pileup_cutoff)\
            .assign(peak_count = 1)\
            .groupby(['sense_gtype', 'pileup'], as_index=False)\
            .agg({'peak_count':'sum'}) \
            .sort_values('pileup')\
            .reset_index() \
            .assign(cum_count = lambda d: d.groupby('sense_gtype').peak_count.cumsum())\
            .assign(log_pile = lambda d: d.pileup.transform(np.log10))\
            .groupby('sense_gtype'):
        ax.plot(pdf['log_pile'],pdf['cum_count'], label = pt, color = ce.encoder[pt])
    ax.legend(bbox_to_anchor=(1,1), frameon=False,fontsize=13)
    ax.vlines(ymin=0, ymax=  1200, x = np.log10(pileup_cutoff), color = 'red')
    ax.set_ylabel('Cumulative number of peaks')
    ax.set_xlabel('Coverage (number of fragments)')
    xrange = np.arange(5) 
    ax.set_xticks(xrange)
    xt = ax.set_xticklabels(['$10^{%i}$' %i for i in xrange])
    sns.despine()
    
def plot_peak_coverage(peaks,ax, log=True):
    xrange = np.arange(5)

    xs = peaks.query('pileup >= 0').pileup
    xcut = pileup_cutoff

    if log:
        xs = np.log10(xs)
        xcut = np.log10(xcut)

    sns.distplot(xs, hist=False, ax = ax)
    ax.vlines(ymin=0, ymax=  10, x = np.log10(pileup_cutoff), color = 'red')
    ax.set_ylabel('Peaks (%)')
    ax.set_xlabel('Coverage (number of fragments)')

    if log:
        ax.set_xticks(xrange)
        ax.set_xlim(xrange[0], xrange[-1])
        xt = ax.set_xticklabels(['$10^%i$' %i for i in xrange])
    sns.despine()   
    
def plot_cov_density(peaks, ax):
    
    for strand, strand_df in peaks\
                .assign(log_pile = lambda d: np.log10(d.pileup))\
                .query('pileup >= %i' %(pileup_cutoff))\
                .groupby('is_sense'):
        sns.distplot(strand_df['log_pile'], 
                     hist=False, 
                     ax = ax,
                    label = strand)
    ax.legend(bbox_to_anchor = (0.6,0.6), title='',fontsize=15)
    
    ax.vlines(ymin = -10, 
              ymax = cum_peak.cum_count.max() + 100,
              x = np.log10(pileup_cutoff),
             color = 'red')



    xrange = np.arange(5) 
    ax.set_xticks(xrange)
    xt = ax.set_xticklabels(['$10^{%i}$' %i for i in xrange])
    ax.set_xlim(xrange.min(),xrange.max() + 1)
    ax.set_xlabel('Coverage (number of fragments)')
    ax.set_ylabel('% peaks')
    
    
def plot_peak_cum_cov(peaks, ax):
    cum_peak = peaks\
        .assign(peak_count = 1)\
        .groupby(['is_sense','pileup'], as_index=False)\
        .agg({'peak_count':'sum'}) \
        .assign(peak_count = lambda d: np.where(d.pileup <= pileup_cutoff,
                                 0,
                                 d.peak_count))\
        .sort_values('pileup') \
        .assign(cum_count = lambda d: d.groupby('is_sense')['peak_count'].cumsum())\
        .assign(cum_count = lambda d: d.groupby('is_sense')['cum_count'].transform(lambda x: x/x.max()))

    for strand, sdf in cum_peak.groupby('is_sense'):
        ax.plot(sdf['pileup'], sdf['cum_count'], label=strand)
    ax.set_xscale('log')
    ax.set_xlim(0, 1e4)
    ax.set_xlabel('')
    ax.set_ylabel('ECDF')
    ax.legend(title = '', fontsize=15, frameon = False)
    sns.despine()
    ax.set_xlabel('Coverage (number of fragments)')


def plot_peak_size(peak_df, ax):
    for strand, strand_df in peak_df\
            .assign(psize = lambda d: d.end - d.start)\
            .groupby('sense_gtype'):
        sns.distplot(strand_df.psize, hist=False, label = strand, ax = ax)
    ax.set_xlabel('Peak size')
    ax.set_ylabel('Density')
    




def is_mt(seq, rnr=False):
    is_chrM = 'not_MT'
    chrom_path = '/stor/work/Lambowitz/ref/hg19'
    if rnr:
        genome = chrom_path + '/new_genes/mt_rnr.fa'
    else:
        genome = chrom_path + '/genome/chrM.minimap2_idx'

    aligner = mappy.Aligner(genome,preset='sr')
    if list(aligner.map(seq)):
        is_chrM = 'is_MT'
    return is_chrM

fa = pysam.Fastafile('/stor/work/Lambowitz/ref/hg19_ref/genome/hg19_genome.fa')
def fetch_seq(chrom, start, end, strand):
    seq = fa.fetch(chrom, int(start), int(end))
    return seq if strand == "+" else reverse_complement(seq)


ce = color_encoder()
colors = simpsons_palette()
ce.encoder = {
    'Long RNA': '#370335',
     'RBP': '#91331F',
     'Repeats': '#197EC0',
     'Unannotated': '#46732E',
     'miRNA': '#FD7446',
     'misc RNA': '#FD8CC1',
     'tRF3':'black',
     'tRF5':'black',
     '.':'black',
     'piRNA': '#D5E4A2',
     'snRNA': '#8A9197',
     'snoRNA': '#FED439'
}


HB_genes = '''
chr7,106809406,106842974,HBP1
chr11,5269309,5271089,HBG1
chr11,5274420,5526835,HBG2
chr11,5289575,5526882,HBE1
chr16,222875,223709,HBA2
chr16,226679,227521,HBA1
chr16,203891,216767,HBM
chr16,230452,231180,HBQ1'''
HB_genes = pd.read_csv(io.StringIO(HB_genes),
                      names = ['chrom','start', 'end', 'HB'])
def is_hb(row):
    answer = 'Not HB'
    if row['chrom'] in HB_genes.chrom.tolist():
        hb_chrom = HB_genes.query('chrom =="%s"' %row['chrom'])
        if any(max(hb_row['start'], row['start']) <= min(hb_row['end'],row['end']) for i, hb_row in hb_chrom.iterrows()):
            answer = 'HB'
    return answer


def anti_tblout():
    tblout = read_tbl(peak_path + '/unfragmented.others.tblout') \
            .groupby('query name', as_index=False)\
            .apply(lambda d: d[d.score == d.score.max()])\
            .filter(regex='name') \
            .rename(columns = {'query name':'peakname',
                                'target name':'rfam'})
    if tblout.shape[0] != 0:
        return tblout.assign(peakname = lambda d: d.peakname.str.split('_chr',expand=True).iloc[:,0])


def plot_anti_bar(antisense_peaks, ax):
    tblout = anti_tblout()
    anti_plot = antisense_peaks.nlargest(15, 'log10p')\
        .assign(antisense_gname = lambda d: np.where(d.antisense_gname == ".",
                                                    d.chrom + ':' + d.start.astype(str) + '-' + d.end.astype(str),
                                                    d.antisense_gname))\
        .assign(is_hb = lambda d: [is_hb(row) for i, row in d.iterrows()])\
    
    if tblout:
        anti_plot = anti_plot.merge(tblout,
                on = 'peakname', how = 'left')
    else:
        anti_plot = anti_plot.assign(rfam = None)
    
    anti_plot = anti_plot.assign(rfam = lambda d: d.rfam.fillna('Other unannotated sncRNA'))\
            .assign(rfam = lambda d: np.where(d.is_hb=="HB", 'Hemaglobin', d.rfam))\
            .assign(rfam = lambda d: np.where(d.rfam=="HBM", 'Hemaglobin', d.rfam))\
            .assign(rfam = lambda d: np.where(d.rfam=="FHbp_thermometer", 'Other unannotated sncRNA', d.rfam))\
            .sort_values('log10p', ascending=False)


    anti_plot\
        .plot\
        .bar('antisense_gname', 'log10p', 
            color = anti_plot.antisense_gtype.map(ce.encoder),
            ax = ax)
    ax.legend().set_visible(False)
    ax.set_xlabel(' ')
    ax.set_ylabel('-$log_{10}$ p-value')
    ax.set_xticklabels(ax.get_xticklabels(), 
                    rotation = 70,
                    rotation_mode='anchor',
                    ha = 'right', va = 'center')
    

    used_rfam = []
    for xt, rfam in zip(ax.get_xticklabels(), anti_plot.rfam):
        xt.set_color(rfam_ce.encoder[rfam])
        used_rfam.append(rfam)
    
    plot_ce = color_encoder()
    plot_ce.encoder = Rfam_labs.copy()
    plot_ce.encoder = {k:v for k,v in plot_ce.encoder.items() if k in used_rfam}
    plot_ce.show_legend(ax, frameon=False, fontsize=15,
                        bbox_to_anchor=(-0.1,0))



class ecoli_mapper():
    def __init__(self):
        bam = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/merged_bam/dedup/unfragmented.chrM_filter.bam'
        index = '/stor/work/Lambowitz/ref/Ecoli/BL21_DE3.fa'
        self.bam = pysam.Samfile(bam)
        self.aligner = BwaAligner(index, options = '-k 12')
        self.matched = re.compile('([0-9]+)M')
        self.clipped = re.compile('([0-9]+)S')
        self.alignments = None

    def ecoli_map(self, chrom, start, end):
        aligned = 0.0
        self.alignments = []
        for aln_count, aln in enumerate(self.bam.fetch(chrom, start, end)):
            alns = self.aligner.align_seq(aln.query_sequence)
            self.alignments.append(alns)
            filtered_alignments = filter(self.filter_bad_cigar, alns)
            if list(filtered_alignments) :
                aligned += 1
        return aligned / (aln_count + 1)


    def filter_bad_cigar(self, aln):
        clipped_base = sum(map(int, self.clipped.findall(aln.cigar))) or 0
        mapped_base = sum(map(int, self.matched.findall(aln.cigar)))
        return (float(clipped_base) / mapped_base) < 0.2  and aln.NM < 3



