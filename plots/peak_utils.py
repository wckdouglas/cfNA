import pandas as pd
import sys
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
from transcriptome_filter import peak_analyzer
import pyximport 
pyximport.install()
from junction_func import get_junction 
import matplotlib.pyplot as plt
sys.path.insert(0,'/stor/home/cdw2854/cfNA/peak_callings')
from structural_peaks import PeakAnalyzer, mRNAFilter, GenicIntersect, NameConversion, GeneMapper, TrnaLookAlike
from exon_coverage import ExonFilter
import dask.dataframe as dd
plt.rc('axes', labelsize=20)
plt.rc('xtick', labelsize = 20)
plt.rc('ytick', labelsize = 20)
plt.rc('font', **{'family':'sans-serif',
                  'sans-serif':'Arial'})



pileup_cutoff = 5
sample_cutoff = 5
project_path = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
project_path = '/stor/work/Lambowitz/yaojun/Work/cfNA/tgirt_map'
peak_path = project_path + '/bed_files/merged_bed/MACS2/annotated'
peak_type_ce = color_encoder()
peak_type_ce.encoder = {'mRNA':'purple',
                 'Pseudogene':'darkblue',
                 'Exon': 'purple',
                 'Intron':'#fca237',
                 'Exon-intron':'#7bb73e',
                 'Within intron':'#f9b109', 
                 'Stem-loop':'#f9b109', 
                 'miRNA':'darkgreen',
                 'rRNA':'#15879b',
                 'Mismapped':'#bcbb76',
                 'Others':'black',
                 'Intergenic':'black',
                 'tRNA-lookalike': 'red',
                 'Full-length intron':'#725001',
                 'RBP':'#91331F',
                 'Excised full-length intron':'#725001'}
gene_mapper = GeneMapper()


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


def load_peaks(filename):
    EPSILON = 1e-100
    peak_df = pd.read_csv(filename, sep='\t') 
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
    pie_df = peaks\
        .query('pileup>=%i & sample_count >= %i' %(pileup_cutoff, sample_cutoff))\
        .groupby('is_sense')\
        .agg({'sense_gtype':'count'})\
        .reset_index() \
        .assign(index = lambda d: d['is_sense'] + '\n(' + d.sense_gtype.astype(str)+')')\
        .set_index('index')\
        .sort_values('sense_gtype', ascending=True)
    
    
    wedges, texts = ax.pie(pie_df.sense_gtype, 
                explode = [0,0.05, 0.2], 
                startangle = 180, 
                colors = ['#2a7a1c','#40649e','#f7b707'])
    ax.set_ylabel('')
    ax.legend().set_visible(False)

    yoffset = [-4,1,1]
    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        ax.annotate(pie_df.index[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y*yoffset[i]),
                     horizontalalignment=horizontalalignment, fontsize=20,
                     arrowprops={'arrowstyle':'-'})


def change_annotation(lab):
    if 'RBP' in lab:
        return lab.replace('RBP','Long RNA\n(RBP-binding sites)')
    elif 'Long RNA' in lab:
        return lab.replace('Long RNA', 'Long RNA\n(Other narrow peaks)')
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
              labeldistance=1.15, 
              colors = colors)
    
    index = peak_pie.index
    ax.legend(bbox_to_anchor = (1.5,0.7), 
              labels = list(map(lambda x: x.split(' ')[0], index)), 
              ncol=2, 
              fontsize=20)
    
    
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
    plot_df = peaks\
        .query('sense_gtype == "%s"' %rnatype)\
        .query('pileup >= %i & sample_count >= %i' %(pileup_cutoff, sample_cutoff)) \
        .assign(RNA_name = lambda d: d.sense_gname.map(assign_rna_name)) \
        .groupby('RNA_name', as_index=False)\
        .agg({'chrom':'count'}) \
        .nlargest(top_n, 'chrom')\
        .assign(color = lambda d: d.RNA_name.map(repeat_color))
    sns.barplot(data=plot_df,x='RNA_name',y='chrom', palette = plot_df.color, ax = ax)
    ax.set_xlabel('')
    ax.set_ylabel('Peak count')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=70, rotation_mode='anchor', ha = 'right')
    ax.legend().set_visible(False)


def color_rbp(x):
    '''
    example: ','.join(rbp_df.head(15).index)
    http://plasmaproteomedatabase.org/
    '''
    color = 'black'
    if x not in ['IGF2BP1','LARP4','LIN28B']:
        color = 'red'
    return color


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
    
    rbp_df = rbp_df.head(top_n)
    colors = list(map(color_rbp, rbp_df.index.values))
    sns.barplot(rbp_df.index, rbp_df[0], palette = colors, ax=ax)
    ax.legend().set_visible(False)
    ax.set_xlabel('')#RNA-binding protein')
    ax.set_ylabel('Number of protected\nRNA binding site')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=70, rotation_mode='anchor', ha = 'right')
    sns.despine()
    for xt in ax.get_xticklabels():
        color = color_rbp(xt.get_text())
        xt.set_color(color)
    return rbp_df


Rfam_labs = {'RnaseP':'black',
            'tRNA': "#efa002",
            'snoRNA':"#CC79A7", 
            'IsrR': "#D55E00",
            'miRNA': "#0072B2",
            'vRFE': "#F0E442",
            'Others':"#009E73", 
            'Unannotated sncRNA':"#009E73", 
            'ToxI':"#56B4E9",
            'KRAS_3UTR':"#E69F00",
            'Hemoglobin':'red',
            'tRNA-lookalike': '#ad1b34',
            'rRNA':'#030544',
            'miRNA-like':"#0072B2",
            'Pseudogene':'#f4162b',
            'Excised structured intron RNA':'#f78d02'}
rfam_ce = color_encoder()
rfam_ce.encoder = Rfam_labs
def group_annotation(x):
    lab = 'Unannotated sncRNA'
    if re.search('tRNA', x):
        lab = 'tRNA'
        lab = 'tRNA-lookalike'
#    elif re.search('RNaseP',x):
#        lab = Rfam_labs[0]
#    elif re.search('[sS][nN][oO]|[sS][nN][rR]|HACA', x):
#        lab = 'snoRNA'
#    elif x == 'IsrR':
#        lab = 'IsrR'
    elif re.search('mir|MIR', x):
        lab = 'miRNA'
        lab = 'Excised structured intron RNA'
    elif x == 'veev_FSE':
        lab = 'vRFE'
    return lab

def get_peak_rfam_annotation(peaks):
    cmscan_df = read_tbl(peak_path + '/unfragmented.tblout') \
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


def long_rna_type(x):
    if x in {'AB019441.29','RPS2P55','RPL13AP25'}:
        rt = 'Pseudogene'
    elif x in 'DAPK1':
        rt = 'miRNA'
    elif re.search('PKD|ARHG|CASK',x):
        rt = 'Full-length intron'
    elif re.search('CPN1|CACNA', x):
        rt = 'tRNA-lookalike'
    elif re.search('PF4', x):
        rt = 'Exon-intron'
    else:
        rt = 'Within intron'
    return rt


trna = TrnaLookAlike()
def trna_lookalike(row):
    trnalookalike = trna.search(row['chrom'], row['start'], row['end'], row['strand'])
    if trnalookalike != ".":
        return 'tRNA-lookalike'
    else:
        return row['rt']


def cat_long_rna_type(d):
    gi = GenicIntersect()
    return d \
        .assign(rt = lambda d: d.picked_RNA_sense.map(long_rna_type))\
        .assign(rt = lambda d: [trna_lookalike(row) for i, row in d.iterrows()])\
        .pipe(gi.fulllength_intron)\
        .assign(rt = lambda d: np.where(d['fulllength_intron']!='.',
                                        'Full-length intron',
                                        d.rt))


    
def plot_long_RNA_peak(peaks, ax, ce, top_n = 10, y_val = 'log10p'):
    lp = peaks[peaks.sense_gtype.str.contains('Long RNA')] \
        .query('sample_count >= %i' %sample_cutoff)\
        .groupby('sense_gname', as_index=False)\
        .apply(pick_lp) 
    rfam_labs =  defaultdict(lambda: 'Others') #get_peak_rfam_annotation(lp)
    rfam_labs['CPN1'] = 'tRNA-lookalike'
    rfam_labs['CASKIN2'] = 'Excised structured intron RNA'
    rfam_labs['DAPK1'] = 'miRNA-like'
    rfam_labs['RP11-51O6.1'] = 'Pseudogene'

    assert(y_val in ['log10p','pileup'])
    name_conversion = NameConversion()
    rev_name_conversion = {v:k for k,v in name_conversion.encoder.items()}

    lp = lp\
        .assign(picked_RNA_sense = lambda d: d.sense_gname.map(name_conversion.convert).str.replace('-NPIPA8','')) \
        .groupby('picked_RNA_sense')\
        .apply(lambda d: d.nlargest(1, y_val))\
        .nlargest(top_n, y_val) \
        .pipe(cat_long_rna_type)\
        .sort_values(y_val, ascending=False)
    colors = lp.rt.map(peak_type_ce.encoder).values
    sns.barplot(data=lp,
                x='picked_RNA_sense',
                y=y_val,
                palette = colors,
                ax = ax)
    ax.legend().set_visible(False)
    ax.set_xlabel('')
    if y_val == 'log10p':
        ax.set_ylabel('-$log_{10}$ p-value', fontsize=20)
    else:
        ax.set_ylabel('Coverage', fontsize=20)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=70, rotation_mode='anchor', ha = 'right')
    
    used_rfam = []
    for i, xt in enumerate(ax.get_xticklabels()):
        gn = xt.get_text()
        if gn in rev_name_conversion.keys():
            gn = rev_name_conversion[gn]
        rfam = rfam_labs[gn]
        used_rfam.append(rfam)

    used = lp.rt.unique()
    cc_ce = color_encoder()
    cc_ce.encoder = {k:v for k,v in peak_type_ce.encoder.items() if k in used}
    cc_ce.show_legend(ax = ax, frameon=False, fontsize=20)

    for col,xt in  zip(colors,ax.get_xticklabels()):
        xt.set_color(col)



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
    ax.vlines(ymin=0, ymax=  8, x = np.log10(pileup_cutoff), color = 'red')
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
        sns.distplot(strand_df.psize, hist=False, label = strand, ax = ax, color = ce.encoder[strand])
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





def anti_tblout():
    tblout = read_tbl(peak_path + '/unfragmented.tblout') \
            .query('strand == "+"')\
            .pipe(lambda d: d[d['E-value']< 0.01])\
            .groupby('query name', as_index=False)\
            .apply(lambda d: d[d.score == d.score.max()])\
            .filter(regex='name') \
            .rename(columns = {'query name':'peakname',
                                'target name':'rfam'})
    if tblout.shape[0] != 0:
        return tblout.assign(peakname = lambda d: d.peakname.str.split('_chr',expand=True).iloc[:,0])

def rename_hb(row):
    if row['hb'] == 'HB':
        gn = gene_mapper.test_gene(row['chrom'], row['start'], row['end'], return_name=True) +\
            '\n(' + row['chrom'] + \
            ':' + str(row['start']) + \
            '-' + str(row['end']) + ')'
    else:
        gn = row['antisense_gname']
    return gn




def plot_anti_bar(antisense_peaks, ax, bbox = (1.2,-0.3)):
    tblout = anti_tblout()
    anti_plot = antisense_peaks.nlargest(30, 'log10p')\
        .assign(antisense_gname = lambda d: np.where(d.antisense_gname == ".",
                                                    d.chrom + ':' + d.start.astype(str) + '-' + d.end.astype(str),
                                                    d.antisense_gname))\
        .assign(hb = lambda d: [gene_mapper.test_gene(chrom, start, end, return_name=False) for chrom,start,end in zip(d.chrom, d.start, d.end)] )
    
    if tblout is not None:
        anti_plot = anti_plot.merge(tblout,
                on = 'peakname', how = 'left')
    else:
        anti_plot = anti_plot.assign(rfam = None)
    
    anti_plot = anti_plot.assign(rfam = lambda d: d.rfam.fillna('Unannotated sncRNA'))\
            .assign(rfam = lambda d: np.where(d.hb=="HB", 'Hemoglobin', d.rfam))\
            .assign(antisense_gname = lambda d: [rename_hb(row) for i, row in d.iterrows()])\
            .assign(rfam = lambda d: np.where(d.rfam=="HBM", 'Hemoglobin', d.rfam))\
            .assign(rfam = lambda d: np.where(d.rfam=="FHbp_thermometer", 'Unannotated sncRNA', d.rfam))\
            .assign(rfam = lambda d: np.where((d.chrom == 'chr13') & (d.start > 57262600) & (d.end < 57262700),
                                                'rRNA',
                                                d.rfam ))\
            .assign(antisense_gname = lambda d: np.where((d.chrom == 'chrX') & (d.start == 12994906),
                                                        'TMSB4X\n(' + d.chrom + ':' + d.start.astype(str) + '-' + d.end.astype(str) + ')',
                                                        d.antisense_gname)) \
            .sort_values('log10p', ascending=False)
    
    if any('HBQ1' in x for x in anti_plot.antisense_gname.tolist()):
        anti_plot = anti_plot.pipe(lambda d: d[d.antisense_gname.str.contains('^HB|TMSB4X')])
    #    anti_plot = anti_plot.pipe(lambda d: d[~d.antisense_gname.str.contains('TLE|SCHLAP1')])


    anti_plot\
        .plot\
        .bar('antisense_gname', 'log10p', 
            color = 'steelblue',
#            color = anti_plot.antisense_gtype.map(ce.encoder),
            ax = ax)
    ax.legend().set_visible(False)
    ax.set_xlabel(' ')
    ax.set_ylabel('-$log_{10}$ p-value')
    ax.set_xticklabels(ax.get_xticklabels(), 
                    rotation = 70,
                    rotation_mode='anchor',
                    ha = 'right', va = 'center')
    

    #used_rfam = []
    #for xt, rfam in zip(ax.get_xticklabels(), anti_plot.rfam):
    #    xt.set_color(rfam_ce.encoder[rfam])
    #    used_rfam.append(rfam)
    
    #plot_ce = color_encoder()
    #plot_ce.encoder = Rfam_labs.copy()
    #plot_ce.encoder = {k:v for k,v in plot_ce.encoder.items() if k in used_rfam}
    #plot_ce.show_legend(ax, frameon=False, fontsize=15,
    #                    bbox_to_anchor=bbox)



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



PEAK_ANALYZER = peak_analyzer(project_path + '/merged_bam/dedup/unfragmented.chrM_filter.dedup.bam',
                    '/stor/work/Lambowitz/ref/hg19_ref/genes/transcriptome.minimap2_idx')

def transcriptome_map(chrom, start, end, strand):
    mapped, num_pairs, transcript = PEAK_ANALYZER.filter_alignments(chrom, int(start), int(end), strand)
    return mapped/num_pairs




class mRNAFilter():
    '''
    if the peak is on exon?
    is the peak also called in transcriptome?
    '''
    def __init__(self):
        ref_path = '/stor/work/Lambowitz/ref/hg19_ref/genes'
        exons = ref_path + '/gencode.exon.bed.gz'
        self.exons = pysam.Tabixfile(exons)
        transcriptom_peaks = project_path + '/transcriptome/macs2/unfragmented.fwd_peaks_genomics.narrowPeak.gz'
        self.transcriptome_peaks = pysam.Tabixfile(transcriptom_peaks)
        self.bam = pysam.Samfile(project_path + '/merged_bam/dedup/unfragmented.chrM_filter.dedup.bam')
        self.bed = pysam.Tabixfile(project_path + '/bed_files/merged_bed/unfragmented.bed.gz')

    def search(self, chrom, start, end, attribute = 'exon'):
        if attribute == 'exon':
            it = self.exons
        elif attribute == 'transcriptome':
            it = self.transcriptome_peaks
        return 'yes' if any(it.fetch(chrom, start, end)) else 'no'


    def spliced(self, chrom, start, end):
        spliced = 0
        for read_count, read in enumerate(self.bam.fetch(chrom, start, end)):
            if 'N' in read.cigarstring:
                spliced += 1
        return spliced/(read_count+1)


    def fragment_test(self, chrom, start, end, strand):
        frag_count = 0
        fulllength = 0
        for frag in self.bed.fetch(chrom, start, end):
            fields = frag.split('\t')
            frag_strand = fields[5]
            if frag_strand == strand:
                frag_count += 1
                if start -5 < int(fields[1]) < start + 5 and end -5 < int(fields[2]) < end + 5:
                    fulllength += 1
        
        if frag_count == 0:
            return 0
        return fulllength



def long_rna_df():
    gi = GenicIntersect()
    peak_file = peak_path + '/unfragmented.filtered.tsv'
    bed = load_peaks(peak_file) \
        .query('sample_count >= %i & pileup >= %i' %(sample_cutoff, pileup_cutoff))\
        .query('sense_gtype == "Long RNA"') \
        .assign(psi = lambda d: list(map(gi.compute_psi, d.chrom, d.start, d.end, d.strand))) 
    columns = bed.columns.tolist()
    gbed = bed\
        .pipe(gi.intersect)\
        .pipe(gi.resolve_genic)
    return gbed

def exonic_filtered_df(peak_df):
    exon_filter = ExonFilter()
    needed_columns = ['chrom','start','end','peakname','score','strand']
    needed_columns.extend(list(set(peak_df.columns.tolist()) - set(needed_columns)))
    exon_filtered_peaks = peak_df\
        .filter(needed_columns)
    exon_filtered_peaks = exon_filter.filter(exon_filtered_peaks, f = 0.1)
    return exon_filtered_peaks


def repeat_color(x):
    color = 'black'
    if not re.search('\)n$|rich$', x):
        color=  'red'
    return color


def plot_repeat_peaks(ax):
#    project_path = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
    peak_path = project_path + '/bed_files/merged_bed/MACS2/annotated'
    #peak_path = project_path + '/CLAM//BED_files/peaks/annotation'
    peak_tsv = peak_path + '/unfragmented.filtered.tsv'
    peak_df = load_peaks(peak_tsv)  \
        .assign(sense_gtype = lambda d: np.where(d.sense_gtype == ".", 'Unannotated', d.sense_gtype))\
        .assign(antisense_gtype = lambda d: np.where(d.antisense_gtype == ".", 'Unannotated', d.antisense_gtype)) \
        .sort_values('pileup', ascending=False)  
    sense_peaks = peak_df.query('is_sense == "Sense"')
    plot_repeats_RNA(sense_peaks, ax, ce, rnatype='Repeats', top_n = 15)
    for xt in ax.get_xticklabels():
        color = repeat_color(xt.get_text())
        xt.set_color(color)


def read_peak_type():
    return pd.read_csv(project_path + '/bed_files/merged_bed/MACS2/annotated/unfragmented.filtered.tsv',
                        sep='\t', usecols = [0,1,2,11]) \
            .assign(coordinate = lambda d: d.chrom + ':' + d.start.astype(str) + '-' + d.end.astype(str)) \
            .drop(['start','chrom','end'], axis=1) \
            .rename(columns={'sense_gtype':'gtype'})
            

        
