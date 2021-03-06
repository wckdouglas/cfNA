from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches 
import seaborn as sns
import pandas as pd
import numpy as np
from tgirt_map.table_tools import change_gene_type
from collections import defaultdict
from sequencing_tools.viz_tools import okabeito_palette, \
                        simpsons_palette, \
                        color_encoder
from sequencing_tools.io_tools import ReadPicardRNA
from collections import defaultdict
import re
import glob
import os
from plotting_utils import label_sample, rename_sample, \
                        label_ce, rna_type_ce, \
                        figure_path, work_path
from functools import lru_cache
plt.rc('font', **{'family':'sans-serif',
                    'sans-serif':'Arial'})

small_RNA_ce = color_encoder()
label_order = ['Untreated','NaOH', 'WGS-sim', 'DNase I', 'DNase I + Exo I',"DNase I - 3'P"]


metric_path = work_path + '/cfNA/tgirt_map/merged_bam/filtered_bam'
metrics = glob.glob(metric_path + '/*.RNA_Metrics')
metrics = list(filter(lambda x: 'sense' not in x, metrics))
def read_metric(metric):
    return pd.read_table(metric, skiprows=6, nrows=1)\
        .pipe(pd.melt) \
        .pipe(lambda d: d[d.variable.str.contains('TRANSCRIPT_STRAND_')])\
        .pipe(lambda d: d[d.variable.str.contains('PCT')]) 

def plot_strand(ax):
    strand_df = {os.path.basename(metric).split('.')[0]: read_metric(metric) for metric in metrics}
    strand_df = pd.concat([d.assign(samplename = k) for k, d in strand_df.items()])\
        .assign(samplename = lambda d: d.samplename.map(label_sample))\
        .assign(variable = lambda d: np.where(d.variable.str.contains('R1'), 'Sense','Antisense'))\
        .assign(value = lambda d: d['value'] * 100)\
        .pipe(pd.pivot_table, index='samplename', columns = 'variable', values = 'value')\
        .reindex(index=label_order)\
        .plot.bar(stacked=True, ax = ax, width = 0.8)
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 70, ha = 'right',rotation_mode="anchor")
    ax.legend()
    ax.set_xlabel('')
    ax.set_ylabel('Protein-coding reads (%)')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], 
              bbox_to_anchor = (.9,1), fontsize = 20,
             frameon=False)
    
    
def plot_coding_bases(ax):
    RNA_base_from_picard = ReadPicardRNA(attr = 'base')
    RNA_base_from_picard.read(metrics) \
        .assign(var_count = lambda d: d.var_count*100)\
        .assign(samplename = lambda d: d.samplename.str.split('.',expand=True).iloc[:,0].map(label_sample))\
        .assign(variable = lambda d: d.variable.str.replace('Utr','UTR'))\
        .assign(variable = lambda d: d.variable.str.replace(' bases',''))\
        .pipe(pd.pivot_table, columns = 'variable', index='samplename', values = 'var_count')\
        .reindex(index=label_order)\
        .plot.bar(stacked=True, width = 0.8, color = simpsons_palette(), ax = ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 70, ha = 'right',rotation_mode="anchor")
    ax.set_ylabel('Protein-coding bases (%)')
    ax.set_xlabel('')
    ax.legend(title = '', bbox_to_anchor = (1,1))
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], bbox_to_anchor = (0.9,1), 
              fontsize = 20, frameon=False)


def plot_insert(ax, samples=['DNase I']):
    insert_path = work_path + '/cfNA/tgirt_map/fragment_sizes'
    data_files = glob.glob(insert_path + '/*.feather')
    df = {os.path.basename(data_file):pd.read_feather(data_file) for data_file in data_files}
    df = pd.concat([val.assign(label = key) for key, val in df.items()]) \
        .assign(label = lambda d: d.label.str.replace('.feather','').str.capitalize())\
        .assign(label = lambda d: np.where(d.label == 'Polya','PolyA-selected', d.label))\
        .groupby(['isize','label','bed'], as_index=False)\
        .agg({'size_count':'sum'})\
        .sort_values('isize') \
        .query('isize < 300') \
        .assign(size_fraction = lambda d: d.groupby(['bed','label'])['size_count'].transform(lambda x: 100* x/ x.sum()))\
        .groupby(['label','isize'], as_index=False)\
        .agg({'size_fraction':'median'})\
        .assign(label = lambda d: d.label.map(label_sample)) \
        .pipe(lambda d: d[d.label.isin(samples)])

#        .pipe(lambda d: d[d.label.str.contains('Alk|Un|[Ee]xo|All')]) \
    
    for lab, lab_df in df.groupby('label'):
        ax.plot(lab_df['isize'], 
                 lab_df['size_fraction'], 
                 linewidth=2,
                 label = lab,
                 alpha=0.7,
                 color = label_ce.encoder[lab])
        
    ax.legend(title= ' ', 
             fontsize = 20, 
             frameon=False, 
             bbox_to_anchor = (0.65,0.8))
    ax.set_xlabel('Read span (nt)')
    ax.set_ylabel('Read pairs (%)')
    
    
    for i, label in enumerate(ax.get_xticklabels()):
        label.set_x(label.get_position()[0] + 1)
    
    return df

def recat_rRNA(gname, gtype):
    if '18S' in gname or '28S' in gname:
        return '18/28S rRNA'
    elif '5.8S' in gname or '5S' in gname or '5_8S' in gname:
        return '5/5.8S rRNA'
    else:
        return gtype

@lru_cache(maxsize=128, typed=False)
def read_count(feature_only=True, dedup=True, rna_group_type = 'grouped_type', col_regex='type|Qcf|QCF|sim'):
    count_file = work_path + '/cfNA/tgirt_map/Counts/all_counts/spreaded_all_counts.feather'
    dedup_df = pd.read_feather(count_file)

    filter_feature = 'Unannotated' if feature_only else ''
    dedup_regex = ':dedup:' if dedup else ':all:'
    countplot_df = dedup_df \
        .assign(grouped_type = lambda d: np.where(d.gene_name.str.contains('^MT-'),'Mt', d[rna_group_type]))\
        .assign(grouped_type = lambda d: np.where(d.gene_name.str.contains('^MT-T'),'Mt-tRNA', d[rna_group_type]))\
        .filter(regex = col_regex)\
        .assign(grouped_type = lambda d: np.where(d[rna_group_type] == "No features", 'Unannotated', d[rna_group_type]))\
        .assign(grouped_type = lambda d: np.where(d[rna_group_type] == "rDNA", 'rRNA', d[rna_group_type]))\
        .assign(grouped_type = lambda d: np.where(d[rna_group_type].str.contains('snoRNA|Y-RNA'), 'Other sncRNA', d[rna_group_type]))\
        .assign(grouped_type = lambda d: np.where(d[rna_group_type].str.contains("vaultRNA|VT|vt"),'Vault RNA', d[rna_group_type]))\
        .assign(grouped_type = lambda d: list(map(recat_rRNA, d.gene_type, d.grouped_type)))\
        .query('grouped_type != "rRNA"')\
        .groupby(rna_group_type)\
        .sum() \
        .pipe(lambda d: d[d.columns[~d.columns.str.contains('anti')]])\
        .pipe(lambda d: d[d.columns[d.columns.str.contains(dedup_regex)]]) \
        .reset_index()\
        .pipe(pd.melt, id_vars = [rna_group_type]) \
        .assign(variable = lambda d: d.variable.str.split(':', expand=True).iloc[:,0])\
        .assign(treatment = lambda d: d.variable.map(label_sample)) \
        .groupby([rna_group_type,'treatment'], as_index=False)\
        .agg({'value':'median'}) \
        .query('%s != "%s"' %(rna_group_type, filter_feature))\
        .assign(value = lambda d: d.groupby('treatment')['value'].transform(lambda x: 100*x/x.sum()))\
        .pipe(pd.pivot_table, index = 'treatment', 
             columns = rna_group_type,
             values = 'value')\
        .reset_index() \
        .sort_values('treatment')\
        .set_index('treatment')\
        .pipe(lambda d: d.reindex(sorted(d.columns), axis=1))
    return countplot_df

def rename_rRNA(x):
    if '5S' in x:
        return '5S_rRNA'
    elif '5-8S' in x:
        return '5.8S_rRNA'
    else:
        return x

def plot_small_count_bar(ax, prep_regex = 'DNase|WGS-sim|NaOH|Untreated', label_order = label_order):
    #.pipe(lambda d: d[d.grouped_type.str.contains('sncRNA|snoRNA|tRNA|miRNA|rRNA|rDNA')])\
    count_file = work_path + '/cfNA/tgirt_map/Counts/all_counts/all_counts.feather'
    small_df = pd.read_feather(count_file) \
        .query("dedup == 'dedup'")\
        .query('gene_type != "No features"')\
        .assign(grouped_type = lambda d: d.gene_type.map(change_gene_type)) \
        .pipe(lambda d: d[d.grouped_type.str.contains('Other sncRNA')])\
        .pipe(lambda d: d[~d.gene_id.str.contains('^[21]8S')])\
        .assign(gene_id = lambda d: np.where(d.grouped_type.str.contains('rRNA'), 
                                            d.gene_id.map(rename_rRNA),
                                            d.gene_id))\
        .assign(gene_type = lambda d: np.where(d.grouped_type.str.contains('rRNA'), 
                                            d.gene_id.map(rename_rRNA),
                                            d.gene_type))\
        .assign(gene_type = lambda d: np.where(d.gene_type.str.contains('Y-RNA'), 
                                            'Y RNA',
                                            d.gene_type))\
        .assign(gene_id = lambda d: np.where(d.gene_type=="Mt_tRNA",d.gene_name, d.gene_id))\
        .assign(gene_type = lambda d: np.where((d.gene_type=="tRNA") & (d.gene_id.str.contains('^MT')),
                                                'Mt_tRNA',
                                                d.gene_type))\
        .assign(gene_type = lambda d: np.where(d.gene_type=="vaultRNA",
                                                'vault RNA',
                                                d.gene_type))\
        .assign(gene_type = lambda d: d.gene_type.str.replace('^tRNA','Nucleo-tRNA'))\
        .pipe(lambda d: d[~d.gene_type.str.contains('^ENSG')])\
        .pipe(lambda d: d[~d.gene_id.str.contains('tRNA')])\
        .assign(treatment = lambda d: d.samplename.map(label_sample))\
        .assign(gene_id = lambda d: d.gene_id.str.replace('-[0-9]+-[0-9]+',''))\
        .groupby(['treatment','gene_id', 'gene_type'], as_index=False)\
        .agg({'read_count':'sum'})\
        .groupby(['treatment','gene_type'], as_index=False)\
        .agg({'read_count':'sum',
            'gene_id':'count'})\
        .assign(value = lambda d: d.groupby('treatment').read_count.transform(lambda x: 100*x/x.sum()))\
        .assign(gene_type = lambda d: d.gene_type.str.replace('Mt_','MT-').str.replace('_',' '))\
        .pipe(lambda d: d[d.treatment.str.contains(prep_regex)])\
        .pipe(pd.pivot_table, index=['treatment'], columns = 'gene_type', values = 'value')\
        .reindex(label_order)

    #(7sk,7sl,Y-RNA, miscRNA,piRNA,snRNA,vaultRNA)
    small_RNA_color= ['#f58231','#e6194b', '#3cb44b', '#ffe119', '#ba0ace', '#6a3d9a', 
                        '#03A8FB', '#F8BF6C', '#CAF5CB', '#fabebe', 
                        '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#F958FC', 
                        '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000']
    small_RNA_ce.encoder = {t:c for t, c in zip(small_df.columns, small_RNA_color)}
    small_df\
        .plot.bar(stacked=True,
                ax = ax, 
                color = small_RNA_color,
                width = 0.8)
    xt = [xt.set_text(xt.get_text()) for xt in ax.get_xticklabels()]
    ax.set_xticklabels(ax.get_xticklabels(), 
                       rotation = 70, 
                       ha = 'left', va = 'top',
                      rotation_mode="anchor")
    ax.legend()
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], 
              bbox_to_anchor = (0.9,1), fontsize = 20,
             frameon=False)
    ax.set_xlabel('')
    ax.set_ylabel('Small RNA read pairs (%)')
    sns.despine()
    return small_RNA_ce


 
def plot_small_count_pie(ax):
    count_file = work_path + '/cfNA/tgirt_map/Counts/all_counts/all_counts.feather'
    small_df = pd.read_feather(count_file) \
        .query("dedup == 'dedup'")\
        .query('gene_type != "No features"')\
        .assign(grouped_type = lambda d: d.gene_type.map(change_gene_type)) \
        .pipe(lambda d: d[d.grouped_type.str.contains('sncRNA|snoRNA|tRNA|miRNA|rRNA|rDNA')])\
        .pipe(lambda d: d[~d.gene_id.str.contains('^[21]8S')])\
        .assign(gene_id = lambda d: np.where(d.grouped_type.str.contains('rRNA'), 
                                            d.gene_id.map(rename_rRNA),
                                            d.gene_id))\
        .assign(gene_type = lambda d: np.where(d.grouped_type.str.contains('rRNA'), 
                                            d.gene_id.map(rename_rRNA),
                                            d.gene_type))\
        .assign(gene_id = lambda d: np.where(d.gene_type=="Mt_tRNA",d.gene_name, d.gene_id))\
        .assign(gene_type = lambda d: np.where((d.gene_type=="tRNA") & (d.gene_id.str.contains('^MT')),
                                                'Mt_tRNA',
                                                d.gene_type))\
        .assign(gene_type = lambda d: d.gene_type.str.replace('^tRNA','Nucleo-tRNA'))\
        .pipe(lambda d: d[~d.gene_id.str.contains('tRNA')])\
        .assign(treatment = lambda d: d.samplename.map(label_sample))\
        .assign(gene_id = lambda d: d.gene_id.str.replace('-[0-9]+-[0-9]+',''))\
        .groupby(['treatment','gene_id', 'gene_type'], as_index=False)\
        .agg({'read_count':'sum'})\
        .groupby(['treatment','gene_type'], as_index=False)\
        .agg({'read_count':'sum',
            'gene_id':'count'})\
        .query('treatment == "DNase I"')\
        .assign(value = lambda d: d.groupby('treatment').read_count.transform(lambda x: 100*x/x.sum()))\
        .assign(gene_type = lambda d: d.gene_type.str.replace('Mt_','MT-').str.replace('_',' '))\
        .assign(gene_type = lambda d: d.gene_type + ' '+ \
                                      d['value'].round(2).astype(str)+ '% (' +\
                                      d.gene_id.astype(str) +')')\
        .set_index('gene_type')\
        .assign(explode = lambda d: (100-d['value'])/100) \
        .assign(explode = lambda d: np.where(d.explode < 0.95,0, 
                                            np.where(d.explode < 0.99, 0.2,
                                                    np.where(d.explode < 0.994, 0.4, 0.8))))\
        .sort_values('value', ascending=False) \
        .assign(label = lambda d: np.where(d['value']>2,d.index, ''))
    ax.pie(small_df['value'], 
       labels=small_df.label,
       startangle = 40,
        wedgeprops=dict(width=0.5), colors = simpsons_palette())
    ax.legend(small_df.index, bbox_to_anchor=(1.2,0.9),
            frameon=False, fontsize=20)


def plot_count(ax, feature_only=True, dedup=True):
    countplot_df = read_count(feature_only=feature_only, dedup=dedup, rna_group_type = 'grouped_type')
   
    colors = rna_type_ce.transform(countplot_df.columns)
    countplot_df\
        .reindex(index=label_order)\
        .plot\
        .bar(stacked=True, 
             color = colors, 
              ax = ax,
             width = 0.8)
    xt = [xt.set_text(xt.get_text()) for xt in ax.get_xticklabels()]
    ax.set_xticklabels(ax.get_xticklabels(), 
                       rotation = 70, 
                       ha = 'left', va = 'top',
                      rotation_mode="anchor")
    ax.legend()
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], 
              bbox_to_anchor = (0.9,1), fontsize = 18,
             frameon=False)
    ax.set_xlabel('')
    ax.set_ylabel('Read pairs (%)')
    sns.despine()
    return countplot_df



def sample_wise_fraction():
    return pd.read_feather(work_path + '/cfNA/tgirt_map/Counts/all_counts/spreaded_all_counts.feather')\
        .filter(regex='group|dedup:sense')\
        .groupby('grouped_type', as_index=False)\
        .sum() \
        .pipe(pd.melt, id_vars = ['grouped_type'],
             var_name = 'samplename', value_name = 'rcount') \
        .assign(rfrac = lambda d: d.groupby('samplename').rcount.transform(lambda x: 100 *x/x.sum()))\
        .assign(label = lambda d: d.samplename.map(label_sample))\
        .pipe(lambda d: d[~pd.isnull(d.label)]) \
        .groupby(['label','grouped_type'], as_index=False)\
        .agg({'rfrac':['min','max']})
