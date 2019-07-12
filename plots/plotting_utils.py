import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
import re
from sequencing_tools.viz_tools import color_encoder, okabeito_palette
import numpy as np
import seaborn as sns
from collections import defaultdict
import pandas as pd
plt.rc('axes', labelsize=25)
plt.rc('xtick', labelsize =25)
plt.rc('ytick', labelsize = 25)
plt.rc('font', **{'family':'sans-serif',
                  'sans-serif':'Arial'})


class cpm_total():
    def __init__(self, dedup=True):
        count_df = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/Counts/all_counts/all_counts.feather'
        dedup_variable = 'dedup' if dedup else 'all'
        self.count_df = pd.read_feather(count_df)\
                .query('dedup == "%s" & strand == "sense"' %dedup_variable)\
                .groupby('samplename', as_index=False)\
                .agg({'read_count':'sum'})

        self.sample_cpm = {}
        for i, row in self.count_df.iterrows():
            self.sample_cpm[row['samplename']] = row['read_count']

        self.sample_df = self.count_df\
            .assign(prep = lambda d: d.samplename.map(label_sample))\
            .groupby('prep', as_index=False)\
            .agg({'read_count':'sum'})

        self.prep_cpm = {}
        for i, row in self.sample_df.iterrows():
            self.prep_cpm[row['prep']] = row['read_count']
        
            



def label_sample(x, salt = False, new_label=False):
    if 'HS' in x:
        return 'High salt (450mM)'
    elif re.search('^[fF]rag|_[fF]rag',x):
        return 'Fragmented'
    elif re.search('[-_]sim',x):
        return 'WGS-sim'
    elif re.search('N[aA][0-9]+|[Aa]lk', x):
        #return 'Alkaline hydrolysis'
        if new_label:
            return "DNA (NaOH)"
        return 'NaOH'
    elif re.search('_L[0-9E]+',x):
        return 'Poly(A)-selected'
    elif re.search('[eE]xo|ED|DE', x):
        return 'DNase I + Exo I'
    elif re.search('[aA]ll|[Uu]nt', x):
        return 'Untreated'
    elif re.search('[pP]hos|3\'P', x):
        if new_label:
            return "RNA (DNase I - 3'P)"
        return "DNase I - 3'P"
    elif re.search('MPF4', x):
        return 'EV'
    elif re.search('MPF10', x):
        return 'RNP'
    elif re.search('MPCEV', x):
        return 'Crude'
    elif re.search('^GC', x):
        return 'HEK293'
    elif re.search('PPF4', x):
        return 'EV (MNase)'
    elif re.search('PPF10', x):
        return 'RNP (MNase)'
    elif re.search('PPCEV', x):
        return 'Crude (MNase)'
    elif re.search('[Qq][cC][Ff][0-9]+|[uU]nf', x):
        if salt:
            return 'Low salt (200mM)'
        elif new_label:
            return "RNA (DNase I)"
        else:
            return 'DNase I'

prep_order = ['DNase I', 'DNase I + Frag', 
              "DNase I - 3'P", 'SMART-Seq']
def label_prep(x):
    if re.search('[Uu]nfrag', x):
        return prep_order[0]
    elif re.search('[fF]ragme', x):
        return prep_order[1]
    elif re.search('[pP]hos', x):
        return prep_order[2]
    elif re.search('[Pp]oly', x):
        return prep_order[3]
    else:
        return x
    
    
def rename_sample(xs):
    sample_dict = defaultdict(int)
    out_name = []
    for x in xs:
        prep = label_sample(x)
        sample_dict[prep] += 1
        out_name.append('%s %i' %(prep, sample_dict[prep]))
    return out_name
    

label_ce = color_encoder()
label_ce.encoder = {}
for label, color in zip(['DNase I', 'DNase I + Exo I',
                         'DNase I + NaOH', 'DNase I + Exo I + NaOH',
                         'NaOH','Untreated','Ladder','Fragmented',"DNase I - 3'P",
                         'HEK293'],
                         ['#d12604','#ff96cb',
                          '#964b06','#f2a157',
                          '#4286f4','black', 'grey',
                          '#592782','#870c47','black']):
    label_ce.encoder[label] = color


RNA_type = ['Antisense', 'Mt', 'Other ncRNA', 'Other sncRNA', 'Protein coding',
            'Repeats', 'miRNA', 'rRNA', 'snoRNA', 'tRNA', 'Vault RNA','Unannotated',
            '5/5.8S rRNA', '18/28S rRNA','Mt-tRNA']
colors = okabeito_palette()
colors = sns.color_palette("Paired", 10)
colors.extend(['gray','black','#f2cf5c','#f29c07','#1797cf'])
rna_type_ce = color_encoder()
rna_type_ce.fit(RNA_type, colors)
rna_type_ce.encoder = {rna:col for rna, col in zip(RNA_type, colors)}




figure_path = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/figure'
if not os.path.isdir(figure_path):
    os.mkdir(figure_path)


def plot_upset(fig, upset_df, ylab = 'Number of full-length intron'):
    '''
    upset_df:
        pandas dataframe, contain column: sample matrix (binary), "count", 'index'

    input example:
    	index	HeLa	K562	UHRR	Plasma	count
        3	0	1	1	1	1
        4	1	0	0	0	2
        5	1	1	0	0	3
        1	0	1	0	0	11
        6	1	1	0	1	12
        7	1	1	1	1	12
        2	0	1	0	1	13
        0	0	0	0	1	30


    '''
    bar_ax = fig.add_axes([0,0.4,1,1])
    heat_ax = fig.add_axes([0,0,1,0.4])
    upset_df.plot.bar('index','count', width=0.9,ax=bar_ax)
    bar_ax.xaxis.set_visible(False)
    bar_ax.set_xlim(0, upset_df.shape[0]-0.5)
    bar_ax.legend().set_visible(False)

    ymax = round(upset_df['count'].max(), -1)
    ticks = np.linspace(0, ymax, 5)
    ticks = np.round(ticks, -1)[1:]
    ticks = np.array(ticks, dtype='int')
    bar_ax.set_yticks(ticks)
    bar_ax.set_yticklabels(ticks)
    [bar_ax.spines[s].set_visible(False) for s in ['top','right']]
    bar_ax.set_ylabel(ylab)
    bar_ax.set_xlim(-0.5,upset_df.shape[0]-0.5)


    matrix = upset_df.drop(['index','count'],axis=1)
    sample_number = matrix.shape[1] 
    heat_ax.imshow(matrix.transpose(), 
                aspect='auto', cmap = 'binary', alpha=0.4)
    heat_ax.set_yticks(range(sample_number))
    heat_ax.set_yticklabels(matrix.columns)
    heat_ax.xaxis.set_visible(False)
    heat_ax.hlines(y = np.arange(sample_number)+0.5, xmin=-0.5, xmax=upset_df.shape[0]-0.5)
    heat_ax.vlines(x = np.arange(upset_df.shape[0])+0.5, ymax = sample_number+0.5, ymin=-0.5)
    heat_ax.set_ylim(-0.5,sample_number-0.5)



#import networkx as nx
#from scipy.spatial import Delaunay
#
#def spring_layout(ax, data, annotations, colors, iterations = 50, k=None):
#    """
#    - data: coordinates of your points [(x1,y1), (x2,y2), ..., (xn,yn)]
#    - annotations: text for your annotation ['str_1', 'str_2', ..., 'str_n']
#    - colors: colors for each annotation ['color1', 'color2',... 'color_n']
#    - iterations: number of iterations for spring layout
#    - k: optimal distance between nodes
#    """
#    if k is None:
#        k = 1 / np.sqrt(len(data))
#    G = nx.Graph()
#    init_pos = {} # initial position
#    x, y = [e[0] for e in data], [e[1] for e in data]
#    xmin, xmax = min(x), max(x)
#    ymin, ymax = min(y), max(y)
#    tri = Delaunay(data)
#    for i, text in enumerate(annotations):
#        xy = data[i]
#        G.add_node(xy)
#        G.add_node(text)
#        G.add_edge(xy, text)
#        init_pos[xy] = xy
#        init_pos[text] = xy
#    for ijk in tri.simplices.copy():
#        edges = zip(ijk, ijk[1:])
#        for edge in edges:
#            G.add_edge(annotations[edge[0]], annotations[edge[1]])
#    pos = nx.spring_layout(G ,pos=init_pos, fixed=data, iterations=iterations,\
#    k=k)
#    if not colors:
#        colors = ['black'] * len(annotation)
#    for i, (name,color) in enumerate(zip(annotations, colors)):
#        if name:
#            xy = data[i]
#            xytext = pos[name] * [xmax-xmin, ymax-ymin] + [xmin, ymin]
#            ax.annotate(name, xy, 
#                xycoords='data', xytext=xytext, textcoords='data', color = color,\
#                arrowprops=dict(arrowstyle="->", connectionstyle="arc3", \
#                color=color))
