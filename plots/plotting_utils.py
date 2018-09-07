import matplotlib.pyplot as plt
from matplotlib import rcParams
import re
from sequencing_tools.viz_tools import color_encoder, okabeito_palette
import numpy as np
import seaborn as sns
from collections import defaultdict


plt.rc('axes', labelsize=15)
plt.rc('xtick', labelsize =15)
plt.rc('ytick', labelsize = 15)
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

def label_sample(x, salt = False):
    if 'HS' in x:
        return 'High salt (450mM)'
    elif 'Frag' in x:
        return 'Fragmented'
    elif re.search('[-_]sim',x):
        return 'WGS-sim'
    elif re.search('N[aA]|[Aa]lk', x):
        #return 'Alkaline hydrolysis'
        return 'NaOH'
    elif re.search('L[12]',x):
        return 'PolyA-selected'
    elif re.search('[eE]xo|ED|DE', x):
        return 'DNase I + Exo I'
    elif re.search('[aA]ll|[Uu]nt', x):
        return 'Untreated'
    elif re.search('[Qq][cC][Ff][0-9]+|[uU]nf', x):
        if salt:
            return 'Low salt (200mM)'
        else:
            return 'DNase I'
    
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
                         'NaOH','Untreated','Ladder'],
                         ['#d12604','#ff96cb',
                          '#964b06','#f2a157',
                          '#4286f4','black', 'grey']):
    label_ce.encoder[label] = color


RNA_type = ['Antisense', 'Mt', 'Other ncRNA', 'Other sncRNA', 'Protein coding',
            'Repeats', 'miRNA', 'rRNA', 'snoRNA', 'tRNA', 'Vault RNA','No features']
colors = okabeito_palette()
colors = sns.color_palette("Paired", 10)
colors.extend(['gray','black'])
rna_type_ce = color_encoder()
rna_type_ce.fit(RNA_type, colors)
rna_type_ce.encoder = {rna:col for rna, col in zip(RNA_type, colors)}


label_order = ['Untreated','NaOH', 'DNase I', 'DNase I + Exo I','WGS-sim']


figure_path = '/stor/work/Lambowitz/cdw2854/cell_Free_nucleotides/tgirt_map/figure'




import networkx as nx
from scipy.spatial import Delaunay

def spring_layout(ax, data, annotations, colors, iterations = 50, k=None):
    """
    - data: coordinates of your points [(x1,y1), (x2,y2), ..., (xn,yn)]
    - annotations: text for your annotation ['str_1', 'str_2', ..., 'str_n']
    - colors: colors for each annotation ['color1', 'color2',... 'color_n']
    - iterations: number of iterations for spring layout
    - k: optimal distance between nodes
    """
    if k is None:
        k = 1 / np.sqrt(len(data))
    G = nx.Graph()
    init_pos = {} # initial position
    x, y = [e[0] for e in data], [e[1] for e in data]
    xmin, xmax = min(x), max(x)
    ymin, ymax = min(y), max(y)
    tri = Delaunay(data)
    for i, text in enumerate(annotations):
        xy = data[i]
        G.add_node(xy)
        G.add_node(text)
        G.add_edge(xy, text)
        init_pos[xy] = xy
        init_pos[text] = xy
    for ijk in tri.simplices.copy():
        edges = zip(ijk, ijk[1:])
        for edge in edges:
            G.add_edge(annotations[edge[0]], annotations[edge[1]])
    pos = nx.spring_layout(G ,pos=init_pos, fixed=data, iterations=iterations,\
    k=k)
    if not colors:
        colors = ['black'] * len(annotation)
    for i, (name,color) in enumerate(zip(annotations, colors)):
        if name:
            xy = data[i]
            xytext = pos[name] * [xmax-xmin, ymax-ymin] + [xmin, ymin]
            ax.annotate(name, xy, 
                xycoords='data', xytext=xytext, textcoords='data', color = color,\
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3", \
                color=color))
