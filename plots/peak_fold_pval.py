from matplotlib import use as mpl_use
mpl_use('agg')
import random
from collections import Counter, defaultdict
from scipy.stats import ranksums
from scipy.special import ndtr
from peak_utils import *
import RNA
from multiprocessing import Pool


def get_pvalues(args):
    i, row = args
    simulation = 10000
    bases = list('ACTG')
    seq = row['seq']
    sense = row['is_sense']
    fold, energy = RNA.fold(seq)
    b_counts = Counter(seq.upper())
    weights = [b_counts[n] for n in bases]
    random_energy = []
    for i in range(simulation):
        random_seq = ''.join(random.choices(bases, k=len(seq), weights=weights))
        s, e = RNA.fold(random_seq)
        random_energy.append(e)
    p_val = sum(1 for e in random_energy if e <= energy)
    return p_val/float(simulation)



def main():
    project_path = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map'
    peak_path = project_path + '/bed_files/merged_bed/MACS2/annotated'
    peak_tsv = peak_path + '/unfragmented.filtered.tsv'

    p = Pool(24)
    peak_df = load_peaks(peak_tsv)  \
        .query('sample_count >= %i &  pileup > %i' %(sample_cutoff, pileup_cutoff))\
        .assign(sense_gtype = lambda d: np.where(d.sense_gtype == ".", 'Unannotated', d.sense_gtype))\
        .assign(antisense_gtype = lambda d: np.where(d.antisense_gtype == ".", 'Unannotated', d.antisense_gtype)) \
        .sort_values('pileup', ascending=False) \
        .reset_index(drop=True)\
        .assign(is_hb = lambda d: [is_hb(row) for i, row in d.iterrows()])\
        .assign(seq = lambda d: list(map(fetch_seq, d.chrom, d.start, d.end, d.strand)))\
        .assign(fold_pval = lambda d: p.map(get_pvalues, d.iterrows()))
    p.close()
    p.join()
    peak_df.to_feather(peak_path + '/fold_pval_peak.feather')


if __name__ == '__main__':
    main()
