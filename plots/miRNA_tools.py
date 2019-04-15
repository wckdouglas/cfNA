import pysam
import RNA
import pandas as pd
import mygene

class miRNA_test():
    def __init__(self):
        miRBase = '/stor/work/Lambowitz/ref/hg19_ref/genes/mirBase.gff2.gz'
        ensembl = '/stor/work/Lambowitz/ref/hg19_ref/genes/miRNA.bed.gz'
        self.tabix = pysam.TabixFile(miRBase)
        self.ensembl = pysam.TabixFile(ensembl)

    def test(self, chrom, start, end):
        boolean = False
        try:
            boolean = any(self.tabix.fetch(chrom, start, end)) or any(self.ensembl.fetch(chrom,start,end))
        except:
            boolean = False

        return 'Yes' if boolean else 'No'



def fold_intron(seq):
    return RNA.fold(seq)[0]


def test_mirtron(fold):
    three_trim = fold.endswith('..')
    five_trim = fold.startswith('..')
    if three_trim and five_trim:
        return 'Both trim'
    elif three_trim:
        return '3-trim'
    elif five_trim:
        return '5-trim'
    else:
        return 'Conventional'




class gene_conversion():
    def __init__(self):
        self.mg = mygene.MyGeneInfo()
        ref_path = '/stor/work/Lambowitz/ref/hg19_ref/genes'
        self.trans = ref_path + '/trans.txt' 
        self.trans = pd.read_csv(self.trans, sep='\t',
                                 names = ['gsymbol','tid']) 
    
    def get_gid(self, genes):
        for res in self.mg.getgenes(genes, 
                                    scope='ensembl.gene', 
                                    fields=['symbol']):
            if 'notfound' in res:
                yield 'None'
            else:
                return res['symbol']

