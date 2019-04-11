import pandas as pd
from collections import defaultdict
import numpy as np
from itertools import groupby
import pysam
from operator import itemgetter
from multiprocessing import Pool


class JunctionIndex():
    '''
    Index the junction sites for each transcript
    '''
    def __init__(self, junction_bed):
        self.junctions = defaultdict(list)
        with open(junction_bed) as bed:
            for line in bed:
                fields = line.strip().split('\t')
                self.junctions[fields[0]].append( int(fields[1]) + 1)
                
    def get_junctions(self, tid):
        return np.array(self.junctions[tid]) -1
                

class Transcripts():
    '''
    index transcript length
    '''
    def __init__(self, transcriptome):
        self.transcripts = {} 
        with open(transcriptome) as Ts:
            for t in Ts:
                fields = t.split('\t')
                self.transcripts[fields[0]] = (int(fields[1]), int(fields[2]))
                
    def iter_transcript(self):
        for tid, (ts, te) in self.transcripts.items():
            yield tid, ts, te
                
    def get_transcript(self, tid):
        return self.transcripts[tid]


class Pileup():
    def __init__(self):
        self.bed= pysam.Tabixfile('/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/transcriptome/unfragmented.fwd.bed.gz')
        
    def add_pileup(self, tid, start, end):
        return sum(1 for i in  self.bed.fetch(tid, start,end))
    
    def add_sample_count(self, tid, start, end):
        sample_count = set()
        for f in self.bed.fetch(tid, start, end):
            sample_count.add(f.split('\t')[3].split(':')[0])
        return len(sample_count)


class TestPeak():
    '''
    For each peak in each transcript, try to look at how likely a splice junction lies within the peak
    '''
    def __init__(self, peak_file, junction_file, transcriptome_bed, 
                 num_sim = 1000, threads = 24, filter=False):
        self.junctions = JunctionIndex(junction_file)
        self.transcripts = Transcripts(transcriptome_bed)
        self.peaks = pd.read_csv(peak_file, 
                    usecols = [0,1,2],
                    names = ['tid','peak_start','peak_end'],
                    sep='\t') \
            .drop_duplicates() 
        if filter:
            pileup = Pileup()
            self.peaks = self.peaks\
                .assign(pileup = lambda d: list(map(pileup.add_pileup, d.tid, d.peak_start, d.peak_end)))\
                .assign(sample_count = lambda d: list(map(pileup.add_sample_count, d.tid, d.peak_start, d.peak_end))) \
                .query('sample_count >= 5 & pileup > 4')
        self.num_sim = num_sim
        self.threads = threads
    
    def __extract_peak__(self, peak):
        '''
        transcript id, peak start and peak end from bed line
        '''
        fields = peak.strip().split('\t')
        return itemgetter(0,1,2)(fields)
    
    def __test_overlap__(self, peak_loc, pos):
        '''
        check with a position is within the peak region
        '''
        for ps, pe in peak_loc:
            if ps <= pos <=pe:
                return 1
        return 0
    
    def __simulate_positions__(self, tid):
        '''
        random select position to check against splice junction
        '''
        num_pos = len(self.junctions.get_junctions(tid))
        _, total_positions = self.transcripts.get_transcript(tid)
        positions = np.arange(total_positions)
        for sim in range(self.num_sim):
            yield np.random.choice(positions, num_pos, replace=False)
        

    def __test_transcripts__(self, args):
        '''
        collect all peaks from one transcript
        and test junctions
        '''
        tid, tid_df = args
        junctions = self.junctions.get_junctions(tid)
        out_row = {}
        peak_loc = list(zip(tid_df.peak_start, tid_df.peak_end))
        out_row['obs_junc_count'] = sum(self.__test_overlap__(peak_loc, pos) for pos in junctions)
        simulated_positions = self.__simulate_positions__(tid)
        for i, rand_pos in enumerate(simulated_positions):
            out_row['sim_%i' %i] = sum(self.__test_overlap__(peak_loc, pos) for pos in rand_pos)
        out_row['peak_count'] = tid_df.shape[0]
        return pd.DataFrame().from_records(out_row, index=[tid])
    
    
    def test_junction(self):
        '''
        run simulation
        '''
        names = ['sim_%i' %i  for i in range(self.num_sim)]
        names.extend(['obs_junc_count','peak_count'])
        
        chunksize = int(self.peaks.shape[0]/self.threads)
        p = Pool(self.threads)
        rows = p.map(self.__test_transcripts__, self.peaks.groupby(['tid']), chunksize = chunksize)
        p.close()
        p.join()

        return pd.concat(rows).reset_index()


def main():
    transcriptome_bed = '/stor/work/Lambowitz/ref/hg19_ref/genes/transcriptome.bed'
    junction_bed = '/stor/work/Lambowitz/ref/hg19_ref/genes/transcriptome.junction.bed'
    peaks = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/transcriptome/macs2/unfragmented.fwd_peaks.narrowPeak'
    tp = TestPeak(peaks, junction_bed, transcriptome_bed, filter=True)
    feather = '/stor/work/Lambowitz/cdw2854/cfNA/tgirt_map/transcriptome/macs2/simulation_filtered.feather'
    tp_test_df = tp.test_junction() \
        .to_feather(feather)
    print('Written %s' %feather)


if __name__ == '__main__':
    main()
