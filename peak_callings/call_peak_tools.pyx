from __future__ import division, print_function
from builtins import range, zip
import numpy as np
cimport numpy as np
import sys
from functools import partial
from pandas import Series
from libc.math cimport sqrt
from operator import itemgetter
import cython
from cpython cimport bool
import six
import pyBigWig as pbw

cdef:
    long WINDOW_SIZE = 50000
    long HALF_WINDOW_SIZE = int(WINDOW_SIZE/2)

@cython.boundscheck(False)
@cython.wraparound(False)
def cal_zscore(wps, control_bw, long start, long end):
    '''
    using z score algorithm as https://github.com/rthurman/hotspot/tree/master/hotspot-distr
    '''
    cdef:
        long center, window_start, window_end
        int peak_window, large_window_wps_sum
        double peak_wps
        double p, expected, sigma, z_score

    center = long((start + end) / 2)
    window_start = center - HALF_WINDOW_SIZE
    window_end = center + HALF_WINDOW_SIZE

    if window_start < 0:
        window_start = 1
        window_end = WINDOW_SIZE
    
    elif window_end > len(wps):
        window_end = len(wps)
        window_start = window_end - window_end


    peak_window = end - start 
    p = peak_window / WINDOW_SIZE

    window_wps = wps[window_start:window_end]
    large_window_wps_sum = window_wps[window_wps>0].sum()
    peak_wps = wps[start:end].max() 

    expected = large_window_wps_sum * p
    
    if expected > 0:
        sigma = sqrt(expected * (1-p))
        z_score = (peak_wps - expected) / sigma 
    else:
        z_score = 0.0
    return z_score, peak_wps


@cython.boundscheck(False)
@cython.wraparound(False)
def cal_binomial_control(chrom, wps, control_bw, long start, long end):
    '''
    If control file exist, use control distribution as background and use binomial
    distribution to look at if the peak is hotspot
    '''
    cdef:
        long center, window_start, window_end
        double control_background, control_peak, test_background, test_peak
        double p, expected, sigma, z_score

    center = long((start + end) / 2)
    window_start = center - HALF_WINDOW_SIZE
    window_end = center + HALF_WINDOW_SIZE

    if window_start < 0:
        window_start = 1
        window_end = WINDOW_SIZE
    
    elif window_end > len(wps):
        window_end = len(wps)
        window_start = window_end - window_end

    control_wps = control_bw.values(chrom, window_start, window_end, numpy=True)
    control_background = control_wps[control_wps>0].sum()
    control_peak = control_bw.values(chrom, start, end, numpy=True).max()

    test_wps = wps[window_start:window_end]
    test_background = test_wps[test_wps>0].sum()
    test_peak = wps[start:end].max() 
    
    p = control_peak / control_background
    expected = test_background * p
    
    if expected > 0:
        sigma = sqrt(expected * (1-p))
        z_score = (test_peak - expected) / sigma 
    else:
        z_score = 0.0
    return z_score, test_peak




@cython.boundscheck(False)
@cython.wraparound(False)
def peak_iterator(wps_array):
    '''
    from a wps array, find locations where wps passed through 0,
    define them as peak starts and ends
    '''
    cdef:
        long start = 0
        long end = 0
        int previous_wps = 0
        long i
        int wps

    assert wps_array.ndim == 1, 'WPS array is not 1d?'

    #wps_array[wps_array > 0] = 1
    #wps_array[wps_array != 1] = -1

    for i, wps in enumerate(wps_array):
        wps = 1 if wps > 0 else -1

        if wps > previous_wps:
            start = i
        elif wps < previous_wps:
            end = i
            if start < end and start != 0:
                yield start, end
        previous_wps = wps


@cython.boundscheck(False)
@cython.wraparound(False)
def merge_and_find_peaks(peaks, int tolerance_unprotected=10):
    '''
    from an iterator of peak starts and ends,
    check if they are located close enough on the genomic axis,
    merge them if they are
    '''
    cdef:
        long previous_start, previous_end 
        long new_start, new_end
        int i

    previous_start, previous_end = six.next(peaks)
    for i, (new_start, new_end) in enumerate(peaks):
        if (new_start - previous_end) <= tolerance_unprotected:
            previous_end = new_end

        else: 
            yield previous_start, previous_end
            previous_end = new_end
            previous_start = new_start
    yield previous_start, previous_end

@cython.boundscheck(False)
@cython.wraparound(False)
def maximum_sum_subarray(ls):
    '''
    #https://gist.github.com/alabid/3734606
    '''
    if len(ls) == 0:
       raise Exception("Array empty") # should be non-empty

    cdef:
        int runSum = ls[0]
        int maxSum = ls[0]
        int i = 0
        int start = 0
        int finish = 0
        int j = 0

    for j in range(1, len(ls)):
        if ls[j] > (runSum + ls[j]):
            runSum = ls[j]
            i = j
        else:
            runSum += ls[j]

        if runSum > maxSum:
            maxSum = runSum
            start = i
            finish = j

    return start, finish + 1


def pick_peak(above_median_starts, above_median_ends, sub_wps):
    '''
        from region that has 50 < size < 150,
        pick best peak (with maximum wps score)
    '''
    cdef:
        int s, e
        
    sub_wps = np.asarray(sub_wps)
    above_median_ends = np.asarray(above_median_ends)
    above_median_starts = np.asarray(above_median_starts)
    max_wps_array = np.array([sub_wps[s:e].max() for s, e in zip(above_median_starts, above_median_ends)])
    maximum_wps = np.where(max_wps_array == max_wps_array.max())[0]
    return above_median_starts[maximum_wps], above_median_ends[maximum_wps]

def adjust_peaks(wps, peak_start, peak_end):
    '''
    getting wps between peak start and end
    subtract the median value
    '''
    sub_wps = wps[peak_start:peak_end]
    median_sub_wps = np.median(sub_wps)
    adjusted_sub_wps = sub_wps - median_sub_wps
    return adjusted_sub_wps

def adjust_median(wpsArray, window=5):
    '''
    subtract running median for the wps array
    '''
    filtered = wpsArray - Series(wpsArray).rolling(window=window).median()
    filtered[np.isnan(filtered)] = np.median(filtered[~np.isnan(filtered)])
    return np.array(filtered, dtype=np.float32)


def filter_called_peaks(wps, out_bed):
    cdef:
        str line
        long start, end
        int padding_region = 30

    with open(out_bed) as peaks:
        for line in peaks:
            fields = line.strip().split('\t')
            start, end = map(long, itemgetter(1,2)(fields))
            wps[(start - padding_region):(end + padding_region)] = 0
    return wps

cdef:
    int Z_FILTER = 10
    int PEAK_SCORE_FILTER = 10
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef int write_short_peaks(wps, control_bigwig, out_bed, chromosome, strand, bool second_pass=False, int peak_count=0):
    '''
    using the input peak locations, pick out wps values,
    assign with z_score
    '''
    cdef:
        int peak_width_threshold = 10
        int peak_iter_count, peak_width
        long peak_start, peak_end, peak_center
        double z_score, peak_score
        str peak_line, peak_name

    peaks = peak_iterator(wps)
    peak_regions = merge_and_find_peaks(peaks, tolerance_unprotected=5) 
    if control_bigwig:
        control_bw = pbw.open(control_bigwig)
        zscore_calculator = partial(cal_binomial_control, chromosome, wps, control_bw)
    else:
        zscore_calculator = partial(cal_zscore, wps)

    # core algorithm
    for peak_iter_count, (peak_start, peak_end) in enumerate(peak_regions):
        peak_width = peak_end - peak_start
        if 300 >= peak_width >= peak_width_threshold:
            peak_name = '%s_peak%i' %(chromosome, peak_count)
            z_score, peak_score = zscore_calculator(peak_start, peak_end)
            peak_center = long((peak_end + peak_start) /2)
            if (z_score >= Z_FILTER and peak_score >= PEAK_SCORE_FILTER) or second_pass:
                peak_line = '{chrom}\t{start}\t{end}\t{name}\t' \
                            '{z_score}\t{strand}\t{peak_center}\t{peak_score}'\
                            .format(chrom=chromosome,
                                    start = peak_start,
                                    end = peak_end,
                                    name = peak_name,
                                    z_score = z_score,
                                    strand = strand,
                                    peak_center = peak_center,
                                    peak_score = peak_score)
                print(peak_line, file = out_bed)
                peak_count += 1
        if peak_iter_count % 10000 == 0:
            print('[%s] Parsed %i peaks: %i-%i' %(out_bed.name, peak_iter_count, peak_start, peak_end), file = sys.stderr)

    if control_bigwig:
        control_bw.close()
    return peak_count
