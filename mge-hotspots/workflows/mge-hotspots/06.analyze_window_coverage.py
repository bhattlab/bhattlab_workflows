#!/usr/bin/python3
"""Hotspot Analysis

This script carries out the analysis described in the section
'Insertion hotspot analysis' described in:

Mobile genetic element insertions drive antibiotic resistance across pathogens
Matthew G. Durrant, Michelle M. Li, Ben Siranosian, Ami S. Bhatt
bioRxiv 527788; doi: https://doi.org/10.1101/527788

There are some differences in the details of the analysis required by
the translation to python. In particular, the R poisson.test function
is not available, neither is biconductor's qvalue package. They are
replaced as follows:

1. poisson.test's pvalue, for a right-tail, greater-than, test boils
   down to the value returned by poisson CDF. Consequently, the scipy
   equivalent function is used.
2. the qvalue module from GDSCTools is copied locally and used directly
   rather than using all of GDSCTools since that package does not yet
   support python 3.8.

Both of the above generate the same values as R code to within about
10 decimal places which was determined by comparing the results of
the R pipeline with this new python implementation.

"""

import argparse
import subprocess
import sys

import numba
import numpy as np
import pandas as pd
import scipy.stats as stats

import paths
import qvalue

parser = argparse.ArgumentParser(
    description="Analyze genome coverage of mobile insertions to determine 'hotspots'")
parser.add_argument("--species", nargs=1, required=True,
    help="species being analyzed for hotspots, a .tsv file")
parser.add_argument("--output-prefix", nargs=1, required=True,
    help="directory for outputs, windows are written to <output-prefix>/window_significance_results")
args = parser.parse_args()

outputs = paths.Outputs(args.output_prefix[0])
species =  pd.read_csv(args.species[0], delimiter='\t')

def read_mobile_insertions(outputs, species):
    # Note that the insertion files do not have a header line, hence the
    # header=None argument.
    insertions = pd.read_csv(outputs.insertions_filename(species, "mobile"),
    header=None, delimiter='\t')
    insertions.columns = ('contig', 'start', 'end', 'cluster')
    return insertions

def read_window_counts(outputs, species):
    # Note that the count files do not have a header line, hence the
    # header=None argument.
    counts = pd.read_csv(outputs.counts_filename(species), header=None, delimiter='\t')
    counts.columns = ('contig', 'start', 'end', 'count')
    return counts

def calc_global_background_coverage_rate(counts, insertions):
    '''
    calculate the background coverage rate, ie. total insertions
    divided by the genome-length.
    '''
    total_insertions = len(insertions.index)
    # The length of a each contig is the maximum end position for its windows;
    # these are then summed to obtain the length of the entire genome.
    genome_length = counts.groupby(by='contig').max().end.sum()
    return float(total_insertions) / float(genome_length)

@numba.jit(nopython=True)
def local_window_numba(starts, w_start, w_end, bg_start, bg_end):
    '''
    Implement the window bounds check for calc_local_background_coverage_rate
    using numba for a large speed up over interpreted python.
    NOTE: this is O(n^2) because of the outer loop in
    calc_local_background_coverage_rate and the inner loop but merging
    them into a single loop is slower because the custom code to do so
    is all in python rather than the jit'ed code or in pandas.DataFrame.apply.
    '''
    nrows = 0
    for i in range(len(starts)):
        if starts[i] < bg_start:
            continue
        if ( ((starts[i] > bg_start) & (starts[i] < w_start)) | 
               ((starts[i] > w_end) & (starts[i] < bg_end))):
            nrows += 1
            continue
        if starts[i] > bg_end:
            break
    return nrows

def calc_local_background_coverage_rate(counts, sorted_starts, width):
    '''
    Calculate the local background coverage rate for every window in counts.
    This is defined as the rate for the portion of the genome immediately
    surrounding the specified window. Graphically, the regions used to compute
    the coverage are shown by +++ below.

       background_start
       +++
       window_start
       ...
       window_end
       +++
       window_end

    NOTE: the insertions must be sorted.
    '''

    # Note that the combination of this apply and iteration over all
    # insertion starts is O(n^2).
    return counts.apply(lambda window:
        local_window_numba(    
            sorted_starts,
            window.start,
            window.end,
            window.start - (width/2),
            window.end + (width/2)) / width, axis='columns')

def analyze_contig(counts, insertions):
    counts['rate_bg'] = calc_global_background_coverage_rate(counts, insertions)
    starts = insertions['start'].to_numpy()
    counts['rate_1kb'] = calc_local_background_coverage_rate(counts, starts, 1000)
    counts['rate_10kb']  = calc_local_background_coverage_rate(counts, starts, 10000)
    counts['max_rate'] = counts.loc[:,('rate_bg', 'rate_1kb', 'rate_10kb')].max(axis='columns')
    counts['pvalue'] = pvalues(counts)
    counts['qval'] = qvalues(counts['pvalue'])
    counts['signif'] = counts['qval'] <= 0.05
    counts['window_order'] = counts.index + 1
    return counts

def analyze_species(outputs, species, genome):
    insertions = read_mobile_insertions(outputs, species)
    insertions = insertions.sort_values(['contig','start', 'end', 'cluster'], ascending=True)
    counts = read_window_counts(outputs, species)

    results = pd.DataFrame()
    for contig, counts_by_contig in counts.groupby('contig'):
        analyze_contig(
            counts_by_contig,
            insertions.loc[insertions['contig']==contig])
        results = pd.concat([results, counts_by_contig])

    return {"counts": results, "species": species, "genome": genome}

def pvalues(counts):#, start, end, max_rate):
    '''
    Calculate the probability that the observed count is greater than the
    average for the interval (hence interval * per base average) assuming a
    poisson distribution of counts. This is a 'greater than'
    or 'right-tail' test, hence 1-cdf. This can be used as a p-value
    with the poisson distribution as the test-statistic when
    computing q-values.
    '''
    return  1 - stats.poisson.cdf(counts['count']-1, (counts['end']-counts['start'])*counts['max_rate'])

def qvalues(pv):
    '''
    Estimate the q-values for the supplied p-values.
    '''
    estimator = qvalue.QValue(pv)
    return estimator.qvalue()

def output_results(outputs, counts, species, genome):
    filename = outputs.significance_results_filename(species, genome)
    counts.to_csv(filename, index=False, sep='\t')

results = species.apply(lambda row: 
    analyze_species(
        outputs,
        row.species,
        row.genome), axis='columns')

for r in results:
    output_results(outputs, r['counts'], r['species'], r['genome'])
