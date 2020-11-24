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

import paths
import qvalue

import numba
import numpy as np
import pandas as pd
import scipy.stats as stats


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

def gen_global_background_coverage_rate(insertions, counts):
    '''
    calculate the background coverage rate, ie. total insertions
    divided by the genome-length.
    '''
    total_insertions = len(insertions.index)
    # The length of a each contig is the maximum end position for its windows;
    # these are then summed to obtain the length of the entire genome.
    genome_length = counts.groupby(by='contig').max().end.sum()
    print(total_insertions, genome_length, total_insertions/genome_length)
    return float(total_insertions) / float(genome_length)

@numba.jit(nopython=True)
def local_window_numba(starts, w_start, w_end, bg_start, bg_end):
    '''
    Implement the window bounds check for gen_local_background_coverage_rate
    using numba for a large speed up over interpreted python.
    '''
    nrows = 0
    for i in range(len(starts)):
        if ( ((starts[i] < w_start) & (starts[i] > bg_start)) | 
               ((starts[i] > w_end) & (starts[i] < bg_end))):
            nrows += 1
    return nrows

def gen_local_background_coverage_rate(insertions, contig, window_start, window_end, width):
    '''
    Calculate the background coverage rate for the portion of
    the genome immediately surrounding the specified window.
    Graphically, the regions used to compute the coverage
    are shown by +++ below.

       background_start
       +++
       window_start
       ...
       window_end
       +++
       window_end
    '''
    bg_start = window_start - (width / 2)
    bg_end = window_end + (width / 2)

    # NOTE: the combination of local_background_coverage_rate and
    # this function is O(n^2) since it iterates over all of the
    # insertion data for every count window even though only
    # local data is required. Once it's all working, revisit this
    # to fold the two loops into one approximately as follows:
    # - sort both tables by start pos
    # - for every row in the counts table, examine the corresponding
    #   bounded set of rows in the insertions table.

# TESTING ONLY...
#    insertions = insertions[insertions.contig==contig]

    hits = local_window_numba(
        insertions['start'].to_numpy(),
        window_start,
        window_end,
        bg_start,
        bg_end,
    )
    return hits / width

def local_background_coverage_rate(counts, insertions, width):
    '''
    apply gen_local_background_coverate_rate to counts for the specified
    insertions and width. See gen_local_background_coverate_rate.
    '''
    return counts.apply(lambda window:
        gen_local_background_coverage_rate(
            insertions,
            window.contig,
            window.start,
            window.end,
            width), axis='columns')

def analyze_species(outputs, species, genome):
    insertions = read_mobile_insertions(outputs, species)
    counts = read_window_counts(outputs, species)

    counts['rate_bg'] = gen_global_background_coverage_rate(
        insertions, counts)

    counts['rate_1kb'] = local_background_coverage_rate(
        counts, insertions, 1000)

    counts['rate_10kb']  = local_background_coverage_rate(
        counts, insertions, 10000)

    counts['max_rate'] = counts.loc[:,('rate_bg', 'rate_1kb', 'rate_10kb')].max(axis='columns')

    counts['pvalue'] = pvalues(
        counts['count'].to_numpy(),
        counts['start'].to_numpy(),
        counts['end'].to_numpy(),
        counts['max_rate'].to_numpy(),
    )

    # NOTE: does this need to be grouped by contig?
    counts['qval'] = qvalues(counts['pvalue'])
    counts['signif'] = counts['qval'] <= 0.05
    counts['window_order'] = counts.index + 1
    return {"counts": counts, "species": species, "genome": genome}

def pvalues(count, start, end, max_rate):
    '''
    Calculate the probability that the observed count is greater than the
    average for the interval (hence interval * per base average) assuming a
    poisson distribution of counts. This is a 'greater than'
    or 'right-tail' test, hence 1-cdf. This can be used as a p-value
    with the poisson distribution as the test-statistic when
    computing q-values.
    '''
    n = len(count)
    pv = np.zeros(n)
    for i in range(n):
        pv[i] = 1 - stats.poisson.cdf(count[i]-1, (end[i]-start[i])*max_rate[i])

    return pv

def qvalues(pv):
    '''
    Estimate the q-values for the supplied p-values.
    '''
    estimator = qvalue.QValue(pv)
    return estimator.qvalue()

results = species.apply(lambda row: 
    analyze_species(
        outputs,
        row.species,
        row.genome), axis='columns')

def output_results(outputs, counts, species, genome):
    filename = outputs.significance_results_filename(species, genome)
    counts.to_csv(filename, index=False, sep='\t')

for r in results:
    output_results(outputs, r['counts'], r['species'], r['genome'])
