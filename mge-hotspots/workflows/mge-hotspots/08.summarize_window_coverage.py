#!/usr/bin/python3
"""
This script finds the windows with the lowest pvalue score that are also
hotspots. These are written out as both a tsv and a bed file. In addition
a tsv file is written with all windows for all species along with the
new hotspot information which can be used to generate a plot with
genome wide insertions and hotspots, see 09.plot_genomewide_insertions.py.
"""

import argparse
import sys

import numba
import numpy as np
import pandas as pd

import paths

parser = argparse.ArgumentParser(
    description="Summarize hotspot insertions")
parser.add_argument("--species", nargs=1, required=True,
    help="species being analyzed for hotspots, a .tsv file")
parser.add_argument("--output-prefix", nargs=1, required=True,
    help="directory for outputs, windows are written to <output-prefix>/{all_top_windows_detailed.tsv,all_top_windows_summarized.bed,  all_top_windows_summarized.tsv}")

args = parser.parse_args()

outputs = paths.Outputs(args.output_prefix[0])
species =  pd.read_csv(args.species[0], delimiter='\t')

@numba.jit(nopython=True)
def group_consecutive_values(signif):
    """
    group_consecutive_values returns an array with runs of group numbers
    that are incremented whenever the signif boolean values changes.
    That is, TTFFTFF will produce 1122344.
    """
    groups = np.zeros(len(signif), dtype=np.int64)
    grp = 1
    prev = signif[0]
    groups[0] = grp
    for i in range(1, len(signif)):
        if prev != signif[i]:
            grp += 1
        groups[i] = grp
        prev = signif[i]
    return groups

def per_species_tables(row):
    """
    per_species_tables creates a table for each species that contains window
    coverage data and a grouping of each run of 'significant' vs
    'non-significant' results as per group_consecutive_values. 
    """
    hotspots = pd.read_csv(
        outputs.significance_results_filename(row.species, row.genome), delimiter='\t')
    hotspots['species'] = row.species
    hotspots['genome'] = row.genome
    hotspots['group'] = group_consecutive_values(hotspots.signif.to_numpy())
    return hotspots

def per_group_min_pvalues(table):
    """
    per_group_min_pvalues adds a column 'top_window' which is True
    for the first row of each group that contains the mininum
    pvalue found in that group.
    """

    # Find the min value in each group.
    per_group_min_pvalues = table.groupby(['group']).pvalue.idxmin()
    table['top_window'] = False
    table.loc[per_group_min_pvalues.values, 'top_window'] = True

    # Find the first entry in each group that has that min value.
    first = table.groupby(['group', 'top_window']).first()
    idx = pd.IndexSlice
    first_top_windows = first.loc[idx[:,True],:]
    table['top_window'] = False
    table.loc[first_top_windows['window_order'].values-1, 'top_window'] = True

def summarize(table):
    top_windows_combined = top_windows.loc[top_windows['signif']]
    pg = top_windows_combined.groupby(['species', 'genome', 'group']).agg(
        {   'species': 'first', 'genome': 'first', 'group': 'first',
            'contig': 'unique', 'start': 'min', 'end': 'max', 'pvalue': 'min'},
    )
    pg.loc[:,'contig'] = pg['contig'].apply(lambda l: l[0])
    return pg

# Create a single table with all species/genomes and a group number that
# identifies runs of significant window counts see group_consecutive_values).
species_tables = species.apply(lambda row: per_species_tables(row),
    axis='columns', result_type='reduce')

# Update each per-species table to identify the first window that has the
# mininum pvalue in each group.
for table in species_tables:
    per_group_min_pvalues(table)

# Combine all windows for all species into a single table.
top_windows = pd.concat(species_tables.values)
# Summarize that table, i.e. keep only hotspot information.
top_windows_summarized = summarize(top_windows)

# Write the complete table.
top_windows.to_csv(outputs.top_windows_detail(), index=False, sep='\t')

# Write the summary table.
top_windows_summarized.to_csv(outputs.top_windows_summary(), index=False, sep='\t')

# Write the summary table in bed form.
bed_file = top_windows_summarized.loc[:,['contig', 'start', 'end', 'pvalue']]
bed_file = bed_file.sort_values(by=['contig', 'start', 'end'])
bed_file.to_csv(outputs.top_windows_summary_bedfile(), header=None, index=False, sep='\t')

