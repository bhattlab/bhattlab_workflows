#!/usr/bin/python3

import argparse
import sys

import pandas as pd

import paths

parser = argparse.ArgumentParser(
    description="Compare hotspot results")
parser.add_argument("--species", nargs=1, required=True,
    help="species hotspots to be compared, a .tsv file")    
parser.add_argument("--output-a", nargs=1, required=True,
    help="directory for outputs-a")
parser.add_argument("--output-b", nargs=1, required=True,
    help="directory for outputs-b")
args = parser.parse_args()

aout = paths.Outputs(args.output_a[0])
bout = paths.Outputs(args.output_b[0])
species =  pd.read_csv(args.species[0], delimiter='\t')


def read_results(outputs, species, genome):
    filename = outputs.significance_results_filename(species, genome)
    return pd.read_csv(filename, sep='\t')

def compare_frames(a, b, columns, species, genome, precision):
    ra = read_results(a, species, genome)
    rb = read_results(b, species, genome)
    if precision == None:
        hsa = ra[columns]
        hsb = rb[columns]
    else:
       hsa = ra[columns].round(precision)
       hsb = rb[columns].round(precision)
    return hsa.compare(hsb)

def compare_columns(species, columns, precision):
    return species.apply(lambda row:
        compare_frames(
            aout,
            bout,
            columns,
            row.species,
            row.genome,
            precision,
        ), axis='columns')

exact_columns = ['contig', 'start', 'end', 'count', 'signif',
    'rate_bg', 'rate_1kb', 'rate_10kb', 'max_rate']
species['exact'] = compare_columns(species, exact_columns, precision=None)
# allow for floating point differences.
dp8 = ['qval', 'pvalue']
species['8dp'] = compare_columns(species, dp8, 8)

def report(msg, column, species, genome):
    if row[column].empty:
        return True
    print("Failed ", msg, species, genome)
    print(row[column])
    return False

errs = []
for _, row in species.iterrows():
    errs.append(report('exact match', 'exact', row.species, row.genome))
    errs.append(report('8 decimal places', '8dp', row.species, row.genome))

if sum(errs) != len(errs):
    sys.exit(1)
else:
    print("ok")

