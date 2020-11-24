#!/usr/bin/python3
import argparse
import subprocess

import paths

import pandas as pd

parser = argparse.ArgumentParser(
    description='Determine genome coverage of insertions. ie. count insertions per sliding window')
parser.add_argument("--species", nargs=1, required=True,
    help="species for which window coverage is to be determined, a .tsv file")
parser.add_argument("--output-prefix", nargs=1, required=True,
    help="directory for outputs, windows are written to <output-prefix>/window_counts")
args = parser.parse_args()

outputs = paths.Outputs(args.output_prefix[0])

species =  pd.read_csv(args.species[0], delimiter='\t')

def coverage(windows_file, insertions_file, coverage_file):
    '''
    generate coverage counts for each sliding window using
    bedtools coverage -counts
    '''
    args = [
        'bedtools', 'coverage', "-counts",
        "-a", windows_file,
        "-b", insertions_file,
        ]
    out = open(coverage_file, "w")
    subprocess.run(args, stdout=out)
    out.close()
    
species.apply(lambda row: 
    coverage(
        outputs.windows_filename(row.genome),
        outputs.insertions_filename(row.species, "mobile"),
        outputs.counts_filename(row.species),
       ), 
       axis='columns')