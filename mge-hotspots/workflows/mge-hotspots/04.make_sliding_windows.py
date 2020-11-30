#!/usr/bin/python3
'''
This scripts uses bedtools makewindows to create per-genome sliding windows.
'''
import argparse
import subprocess

import pandas as pd

import paths

# sliding windows of size 500 at intervals of 50 bp see generateWindows
# below.
windowSize = 500
windowStep = 50

parser = argparse.ArgumentParser(
    description='Build sliding windows for each genome')
parser.add_argument("--species", nargs=1, required=True,
    help="species for which sliding windows are to be generated, a .tsv file")
parser.add_argument("--output-prefix", nargs=1, required=True,
    help="directory for outputs, windows are written to <output-prefix>/sliding_windows")
args = parser.parse_args()

outputs = paths.Outputs(args.output_prefix[0])

species =  pd.read_csv(args.species[0], delimiter='\t')

def generateWindows(lengths_filename, windows_filename, size, step):
    '''
    Create a file named <dir>/<genome>.sliding_windows.bed
    containing sliding windows with the requested size and step.
    '''
    # Create a temporary file with the bedtools configuration, which
    # consists of the regions and associated length in tsv format. In our
    # case there is just one row in this config.
    outFile = open(windows_filename, "w")
    args = [
        'bedtools', 'makewindows',
        '-w', size, '-s', step,
        '-g', lengths_filename,
        ]
    result = subprocess.run(args, stdout=outFile)
    if result.returncode != 0:
        print(result)
        result.check_returncode()
    outFile.close()

species.apply(lambda row: 
    generateWindows(
        outputs.length_filename(row.genome),
        outputs.windows_filename(row.genome),
        str(windowSize),
        str(windowStep),
        ), 
        axis='columns')
