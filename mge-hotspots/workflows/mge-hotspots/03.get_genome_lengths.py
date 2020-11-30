#!/usr/bin/python3
'''
This script creates per-genome input files for the bedtools makewindows command.
'''
import argparse
import os
import subprocess
import sys

import pandas as pd

import paths

parser = argparse.ArgumentParser(
    description='Determine the lengths of the genomes for a given set of species')
parser.add_argument("--species", nargs=1, required=True,
    help="species whose genome lengths are to be determined, a .tsv file")
parser.add_argument("--genomes", nargs=1, required=True,
    help="directory containing reference genomes")
parser.add_argument("--output-prefix", nargs=1, required=True,
    help="directory for outputs, lengths are written to <output-prefix>/genome_lengths")
args = parser.parse_args()
genome_dir = args.genomes[0]

outputs = paths.Outputs(args.output_prefix[0])
lengths_dir = outputs.lengths_dir
lengths_tsv = outputs.lengths_tsv

species = pd.read_csv(args.species[0], delimiter='\t')

genome_info = pd.DataFrame({
    'species': species['species'],
    'genome': species['genome'],
    'genome_file': species['genome'].apply(lambda g: paths.genome_filename(genome_dir, g))})

def genomeLength(genome_file, lengths_file):
    args = ["samtools", "faidx", "-o", "-", genome_file]
    result = subprocess.run(args, encoding='utf-8', capture_output=True)
    fields = result.stdout.split('\t')
    if len(fields) != 2:
        print(result)
        result.check_returncode()
    f = open(lengths_file, "w")
    f.write(result.stdout)
    f.close()
    return [fields[0], int(fields[1])]

# determine the region and length of each genome. Note that the to_list() is
# is required to turn the pd.Series into columns.
regionLengths = genome_info.apply(lambda row:
    genomeLength(
        row.genome_file,
        outputs.length_filename(row.genome),
    ),
    axis='columns',
).to_list()

genome_info = pd.concat([
    genome_info,
    pd.DataFrame(regionLengths, columns=('region', 'length'))
    ], axis='columns')


genome_info.to_csv(lengths_tsv, index=False, sep='\t')
