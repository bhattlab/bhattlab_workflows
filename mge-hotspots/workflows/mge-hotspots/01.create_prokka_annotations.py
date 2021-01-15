#!/usr/bin/python3
"""Prokka Genome Annonations
"""

import argparse
import multiprocessing
import subprocess
import sys

import pandas as pd

import paths

parser = argparse.ArgumentParser(
    description="Generate prokka annotations")
parser.add_argument("--species", nargs=1, required=True,
    help="species to be annotated, a .tsv file")
parser.add_argument("--genomes", nargs=1, required=True,
    help="directory containing reference genomes")
parser.add_argument("--output-prefix", nargs=1, required=True,
    help="directory for outputs, annotated genomes are written to <output>/")
args = parser.parse_args()

outputs = paths.Outputs(args.output_prefix[0])
species =  pd.read_csv(args.species[0], delimiter='\t')
genome_dir = args.genomes[0]

def annotate_species(genome_dir, outdir):
    result = subprocess.run(
        ["prokka",
            "--kingdom", "Bacteria", 
            "--cpus", str(multiprocessing.cpu_count()),
            genome_dir,
            "--force", "--outdir", outdir,
        ],
        encoding='utf-8',
        capture_output=True,
    )
    if result.returncode != 0:
        print(result)
        print(result.stderr)
    result.check_returncode()

species.apply(lambda sp:
    annotate_species(
        genome_dir+"/"+sp.genome+".fna",
        outputs.prokka_genome_dir(sp.species),
    ),
    axis='columns')

