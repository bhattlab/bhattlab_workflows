#!/usr/bin/python3
'''
This script creates separate .bed for all reads as well as for those
categorised as being 'mobile'.
'''

import argparse
import os
import sys

import pandas as pd

import paths

parser = argparse.ArgumentParser(
    description='Process mobile insertions to create bed files')
parser.add_argument("--genotypes", nargs=1, required=True,
    help="genotypes to be include, a .tsv file")
parser.add_argument("--species", nargs=1, required=True,
    help="species whose insertions are to be analyzed, a .tsv file")
parser.add_argument("--clusters", nargs=1, required=True,
    help="cluster summaries to be processed, a .tsv file")
parser.add_argument("--output-prefix", nargs=1, required=True,
    help="directory for outputs, insertions are written to <output-prefix>/insertions/")
parser.add_argument("--verbose", action='store_true',
    help="trace/debugging output")
parser.add_argument("--limit", type=int, nargs=1,
    help="read at most limit rows")
args = parser.parse_args()

limit = None
if args.limit != None:
    limit = args.limit[0]

genotypes = pd.read_csv(args.genotypes[0], delimiter='\t', nrows=limit)
species = pd.read_csv(args.species[0], delimiter='\t', nrows=limit)
cluster_summaries = pd.read_csv(args.clusters[0], delimiter='\t', nrows=limit)
outputs = paths.Outputs(args.output_prefix[0])
insertions_dir = outputs.insertions_dir


if args.verbose:
    print("Output directory: ", insertions_dir, "\n")
    print("Genotypes dataframe info:")
    print(genotypes.info())
    print("\nSpecies dataframe info:")
    print(species.info())
    print("\nCluster summary dataframe info")
    print(cluster_summaries.info())
    print("\n")

columns = ['species', 'contig', 'pos_3p', 'pos_5p', 'group']

all_genotypes = genotypes[columns].drop_duplicates()

# join genotypes with filtered clusters to obtain mobile genotypes that
# are relevant
tg = genotypes
# transform 'cluster1;cluster2;cluster3' etc into a list that can be used
# with explode to obtain a row for each cluster.
tg.loc[:,'cluster'] = tg['cluster'].str.split(";")
tg = tg.explode('cluster')
filtered = cluster_summaries.loc[cluster_summaries['num_unique_sites_all'] > 3]
tg = tg.merge(filtered, on=['species', 'cluster', 'group'], how='inner')
mobile_genotypes = tg[columns].drop_duplicates()

def adjust_pos(genotypes):
    '''
    Add one to the 5p pos and format both the 5p and 3p as ints rather floats.
    '''
    genotypes.loc[:,'pos_5p'] = genotypes.loc[:,'pos_5p'].apply(lambda x: str(int(x+1)))
    genotypes.loc[:,'pos_3p'] = genotypes.loc[:,'pos_3p'].apply(lambda x: str(int(x)))

adjust_pos(all_genotypes)
adjust_pos(mobile_genotypes)

def write_species_bed_file(dir, suffix, species, genotypes):
    out = genotypes.loc[genotypes['species'] == species].drop(columns='species')
    filename = outputs.insertions_filename(sp, suffix)
    out.to_csv(filename, index=False, header=False, sep='\t')

for sp in species['species']:
    write_species_bed_file(insertions_dir, "all", sp, all_genotypes)
    write_species_bed_file(insertions_dir, "mobile", sp, mobile_genotypes)
