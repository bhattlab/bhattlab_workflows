#!/usr/bin/python3
"""
This script finds the closest genes to ...
"""

"""
inputs:

all_top_windows_summarized.bed -> all_top_windows_summarized.merged.bed

insert_beds/$sp_name.insertions.mobile.bed
genome_lengths/all_genomes.genomeFile.txt
prokka_annot/$sp_name.coding_seqs.bed

"""

import argparse
import csv
import io
import subprocess
import sys

import numpy as np
import pandas as pd

import paths

parser = argparse.ArgumentParser(
    description="Extract coding sequences from prokka annotations")
parser.add_argument("--species", nargs=1, required=True,
    help="species being analyzed, a .tsv file")
parser.add_argument("--output-prefix", nargs=1, required=True,
    help="directory for outputs, windows are written to <output-prefix>/closest_genes/<species>.closest_genes.tsv")
args = parser.parse_args()

outputs = paths.Outputs(args.output_prefix[0])
species =  pd.read_csv(args.species[0], delimiter='\t')


def intersectAndMerge(bed_a, bed_b):
    """
    intersectAndMerge intersects the ranges in bedfiles a and b and then
    merges the resulting ranges wherever they overlap. The result is
    a pd.Dataframe containing the intersected ranges.
    """
    intersectProc = subprocess.Popen(
        ['bedtools', 'intersect',
            '-a', bed_a,
            '-b', bed_b,
            '-wa', 
            '-wb',
        ],
        stdout=subprocess.PIPE)
    mergeProc = subprocess.Popen(
        ['bedtools', 'merge',
            '-c', '4,6,7',
            '-o', 'max,min,max',
        ],
        stdin=intersectProc.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True)
    out, errs = mergeProc.communicate()
    intersectProc.wait()

    if len(out) == 0:
        # Catch the case where the two sets of ranges have no
        # intersections.
        return None

    if errs is not None and len(errs) > 0:
        print("merge for %s %s: errors: %s" % (a, b, errs))
        sys.exit(1)

    return pd.read_csv(io.StringIO(out), header=None, delimiter='\t')

def padWindowsAndSort(ranges, genomeLength):
    """
    padWindowsAndSort extends each range taking the length of each chromosome
    into account so that the extended ranges do not extend beyond the end of a
    chromosome. ranges is pd.Dataframe and genomeLength a file containing
    lengths as understood by bedtools. padWindowsAndSort returns a
    a pd.Dataframe.
    """
    ranges.to_csv("./tmp.bed", index=False, header=False, sep='\t')
    slopProc = subprocess.Popen(
        ['bedtools', 'slop',
            '-i', "./tmp.bed",
            '-g', genomeLength,
            '-pct', '-b', '-0.5',
        ],
        stdout=subprocess.PIPE)
    sortProc = subprocess.Popen(
        ['bedtools', 'sort'],
        stdin=slopProc.stdout,
        stdout=subprocess.PIPE,
        text=True)
    out, errs = sortProc.communicate()
    slopProc.wait()

    if errs is not None and len(errs) > 0:
        print("slop and sort for %s: errors: %s" % (genome, errs))
        sys.exit(1)

    return pd.read_csv(io.StringIO(out), header=None, delimiter='\t')

def sortBedWindows(a):
    '''
    sortBedWindows returns a pd.Dataframe containing the sorted ranges in the
    bedtools file a
    '''
    sortProc = subprocess.run(
        ['bedtools', 'sort', '-i', a],
        capture_output=True)
    sortProc.check_returncode()
    return pd.read_csv(io.BytesIO(sortProc.stdout), header=None, delimiter='\t')

def findClosest(a, b, genomeLength):
    """
    findClosest returns a pd.Dataframe containing the closest ranges
    from a and b. Both must be sorted. genomeLength is a bedtools
    file containing genome lengths. a, b and the result are pd.Dataframes.
    """
    a.to_csv("./a.bed", index=False, header=False, sep='\t')
    b.to_csv("./b.bed", index=False, header=False, sep='\t')
    closestProc = subprocess.run(
        ['bedtools', 'closest',
            '-D', 'b', 
            '-a', './a.bed',
            '-b', './b.bed',
        ],   
        capture_output=subprocess.PIPE)
    closestProc.check_returncode()

    closest = pd.read_csv(io.BytesIO(closestProc.stdout), header=None, delimiter='\t')
    '''
    Filter out any rows that failed to match, the output of
    bedtools is as shown below, with the chromosome value
    being 'none' for the case where there is no match between
    the samples.

    $ bedtools closest -a a.bed -b b1.bed b2.bed
    chr1  10  20  a1  1 - 1 chr1  5   6   b1.1  1 -

    '''
    a = closest.iloc[:,0] !='none'
    b = closest.iloc[:,7] !='none'
    closest = closest.loc[np.logical_and(a,b)]
    return closest

def intersectABC(a, b, c):
    '''
    interstABC intersects a and b, and the result of that with c. 
    a is a bed file, b, c and the result are pd.Dataframes.
    '''
    b.to_csv("./b.bed", index=False, header=False, sep='\t')
    intersectProc = subprocess.run(
        ['bedtools', 'intersect',
            '-a', a,
            '-b', './b.bed',
            '-wa', '-wb',
        ],
        capture_output=True,
    )
    intersectProc.check_returncode()
    intersected = pd.read_csv(io.BytesIO(intersectProc.stdout), header=None, delimiter='\t')
    intersected = intersected.iloc[:,[0,1,2,3,12,14]]
    intersectProc = subprocess.Popen(
        ['bedtools', 'intersect',
            '-a', c,
            '-b', '-',
            '-wa', '-wb',
            '-f', '1.0'
        ],
        stdout=subprocess.PIPE,
        stdin=subprocess.PIPE,
    )
    buf = intersected.to_csv(None, index=False, header=False, sep='\t')
    out, errs = intersectProc.communicate(input=buf.encode('utf-8'))
    if errs is not None and len(errs) > 0:
        print("interstAnnotationsAndInsertions: intersect windows: %s" % (errs))
        sys.exit(1)
  
    return pd.read_csv(io.BytesIO(out), header=None, delimiter='\t')

def get_closest_genes(outputs, species, genome):
    """
    """
    all_top_windows = outputs.top_windows_merged_bedfile()
    genome_lengths = outputs.length_filename(genome)
    mobile_insertions = outputs.insertions_filename(species, "mobile")
    prokka_windows = outputs.prokka_bed_file(species)

    mobile_insertion_windows = intersectAndMerge(
        all_top_windows, mobile_insertions)

    if mobile_insertion_windows is None:
        return None

    mobile_insertion_windows = mobile_insertion_windows.iloc[:, [0, 4 ,5, 3]]

    expanded_mobile_insertion_windows = padWindowsAndSort(mobile_insertion_windows, genome_lengths)
    sorted_prokka_windows = sortBedWindows(prokka_windows)

    closestAnnotations = findClosest(
        expanded_mobile_insertion_windows,
        sorted_prokka_windows,
        genome_lengths)

    closest = intersectABC(
        all_top_windows,
        closestAnnotations,
        outputs.windows_filename(genome),
    )

    closest.to_csv(outputs.closest_genes_file(species), index=False, header=False, sep='\t')

    intergenic = closest.loc[closest.iloc[:,8]==0]
    intergenic.to_csv(outputs.intergenic_hotspots_file(species), index=False, header=False, sep='\t')

species.apply(lambda sp: 
    get_closest_genes(outputs, sp.species, sp.genome), axis='columns')
