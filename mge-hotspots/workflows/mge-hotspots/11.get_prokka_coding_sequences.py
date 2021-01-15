#!/usr/bin/python3
"""
This script extracts coding sequence information from the prokka
genome annotations to produce a bed file and tsv file with those
coding sequences for use in downstream analysis.
"""
import argparse
import csv

import pandas as pd

import paths

parser = argparse.ArgumentParser(
    description="Extract coding sequences from prokka annotations")
parser.add_argument("--species", nargs=1, required=True,
    help="species being analyzed, a .tsv file")
parser.add_argument("--output-prefix", nargs=1, required=True,
    help="directory for outputs, windows are written to <output-prefix>/{prokka_annot/<species>.coding_seqs.bed,all_prokka_genomes_name_conversion.tsv}")
args = parser.parse_args()

outputs = paths.Outputs(args.output_prefix[0])
species =  pd.read_csv(args.species[0], delimiter='\t')

class gff_file:
    def __init__(self, filename):
        self.n = filename
        self.f = open(filename, 'r')
        self.insequence = False

    def __iter__(self):
        return self

    def __next__(self):
        line = next(self.f).strip()

        # skip all lines that are not sequence-regions.
        while line.startswith('#') or len(line) == 0 or not self.insequence:
            if line.startswith('##sequence-region'):
                self.insequence = True
            elif line.startswith('##'):
                self.insequence = False
            line = next(self.f).strip()
            pass

        parts = line.split('\t')
        # ignore any non CDS lines.
        while parts[2] != "CDS":
            line = next(self.f).strip()
            parts = line.split('\t')

        comment = parts[8]
        annotations = dict(tmp.split("=") for tmp in parts[8].split(";"))
        if 'gene' not in annotations:
            annotations['gene'] = "NA"

        return {
            'contig': parts[0],
            'start': int(parts[3])-1,
            'end': int(parts[4]),
            'strand': parts[6],
            'annotations': annotations
        }

def annotations_for_species(sp, conversions_tsv):
    """
    Write the bed file for each species and append the
    annotations for this species to conversions_tsv.
    """
    files = outputs.prokka_gff_files(sp)
    gf = gff_file(files[0])
    bedfile = open(outputs.prokka_bed_file(sp), 'w')
    bedfile_tsv = csv.writer(bedfile, delimiter='\t')
    for l in gf:
        conversions_tsv.writerow([
            l['annotations']['ID'],
            l['start'],
            l['end'],
            l['annotations']['gene'],
            l['annotations']['product']])
        bedfile_tsv.writerow([l['contig'], l['start'], l['end'], ".", l['annotations']['ID'], l['strand']])
    bedfile.close()


conversions = open(outputs.prokka_name_conversions(), 'w')
conversions_tsv = csv.writer(conversions, delimiter='\t')
conversions_tsv.writerow(['id', 'start', 'end', 'gene', 'product'])

species.apply(lambda row: annotations_for_species(row.species, conversions_tsv), axis='columns')

conversions.close()

