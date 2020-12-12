#!/usr/bin/python3

import argparse
import sys

import pandas as pd

parser = argparse.ArgumentParser(
    description="Add a custom column and value to a tsv file")

parser.add_argument("--column", nargs=1, required=True,
    help="name of column to add to a .tsv file")    
parser.add_argument("--value", nargs=1, required=True,
    help="value for every row in the new column")
parser.add_argument("--input", nargs=1, required=True,
    help="the file to be used as input")
parser.add_argument("--output", nargs=1, required=True,
    help="the output file")
parser.add_argument("--position", nargs=1, type=int, 
    help="position for the new column, [0..num-columns)")
parser.add_argument("--overwrite", action="store_true",
    help="overwrite existing column values")
args = parser.parse_args()
colname = args.column[0]
colval = args.value[0]
inputFilename = args.input[0]
outputFilename = args.output[0]
overwrite = args.overwrite

input = pd.read_csv(inputFilename, delimiter='\t')

if colname in input and not overwrite:
    print("column '%s' already exists in %s\n" %(colname, inputFilename))
    sys.exit(1)

if args.position is None:
    pos=0
else:
    pos =args.position[0]

input.insert(pos, colname, colval)
input.to_csv(outputFilename, index=False, sep='\t')
