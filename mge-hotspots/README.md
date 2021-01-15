This directory contains a python implementation of the mobile genetic element hotspot analysis original described in:

```
Mobile genetic element insertions drive antibiotic resistance across pathogens
Matthew G. Durrant, Michelle M. Li, Ben Siranosian, Ami S. Bhatt
bioRxiv 527788; doi: https://doi.org/10.1101/527788
```

In addition it also adds prokka annotations and produces tables of
the hotspots closest to a gene as well as those that fall within
a gene.

The pipeline is implemented using `nextflow` and can be run using
`workflows/mge-hotspots/run-nextflow.sh` which takes a single
flag identifyin the parameter to use (`workflows/params`). The parameters
files simple point to the locations of the input files and where
the results should be placed. The inputs are as follows:

1. a species file defining the species to be included in the analysis.
2. a directory with reference genomes for the species being analyzed
3. the genotypes and clusters files produced by mge-finder.

The `conda-environment.yml` file lists the packages required to create
a `minconda` environment suitable for running the analysis.

# Inputs

The inputs are all tsv's, and the columns are as follows. Note that only
those marked with a * are used directly by the python scripts.

```
species tsv file:
  * species
  * genome
    species_abbrev
    genome_abbrev
    name_abbrev
    color
    genome_name
```

```
genotype tsv file:
  * species 
  * sample
    pair_id
  * contig
  * pos_3p
  * pos_5p
    seqid
  * cluster
  * group
    conf

clusters tsv:
  * species
  * cluster
  * group
  * num_unique_sites_all
    num_unique_sites_unambig
    num_unique_seqs
    mean_length
    min_length
    max_length
    num_samples
    mean_copy
    max_copy
    IAwFC
    IAwHC
    IO
    IAwoC
    IDB
    ArSC
    ArMS
    ArML
    A
    repr_seqid
    repr_seq_length
    repr_seq
```

# Outputs

## Plots

genomewide_insertion_hotspots.pdf

window_significance_plots/<species>.<genome>.pdf

## Tables
window_significance_results/<species>window_signif.tsv

closest_genes/<species>{closest_genes,intergenic_insertions}.tsv
