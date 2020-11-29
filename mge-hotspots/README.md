This directory contains a python implementation of the mobile genetic element hotspot analysis original described in:

```
Mobile genetic element insertions drive antibiotic resistance across pathogens
Matthew G. Durrant, Michelle M. Li, Ben Siranosian, Ami S. Bhatt
bioRxiv 527788; doi: https://doi.org/10.1101/527788
```

The pipeline is implemented using `nextdflow`, though its very simple
and could just as easily be run from a shell script. It takes on the order of two minutes to run. The pipeline can be run using
`workflows/mge-hotspots/run-nextflow.sh` which takes a single
flag identifyin the parametert to use (`workflows/params`). The parameters
files simple point to the locations of the input files and where
the results should be placed. The inputs are as follows:

1. a species file defining the species to be included in the analysis.
2. a directory with reference genomes for the species being analyzed
3. the genotypes and clusters files produced by mge-finder.

These are all tsv's, and the headers are as follows:

```
species tsv:
    species
    genome
    species_abbrev
    genome_abbrev
    name_abbrev
    color
    genome_name

genotype tsv:
    species
    sample
    pair_id
    contig
    pos_3p
    pos_5p
    seqid
    cluster
    group
    conf

clusters tsv:
    species
    cluster
    group
    num_unique_sites_all
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