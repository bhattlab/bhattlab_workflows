Example configs and snakemake submits for all, conda install description

# preprocessing 

Jessy

0) multiplex check - Eli

1) qc - fastqc

2) trim - trim galore: all adapters, min read 60

3) deduplication - seakit: fwd and rev separately

4) sync - eli's script

5) post qc

6) multiqc

7) separate reference reads and give ref alignment stats


# assembly

# smart memory allocation, multiplying runtime per attempts

spades - Eli

megahit - jessy

quast

# Classification and taxonomic barplots

Given one or more input datasets, this workflow performs taxonomic classification with kraken, then visualizes 
the resulting compositions in a barplot.

## First time setup

First, install [miniconda3](https://conda.io/miniconda.html).  Then, if you haven't already, set up a Snakemake 
profile for SCG by following the instructions [here](https://github.com/bhattlab/slurm).  Then:

```
conda env create -f /path/to/classification.yaml
```

## Running the workflow

```
source activate classification
snakemake --configfile /path/to/config.yaml -s /path/to/classification/Snakefile --profile scg
```
